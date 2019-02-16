#!/usr/bin/perl
# Aroon Chande
# BIOL-7210 -- Comp Genomics 2016
# Genome Assembly pipeline
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Term::Report;
my $prog = basename($0);
if (@ARGV < 1){print_usage();exit 1;}
my($help,$inDir,$outDir,$i,$base,$out,$assemblyStatus,$currentStatus,$r1,$r2,$R,$ref,@steps,$link,$quast,$contigs);
GetOptions ('h' => \$help, 'o=s' => \$outDir, 'in=s' => \$inDir, 'R=s' => \$ref, 'steps=s{1,}' => \@steps, 'contigs=s' => \$contigs);
if (grep(/,/, @steps)){@steps = split(/,/,join(',',@steps));}
die print_usage() if (defined $help);
die print_usage() unless ((defined $outDir) && (defined $inDir));
# Multiplier for main progress bar
my $stepMult = scalar @steps;
$stepMult++ if (grep(/abyss/i, @steps));
# Remove generalization and make assumptions about naming scheme
my @infiles = glob ( "$inDir/*_1.fastq.gz" );
my $max = @infiles;
my $amax = $max * $stepMult;
# Initialize pipeline and submodule progress bars
initialize_bar();
$assemblyStatus->start;
$currentStatus -> start;
system(`mkdir -p $outDir/contigs`);
#SPAdes
if (grep(/spades/i, @steps)){
	update_bar("SPAdes progress:","1");
	foreach (@infiles){
		get_files($_);
		$currentStatus->subText("Running SPAdes on $out");
		system(`mkdir -p $outDir/spades/$out`);
		system(`spades.py -1 $r1 -2 $r2 -o $outDir/spades/$out -k 99,113,127 --only-assembler &>>$outDir/status.log`);
		$link = join(".","spades",$out,"fa");
		system(`cp $outDir/spades/$out/contigs.fasta $outDir/contigs/$link`);
		system(`rm -rf $outDir/spades/$out/corrected`);
		sleep 1;
		$assemblyStatus->update();
		$currentStatus->update();
	}
}
# Velvet

if (grep(/velvet/i, @steps)){
	update_bar("Velvet progress:","1");
	foreach (@infiles){
		get_files($_);
		$currentStatus->subText("Running Velvet on $out");
		system(`mkdir -p $outDir/velvet/`);
		system(`VelvetOptimiser.pl -d $outDir/velvet/$out/ -s 97 -e 127 -x 10 -f '-fastq.gz -shortPaired -separate $r1 $r2' -t 6 --optFuncKmer 'n50' &>>$outDir/status.log`);
		$link = join(".","velvet",$out,"fa");
		system(`cp $outDir/velvet/$out/contigs.fa $outDir/contigs/$link`);
		system(`rm -rf $outDir/velvet/$out/Sequences`);
		sleep 1;
		$assemblyStatus->update();
		$currentStatus->update();
	}
}
# ABySS
if (grep(/abyss/i, @steps)){
 	update_bar("ABySS progress:","2");
	foreach (@infiles){
		for my $kmer ("121"){
			get_files($_);
			$currentStatus->subText("Running ABySS on $out with k of $kmer");
			system(`mkdir -p $outDir/abyss/k$kmer`);
			system("abyss-pe -C $outDir/abyss/k$kmer k=$kmer name=$out in='$r1 $r2' j=16 &>>$outDir/status.log ");
			$link = join(".","abyss",$out,"fa");
			my $contig = join("-",$out,"contigs.fa");
			system(`cp $outDir/abyss/k$kmer/$contig $outDir/contigs/$link`);
			sleep 1;
			$currentStatus->update();
			$assemblyStatus->update();
		}
	}
}
# Smalt
if (grep(/smalt/i, @steps)){
	update_bar("SMALT progress:","3");
	system(`mkdir -p $outDir/smalt/ && cp $ref $outDir/smalt/reference.fna &>$outDir/status.log`);
	system(`smalt index -k 13 $outDir/smalt/reference  $outDir/smalt/reference.fna &>>$outDir/status.log `);
	foreach (@infiles){
		get_files($_);
		system(`mkdir -p $outDir/smalt/$out`);
		$currentStatus->subText("Mapping $out");
		system(`smalt map -F fastq -f sam -i 1000 -n 6 -o $outDir/smalt/$out/map.sam $outDir/smalt/reference $r1 $r2 &>>$outDir/status.log `);
		$currentStatus->update();
		$currentStatus->subText("Converting and sorting $out");
		system(`sambamba view -t 6 -f bam -S $outDir/smalt/$out/map.sam --output-filename $outDir/smalt/$out/map.bam &>>$outDir/status.log`);
		system(`sambamba sort -t 6 -o $outDir/smalt/$out/sorted.bam $outDir/smalt/$out/map.bam &>>$outDir/status.log`);
		system(`sambamba index -t 6 $outDir/smalt/$out/sorted.bam $outDir/smalt/sorted.index &>>$outDir/status.log`);
		$currentStatus->update();
		$currentStatus->subText("Extracting assembled reads for $out");
		system(`sambamba mpileup -t 6 $outDir/smalt/$out/sorted.bam --output-filename $outDir/smalt/$out/assembly.bcf --samtools -f $outDir/smalt/reference.fna -gu  --bcftools call -c -O b  &>>$outDir/status.log`);
		system(`bcftools view -O v $outDir/smalt/$out/assembly.bcf 2>>$outDir/status.log| vcfutils.pl vcf2fq > $outDir/smalt/$out/assembly.fq 2>>$outDir/status.log`);
		$link = join(".","smalt",$out,"fa");
		system(`seqret -sequence $outDir/smalt/$out/assembly.fq -outseq $outDir/smalt/$out/$link && cp $outDir/smalt/$out/$link $outDir/contigs/$link`);
		system(`rm $outDir/smalt/$out/map.sam $outDir/smalt/$out/map.bam &>/dev/null`);
		$currentStatus->update();
		$assemblyStatus->update();
	}
}


if (grep(/meta/i, @steps)){
	update_bar("metassembler progress:","1");
	system(`mkdir -p $outDir/meta/ 2>$outDir/status.log`);
	foreach (@infiles){
	my $confFile =join(".",$out,"config");
	open OUT, ">$outDir/$out/$confFile" or die;
	my $spades = join('.',"spades",$out,"fa");
	my $velvet = join('.',"velvet",$out,"fa");
	my $abyss = join('.',"abyss","115",$out,"fa");
	my $meta = join(".","meta",$out,"fa");
	print OUT "[global]\nbowtie2_threads=12\nbowtie2_read1=$r1\nbowtie2_read2=$r2\nbowtie2_maxins=3000\nbowtie2_minins=1000\nbowtie2_threads=6\ngenomeLength=1825000\nmateAn_A=1300\nmateAn_B=2300\n[1]\nfasta=$contigs/$spades\nID=Spades\n[2]\nfasta=$contigs/$abyss\nID=Abyss\n[3]\nfasta=$contigs/$velvet\nID=Velvet\n";
	close OUT;
	system("metassemble --conf $outDir/$out/$confFile --outd $outDir/$out 2>>$outDir/status.log");
	system("cp $outDir/meta/$out/Metassembly/QVelvet.Abyss.Spades/M1/QVelvet.Abyss.Spades.fasta $outDir/contigs/$meta");
	}
}


exit 0;

#######
sub initialize_bar(){
	my $currentReport = Term::Report->new(
   		# startRow > Make submodule progress bar start below pipeline bar
   		startRow => 4,
   		numFormat => 1,
   		statusBar => [
			scale => 50,
			label => "SPAdes progress:",
			# Disabled for now
			#showTime => 1,
			startRow => 4,
			subText => 'Running...',
			subTextAlign => 'center'
   		],
	);
	my $assemblyReport = Term::Report->new(
		startRow => 1,
		numFormat => 1,
   		statusBar => [
   			scale => 100,
    		startRow => 1,
    		label => 'Pipeline progress: ',
    		subText => 'Get some coffee',
    		subTextAlign => 'center',
    		# Need to add report generation functions
   			],
	);
	$assemblyStatus = $assemblyReport->{statusBar};
	$assemblyStatus->setItems($amax);
	$currentStatus = $currentReport->{statusBar};
	$currentStatus->setItems($max);
	$currentStatus->label("Starting pipeline:");
}
sub update_bar(){
	my ($label,$multiplier) = @_;
	my $realmax = $max * $multiplier;
	$currentStatus->reset({
    	start=>0,
        setItems => $realmax,
		label => $label,
	});
}
sub get_files(){
	$base = shift;
	$base  =~ s/_1\.fastq\.gz//g;
	($out) = $base =~ m/(SRR\d*)/;
	$r1 = join('_',$base,"1.fastq.gz");
	$r2 = join('_',$base,"2.fastq.gz");
}

sub print_usage{
    warn <<"EOF";

USAGE
  $prog -in <indir> -o <outdir> -R <reffile> --steps <steps,to,run>

DESCRIPTION
	Spades pipeline

OPTIONS
	-in	dir		Directory with fq.gz
	-o	dir		output folder
	-R	file	Reference genome file
	--steps	list	Comma separated list of steps to run
					Valid steps: abyss, velvet, spades
EXAMPLES
  $prog -in ./scratch/reads -o ./scratch/assemblies --steps velvet,abyss,spades -R ./reference/genomic.fna.gz
  $prog -h

EXIT STATUS
  0     Successful completion
  >0    An error occurred

EOF
}
