# Computational Genomics Genome Assembly Pipeline

This pipeline is designed by team1-group1, to assemble genomes from PE Illumina reads using a number of 
assembly programs. This pipeline is used to generate final result.

## Requirements

This pipeline is based on [python3](https://www.python.org/download/releases/3.0/) and requires pandas and numpy.

This pipeline also requires these software are installed and can be found in PATH(depends on which one you want to use).

[ABYSS](https://github.com/bcgsc/abyss)(optional)

[SPAdes](https://github.com/ablab/spades)(optional)

[SEKSA](https://github.com/ncbi/SKESA)(optional)

[QUAST](https://github.com/ablab/quast)

## GETTING STARTED

#### Quick start:
```
git clone https://github.gatech.edu/compgenomics2019/Team1-GenomeAssembly.git
cd Team1-GenomeAssembly
python3 assemble_pipeline_g1.py -k -i dataset/f1.fq.gz dataset/f2.fq.gz -a skesa -a spades
```
#### Other arguments:

`--trimmomatic` specify trimmomatic jar file here if you do not want to use the one in bin/

`-a`    declare which assembler to use, if given multiple, then the best is used as output

`-i file1 file2`        pair end input files.

`-t`         tmp folder. please be careful because tmp will be cleared when pipeline starts.

`-n`        number of threads to use in each step

`-k`                    if you do not set this, tmp is cleared when pipeline finishes

`--skip-crop`           trimmomatic will not do HEADCROP and CROP if this flag is set

`--trim-only`           assemble steps are skip if this is set, if you want to keep trim result please set `-k`

`--assemble-only`       trim steps are skip if this is set

- Input file must be gzipped fastq file, fastq file can be fetched 
from NCBI using [SRAToolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/).

- This pipeline is designed to process one sample, if you have many samples, please run this pipeline multiple times. 
Please specify different tmp folder for each sample, otherwise tmp files could mess up.

## Contact

Developed by Hanying Pan(hpan48@gatech.edu).