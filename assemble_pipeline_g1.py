#! /projects/home/hpan48/anaconda3/bin/python3
"""
this is genome assembly tool for team1 group1
author: hanying

assumptions:
input files are phred33 PE files
"""
import os
import pandas as pd
import numpy as np
import shutil
import argparse


def assemble_genomes(_tmp_dir, _assemblers, _threads):
    """
    run different assemblers and choose the best result
    :param _tmp_dir: tmp directory
    :param _assemblers: assemblers given by user
    :return: None
    """
    quast_command = "quast.py -t %d -o %s/quast" % (_threads, _tmp_dir)
    print(_assemblers)
    for tool in _assemblers:
        quast_command += eval("run_" + tool)(_tmp_dir)
    # quast_command += " 1> /dev/null"
    os.system(quast_command)
    quast_result = "%s/quast/report.tsv" % _tmp_dir
    result = pd.read_table(quast_result)
    result.loc["score"] = np.log(result.loc["Total length (>= 0 bp)"] * result.loc["N50"] / result.loc["# contigs"])
    best = result.loc["score"].idxmax()
    print("-" * 20 + "quast finished" + "-" * 20)
    print(result)
    os.system("mv %s/%s genome.fasta" % (_tmp_dir, best))
    print("best assembly is %s" % best)

def run_abyss(_tmp_dir):
    """
    run abyss on trimmed file
    :param _tmp_dir: tmp directory
    :return: output contigs file name(s)
    """
    out_files = []
    abyss_command = "abyss_pe in={0}/trimmed_1.fastq {0}/trimmed_2.fastq".format(_tmp_dir)
    # kmers = [21, 77, 99]
    kmers = [77]
    for k in kmers:
        os.system(abyss_command + " name=abyss-%d" % k)
        out_files.append("abyss-%d-contigs.fa" % k)
    return " " + " ".join(out_files)


def run_skesa(_tmp_dir):
    """
    run skesa on trimmed file
    :param _tmp_dir: tmp directory
    :return: output contigs file name
    """
    skesa_cmd = "skesa --fastq {0}/trimmed_1P.fastq,{0}/trimmed_2P.fastq --contigs_out {0}/skesa_contigs_21_500 --kmer 21 --min_contig 500".format(_tmp_dir)
    # print(skesa_cmd)
    os.system(skesa_cmd)
    print("-" * 20 + "skesa finished" + "-" * 20)
    return " %s/skesa_contigs_21_500" % _tmp_dir


def run_spades(_tmp_dir):
    """
    run spades on trimmed file
    :param _tmp_dir: tmp directory
    :return: output contigs file name
    """
    spades_cmd = "spades.py --only-assembler -k 21 --phred-offset 33 -1 {0}/trimmed_1P.fastq -2 {0}/trimmed_2P.fastq -s {0}/trimmed_U.fastq -o {0}/spades".format(_tmp_dir)
    # spades_cmd = "spades.py --phred-offset 33 -1 {0}/trimmed_1P.fastq -2 {0}/trimmed_2P.fastq -s {0}/trimmed_U.fastq -o {0}/spades".format(_tmp_dir)
    print(spades_cmd)
    os.system(spades_cmd)
    print("-" * 20 + "spades finished" + "-" * 20)
    os.system("mv {0}/spades/contigs.fasta {0}/spades/spades_contigs.fasta".format(_tmp_dir))
    return "%s/spades/spades_contigs.fasta" % _tmp_dir


def run_trim(trimmomatic_jar, _input_files, _tmp_dir, _threads, _skip_crop, window, threshold, headcrop, crop):
    """
    run trimmomatic on fastq file
    :param trimmomatic_jar: trimmomatic jar file
    :param _input_files: input fastq file
    :param _tmp_dir: tmp directory
    :param window: windown size of sliding window
    :param threshold: threshold in each window
    :param headcrop: hard crop from beginning
    :param crop: hard crop from end
    :return: drop rate of this trimming
    """
    command = "java -jar %s PE -threads %d %s %s -baseout %s/trimmed.fastq" % (trimmomatic_jar, _threads,_input_files[0], _input_files[1], _tmp_dir)
    if not _skip_crop:
        if headcrop:
            command += " HEADCROP:%d" % headcrop
        if crop:
            command += " CROP:%d" % crop
    command += " SLIDINGWINDOW:%d:%d MINLEN:100 2>&1" % (window, threshold)
    print(command)
    trim_summary = os.popen(command).readlines()
    # print(trim_summary)
    drop_rate = trim_summary[-2].strip().split()[-1][1:-2]
    os.system("cat {0}/trimmed_1U.fastq {0}/trimmed_2U.fastq > {0}/trimmed_U.fastq".format(_tmp_dir))
    return float(drop_rate)


def run_fake_trim(trimmomatic_jar, _input_files, _tmp_dir, _threads):
    """
    generate fastq file for assembler, do not change file content
    :param trimmomatic_jar: trimmomatic jar file
    :param _input_files: input fastq file
    :param _tmp_dir: tmp directory
    :return: None
    """
    command = "java -jar %s PE -threads %d %s %s -baseout %s/trimmed.fastq MINLEN:100 2>&1" % (trimmomatic_jar, _threads,_input_files[0], _input_files[1], _tmp_dir)
    os.popen(command)
    os.system("rm -rf {0}/trimmed_*U.fastq".format(_tmp_dir))


def run_fastqc(_input_file, _tmp_dir):
    """
    run fastqc for each file
    :param _input_file: input fastq file
    :param _tmp_dir: tmp directory
    :return: None
    """
    os.system("fastqc --extract -t 15 -o %s %s 2> /dev/null 1> /dev/null" % (_tmp_dir, _input_file))


def check_crop(_tmp_dir, _fastqc_dirs):
    """
    check if hard crop is needed
    :param _tmp_dir: tmp directory
    :param _fastqc_dirs: location of fastqc report
    :return: [head_crop, crop]
    """
    crops = [0, 0]
    for i, dir in enumerate(_fastqc_dirs):
        qualities = []
        positions = []
        data_file_name = "%s/%s/fastqc_data.txt" % (_tmp_dir, dir)
        with open(data_file_name, "r") as data:
            data_recorded = False
            for line in data:
                if line.startswith("#"):
                    continue
                elif "Per base sequence quality" in line:
                    data_recorded = True
                elif data_recorded and line.startswith(">>"):
                    break
                elif data_recorded:
                    position, quality, *tmp = line.split()
                    positions.append(position)
                    qualities.append(float(quality))
                continue
        i = 0
        while qualities[i] < 20:
            i += 1
        crops[0] = max(int(positions[i].split("-")[-1]), crops[0])
        i = len(qualities) - 1
        while qualities[i] < 20:
            i -= 1
        crops[1] = min(int(positions[i].split("-")[0]), crops[1])
    return crops


def trim_files(input_files, tmp_dir, trimmomatic_jar, threads, skip_trim, skip_crop):
    """
    Trim input files.
    Trim is done on those not passing fastqc per seq quality. sliding window is used, and window increases each step until drop rate is less that 33%.
    For those passes fastqc, a minlen:100 trimming is done for format consistency.
    :param input_files: input fastq file
    :param tmp_dir: tmp directory
    :param trimmomatic_jar: trimmomatic jar file
    :return: None
    """

    if skip_trim:
        run_fake_trim(trimmomatic_jar, input_files, tmp_dir)
        return

    window_steps = [4, 8, 12, 20, 35, 50, 70, 100]
    fastqc_dirs = ["", ""]
    trim_condition = False

    # get fastqc for raw input
    for i, file in enumerate(input_files):
        if file.endswith(".fq.gz"):
            fastqc_dirs[i] = os.path.split(file)[-1].rstrip(".fq.gz") + "_fastqc"
        elif file.endswith(".fastq.gz"):
            fastqc_dirs[i] = os.path.split(file)[-1].rstrip(".fastq.gz") + "_fastqc"
        else:
            print("irregular file name provided, check your input")
            exit(1)
        run_fastqc(file, tmp_dir)
        os.remove("%s/%s.html" % (tmp_dir, fastqc_dirs[i]))
        os.remove("%s/%s.zip" % (tmp_dir, fastqc_dirs[i]))
    print("-" * 20 + "fastqc finished" + "-" * 20)

    if os.popen("grep \"Total Sequences\" %s/%s/fastqc_data.txt" % (tmp_dir, fastqc_dirs[0])).readline().split()[0] == "PASS" and \
            os.popen("grep \"Total Sequences\" %s/%s/fastqc_data.txt" % (tmp_dir, fastqc_dirs[1])).readline().split()[0] == "PASS":
        run_fake_trim(trimmomatic_jar, input_files, tmp_dir, threads)
        return

    trim_condition = [window_steps[0], 20, *check_crop(tmp_dir, fastqc_dirs)]
    while trim_condition is not False:
        os.system("rm -rf {0}/trimmed_*.fastq".format(tmp_dir))
        drop_rate = run_trim(trimmomatic_jar, input_files, tmp_dir, threads, skip_crop, *trim_condition)
        print("-" * 10, drop_rate)
        if drop_rate > 33 and trim_condition != window_steps[-1]:
            trim_condition[0] = window_steps[window_steps.index(trim_condition[0]) + 1]
        else:
            trim_condition = False
        print("-" * 20 + "trim finished" + "-" * 20)


def main():
    supported_assemblers = ["spades", "skesa", "abyss"]
    description="""
    This is the genome assembly pipeline of team1 group1.
    For detail please see github(https://github.gatech.edu/compgenomics2019/Team1-GenomeAssembly)
    """

    usage = """
    assemble_pipeline_g1.py -t /path/to/tmp -i /path/to/file1 /path/to/file2 -a spades -a skesa -a abyss
    """

    parser = argparse.ArgumentParser(description=description, usage=usage)
    parser.add_argument('--trimmomatic', default="bin/trimmomatic.jar", metavar="trimmomatic_jar_file", help='provide trimmomatic file if you want to use your own')
    parser.add_argument('-a', required=True, choices=supported_assemblers, action="append", help='assemblers to use')
    parser.add_argument('-i', required=True, metavar=("file1", "file2"), nargs=2, help='pair end input files. MUST be gzipped fastq file')
    parser.add_argument('-t', default="tmp", metavar="tmp folder", help='tmp folder, if exist, will be cleared')
    parser.add_argument('-n', type=int, default=1, metavar="num_threads", help='number of threads to use. default: 1')
    parser.add_argument('-k', action="store_true", help='set this flag to keep tmp. default:False')
    parser.add_argument('--skip-crop', action="store_true", help='set to true to skip crop step')
    parser.add_argument('--trim-only', action="store_true", help='set to true to skip assembly')
    parser.add_argument('--assemble-only', action="store_true", help='set to true to skip trim')
    args = parser.parse_args()

    if os.path.exists(args.t):
        shutil.rmtree(args.t)
    os.mkdir(args.t)
    print("-"* 20 + "%s cleared. Now start." % args.t +"-"* 20)
    trim_files(args.i, args.t, args.trimmomatic, args.n, args.assemble_only, args.skip_crop)

    if not args.trim_only:
        assemble_genomes(args.t, args.a, args.n)

    if not args.k:
        shutil.rmtree(args.t)
    print("-" * 20 + "all finished" + "-" * 20)
    print("assembled genome is in genome.fastq")


if __name__ == "__main__":
    main()
