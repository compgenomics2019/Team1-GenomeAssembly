#! /projects/home/hpan48/anaconda3/bin/python3
"""
this is genome assembly tool for team1 group1
author: hanying

assumptions:
input files are phred33 PE files
"""
import os
import shutil
import argparse


def assemble_genomes(_tmp_dir, _assemblers):
    quast_command = "quast.py -t 16 -o {0}/quast"
    a = "  1> /dev/null".format(_tmp_dir)
    for tool in _assemblers:
        quast_command += eval("run_" + tool)(_tmp_dir)
    # quast_command += " 1> /dev/null"
    os.system(quast_command)


def run_skesa(_tmp_dir):
    os.system("skesa --fastq {0}/trimmed_1P.fastq,{0}/trimmed_2P.fastq --contigs_out {0}/contigs_500_skesa_21 --kmer 21 --min_contig 500".format(_tmp_dir))
    print("-------------skesa done-------------")
    return ""


def run_spades(_tmp_dir):
    os.system("spades.py --phred-offset 33 -1 {0}/trimmed_1P.fastq -2 {0}/trimmed_2P.fastq -s {0}/trimmed_U.fastq -o {0}/spades".format(_tmp_dir))
    print("-------------spades done-------------")
    return "{0}/spades/" % _tmp_dir


def run_trim(trimmomatic_jar, _input_files, _tmp_dir, window, threshold, headcrop, crop):
    command = "java -jar %s PE -threads 15 %s %s -baseout %s/trimmed.fastq" % (trimmomatic_jar, _input_files[0], _input_files[1], _tmp_dir)
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


def run_fastqc(_input_file, _tmp_dir):
    os.system("fastqc --extract -t 15 -o %s %s 2> /dev/null 1> /dev/null" % (_tmp_dir, _input_file))


def check_crop(_tmp_dir, _fastqc_dirs):
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


def trim_files(input_files, tmp_dir, trimmomatic_jar):
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

    if not (os.popen("grep \"Total Sequences\" %s/%s/fastqc_data.txt" % (tmp_dir, fastqc_dirs[0])).readline().split()[0] == "PASS" and
            os.popen("grep \"Total Sequences\" %s/%s/fastqc_data.txt" % (tmp_dir, fastqc_dirs[1])).readline().split()[0] == "PASS"):
        trim_condition = [window_steps[0], 20, *check_crop(tmp_dir, fastqc_dirs)]

    while trim_condition is not False:
        os.system("rm -rf {0}/trimmed_*.fastq".format(tmp_dir))
        drop_rate = run_trim(trimmomatic_jar, input_files, tmp_dir, *trim_condition)
        print("-" * 10, drop_rate)
        # file name changes after trimming
        # for i, tmp in enumerate(zip(input_files, fastqc_dirs)):
        #     file, dir = tmp
        #     shutil.rmtree(os.path.join(tmp_dir, dir))
        #     run_fastqc(file, tmp_dir)
        #     os.remove("%s/%s.html" % (tmp_dir, fastqc_dirs[i]))
        #     os.remove("%s/%s.zip" % (tmp_dir, fastqc_dirs[i]))
        #
        # is_pass = os.popen("grep \"Total Sequences\" %s/%s/fastqc_data.txt" % (tmp_dir, fastqc_dirs[0])).readline().split()[0] == "PASS" and \
        #           os.popen("grep \"Total Sequences\" %s/%s/fastqc_data.txt" % (tmp_dir, fastqc_dirs[1])).readline().split()[0] == "PASS"

        if drop_rate > 33 and trim_condition != window_steps[-1]:
            trim_condition[0] = window_steps[window_steps.index(trim_condition[0]) + 1]
        else:
            trim_condition = False
    print("-------------trim_done-------------")


def main():
    supported_assemblers = ["spades", "skesa"]
    parser = argparse.ArgumentParser(description='this is description')
    parser.add_argument('-i', required=True, metavar=("file1", "file2"), nargs=2, help='pair end')
    parser.add_argument('-t', default="tmp", help='tmp folder, if exist, will be cleared')
    parser.add_argument('-n', type=int, default=1, help='threads')
    parser.add_argument('-o', default="output", help='out folder, if exist, will be cleared')
    parser.add_argument('-a', default="spades", choices=supported_assemblers, action="append", help='assemblers to use')
    parser.add_argument('-v', action="store_true", help='verbose')
    parser.add_argument('-k', action="store_true", help='set this to keep tmp')
    parser.add_argument('--trimmomatic', default="bin/trimmomatic.jar", help='point to your trimmomatic file')
    parser.add_argument('--skip-crop', action="store_true", help='set to true to skip crop step')
    parser.add_argument('--trim-only', action="store_true", help='set to true to skip assembly')
    args = parser.parse_args()

    if os.path.exists(args.t):
        shutil.rmtree(args.t)
    os.mkdir(args.t)

    trim_files(args.i, args.t, args.trimmomatic)

    if not args.trim_only:
        assemble_genomes(args.t, args.a)

    if not args.k:
        os.removedirs(args.t)


if __name__ == "__main__":
    main()
