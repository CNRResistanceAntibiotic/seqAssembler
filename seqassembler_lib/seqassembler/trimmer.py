import argparse
import logging
import os
import subprocess


def subset(in_f, in_r, output_dir):
    # https://github.com/lh3/seqtk/
    # Subsample 10000 read pairs from two large paired FASTQ files (remember to use the same random seed to keep
    # pairing):
    # seqtk sample -s100 read1.fq 10000 > in_F.fq
    # seqtk sample -s100 read2.fq 10000 > in_R.fq
    subset_list = []
    for index, infile in enumerate([in_f, in_r]):
        if infile:
            outfile = os.path.join(output_dir, 'subset_' + str(index + 1) + '.fastq.gz')
            cmd = f'seqtk sample -s100 {infile} 10000 > {outfile}'
            logging.info(f'\nSub-setting {infile} in {outfile}:\n{cmd}\n')
            out_str = subprocess.check_output(cmd, shell=True)
            logging.info(out_str)
            subset_list.append(outfile)
    return subset_list


def sickle_call(in_f, in_r, output_dir, read_type, t, q, l, g, x, n):
    # https://github.com/najoshi/sickle
    # sickle se -t sanger -q 30 -l 40 -x -n -g -f in.fastq -o tr.fastq.gz
    # sickle pe -t illumina -l 80 -q 20 -g -f in_F.fastq -r in_R.fastq -o tr_F.fastq.gz -p tr_R.fastq.gz -s
    # tr_S.fastq.gz
    if read_type == 'se':
        trim_list = [os.path.join(output_dir, 'sk_s1_se.fastq.gz')]
        cmd = f'sickle se -t {t} -q {q} -l {l}'
        option = (' -g', ' -x', ' -n')
        for index, item in enumerate([g, x, n]):
            if item:
                cmd = cmd + option[index]
        cmd = cmd + f' -f {in_f} -o {trim_list[0]} '
        logging.info(f'\nTrimming unpaired file {in_f} with sickle:\n{cmd}\n')
        out_str = subprocess.check_output(cmd, shell=True)
        logging.info(out_str)
    else:
        trim_list = [os.path.join(output_dir, 'sk_s1_pe.fastq.gz'), os.path.join(output_dir, 'sk_s2_pe.fastq.gz'),
                     os.path.join(output_dir, 'sk_s3_up.fastq.gz')]
        cmd = f'sickle pe -t {t} -q {q} -l {l}'
        option = (' -g', ' -x', ' -n')
        for index, item in enumerate([g, x, n]):
            if item:
                cmd = cmd + option[index]
        cmd = cmd + f' -f {in_f} -r {in_r} -o {trim_list[0]} -p {trim_list[1]} -s {trim_list[2]}'
        logging.info(f'\nTrimming paired-end files {in_f} and {in_r} with sickle:\n{cmd}\n')
        out_str = subprocess.check_output(cmd, shell=True)
        logging.info(out_str)
    return trim_list


def trimmomatic_call(in_f, in_r, output_dir, read_type, trim_c, trim_l, trim_t, trim_w, trim_q, trim_m):
    # java -jar $TRIM/trimmomatic.jar PE in_f in_r s1_pe s1_se s2_pe s2_se
    # ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10
    # LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40

    # java -jar $TRIM/trimmomatic-0.30.jar SE in_f s1_se.fq.gz
    # ILLUMINACLIP:$TRIM/adapters/TruSeq3-SE:2:30:10
    # LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40

    trimmomatic_dir = '/usr/local/Trimmomatic-0.38'
    adapters_dir = os.path.join(trimmomatic_dir, 'adapters')
    if read_type == 'se':
        trim_list = [os.path.join(output_dir, 'tr_s1_se.fastq.gz')]
        cmd = f"java -jar trimmomatic.jar SE -threads 4 -phred33 {in_f} {trim_list[0]} ILLUMINACLIP:" \
              f"{os.path.join(adapters_dir, 'adapters.fasta')}:{trim_c} LEADING:{trim_l} TRAILING:{trim_t} " \
              f"SLIDINGWINDOW:{trim_w}:{trim_q} MINLEN:{trim_m}"
        logging.info(f'Trimming unpaired file {in_f} with trimmomatic:\n{cmd}\n')
        out_str = subprocess.check_output(cmd, shell=True)
        logging.info(out_str)

    else:
        trim_list = [os.path.join(output_dir, 'tr_s1_pe.fastq.gz'), os.path.join(output_dir, 'tr_s1_up.fastq.gz'),
                     os.path.join(output_dir, 'tr_s2_pe.fastq.gz'), os.path.join(output_dir, 'tr_s2_up.fastq.gz'),
                     os.path.join(output_dir, 'tr_s3_up.fastq.gz')]
        cmd = f'java -jar trimmomatic.jar PE -threads 4 -phred33 {in_f} {in_r} {trim_list[0]} {trim_list[1]}' \
              f' {trim_list[2]} {trim_list[3]} ILLUMINACLIP:{os.path.join(adapters_dir, "adapters.fa")}:{trim_c}' \
              f' LEADING:{trim_l} TRAILING:{trim_t} SLIDINGWINDOW:{trim_w}:{trim_q} MINLEN:{trim_m}'
        logging.info(f'Trimming paired-end files {in_f} and {in_r} with trimmomatic:\n{cmd}\n')
        out_str = subprocess.check_output(cmd, shell=True)
        logging.info(out_str)
        cmd = f'cat {trim_list[1]} {trim_list[3]} > {trim_list[4]}'
        subprocess.check_output(cmd, shell=True)

        # remove files
        os.remove(trim_list[1])
        os.remove(trim_list[3])

        trim_list = [os.path.join(output_dir, 'tr_s1_pe.fastq.gz'),
                     os.path.join(output_dir, 'tr_s2_pe.fastq.gz'),
                     os.path.join(output_dir, 'tr_s3_up.fastq.gz')]
    return trim_list


def pre_main(arguments):
    in_f = arguments.in_f
    in_r = arguments.in_r
    subset_size = arguments.subset_size

    output_dir = arguments.output_dir
    trimmomatic = arguments.trimmomatic
    sickle = arguments.sickle
    trim_c = arguments.trimc
    trim_l = arguments.triml
    trim_t = arguments.trimt
    trim_w = arguments.trimw
    trim_q = arguments.trimq
    trim_m = arguments.trimm
    sickle_t = arguments.sicklet
    sickle_q = arguments.sickleq
    sickle_l = arguments.sicklel
    sickle_g = arguments.sickleg
    sickle_x = arguments.sicklex
    sickle_n = arguments.sicklen
    main(in_f, in_r, output_dir, subset_size, trim_c, trim_l, trim_t, trim_w, trim_q, trim_m, sickle_t, sickle_q,
         sickle_l, sickle_g, sickle_x, sickle_n, trimmomatic, sickle)


def main(in_f, in_r, output_dir, subset_size="all", trim_c="2:30:10", trim_l=3, trim_t=3, trim_w=4, trim_q=15,
         trim_m=36, sickle_t="sanger", sickle_q=20, sickle_l=40, sickle_g=False, sickle_x=False, sickle_n=False,
         trimmomatic=False, sickle=False):
    if in_r == '':
        read_type = 'se'
        subset_size = 'all'
    else:
        read_type = 'pe'

    if output_dir == '':
        output_dir = os.path.abspath(os.path.splitext(in_f)[0])
    if not os.path.exists(output_dir):
        cmd = f'mkdir {output_dir}'
        logging.info(f'Make output directory:\n{cmd}\n')
        out_str = subprocess.check_output(cmd, shell=True)
        logging.info(out_str)

    # PARSE SUBSET SIZE FOR SEQTK
    if subset_size == 'all':
        input_list = [in_f, in_r]
    else:
        input_list = subset(in_f, in_r, output_dir)

    # PARSE TRIMMOMATIC ARGUMENTS FOR READ TRIMMING
    if trimmomatic:
        trimmomatic_call(input_list[0], input_list[1], output_dir, read_type, trim_c, trim_l, trim_t, trim_w, trim_q,
                         trim_m)

    # PARSE SICKLE ARGUMENTS FOR READ TRIMMING
    if sickle:
        sickle_call(input_list[0], input_list[1], output_dir, read_type, sickle_t, sickle_q, sickle_l, sickle_g,
                    sickle_x, sickle_n)


def version():
    return "0.0.1"


def run():
    parser = argparse.ArgumentParser(description='Trimming')
    # READS INPUT FILES
    parser.add_argument('-f', '--in_forward', dest='in_f', action="store", default='input1.fastq.gz',
                        help='Input forward read file (fastq or fastq.gz) [input1.fastq.gz]')
    parser.add_argument('-r', '--in_reverse', dest='in_r', action="store", default='input2.fastq.gz',
                        help='Input reverse read file (fastq or fastq.gz) [input2.fastq.gz]')

    # OUTPUT DIR
    parser.add_argument('-o', '--output_dir', dest='output_dir', action="store", default='', help='Output directory')

    # SUBSET SIZE WITH SEQTK
    parser.add_argument('-s', '--subset', dest='subset_size', action="store", default='all',
                        help='Size of the read subset. default:[all]')

    # TRIMMING WITH TRIMMOMATIC
    parser.add_argument('-tr', '--trimmomatic', dest='trimmomatic', action="store_true",
                        help='Trimming with trimmomatic')
    parser.add_argument('--trimc', dest='trimc', action="store", default='2:30:10',
                        help='sickle option trailing [2:30:10]')
    parser.add_argument('--triml', dest='triml', action="store", default='3', help='sickle option leading [3]')
    parser.add_argument('--trimt', dest='trimt', action="store", default='3', help='sickle option trailing [3]')
    parser.add_argument('--trimw', dest='trimw', action="store", default='4', help='sickle option windows_size [4]')
    parser.add_argument('--trimq', dest='trimq', action="store", default='15', help='sickle option quality [15]')
    parser.add_argument('--trimm', dest='trimm', action="store", default='36', help='sickle option minlen [36]')

    # TRIMMING WITH SICKLE
    parser.add_argument('-sk', '--sickle', dest='sickle', action="store_true", help='Trimming with sickle')
    parser.add_argument('--sicklet', dest='sicklet', action="store", default='sanger',
                        help='sickle option type [sanger]')
    parser.add_argument('--sickleq', dest='sickleq', action="store", default='20', help='sickle option quality [20]')
    parser.add_argument('--sicklel', dest='sicklel', action="store", default='40', help='sickle option lenght [40]')
    parser.add_argument('--sickleg', dest='sickleg', action="store_true", help='sickle option g')
    parser.add_argument('--sicklen', dest='sicklen', action="store_true", help='sickle option n')
    parser.add_argument('--sicklex', dest='sicklex', action="store_true", help='sickle option x')
    # VERSION
    parser.add_argument('-V', '--version', action='version', version='rgi-' + version(), help="Prints version number")
    return parser.parse_args()


if __name__ == "__main__":
    args = run()
    pre_main(args)
