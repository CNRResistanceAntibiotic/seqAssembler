import os
import argparse
import shutil
from subprocess import Popen, PIPE, STDOUT


def alignment_bwa(fasta_file, fastq_file1, fastq_file2, threads=8, force=False):
    sam_file = os.path.splitext(fasta_file)[0] + '.sam'
    bam_file = os.path.splitext(fasta_file)[0] + '.bam'
    bwa_index_dir = os.path.join(os.path.dirname(bam_file), "bwa_index")

    if not os.path.exists(bwa_index_dir):
        os.mkdir(bwa_index_dir)
    else:
        print("\nIndex BWA already exist!\n")
        print(f"At {bwa_index_dir}")

    index_fasta_file = os.path.join(bwa_index_dir, os.path.basename(fasta_file))
    shutil.copy(fasta_file, index_fasta_file)

    if not os.path.exists(bam_file) or force:
        # index reference
        cmd = f"bwa index {index_fasta_file}"
        process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT).stdout.read()

        # make log
        filename_log = "logBWA_index.txt"
        header = f"Command line executed: {cmd}\n\n\n{process.decode('utf-8')}"
        log_process_output(header, bwa_index_dir, filename_log)

        # remove fasta used for index
        os.remove(index_fasta_file)

        # alignment
        cmd = f"bwa mem -t {threads} {index_fasta_file} {fastq_file1} {fastq_file2} > {sam_file}"
        process = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT).stdout.read()

        # make log
        filename_log = "logBWA_MEM.txt"
        header = f"Command line executed: {cmd}\n\n\n{process.decode('utf-8')}"
        log_process_output(header, os.path.dirname(bam_file), filename_log)
    else:
        print(f'\nAlignment file {bam_file} already done\n')
    return sam_file


def convert_sam_to_bam(sam_file, force=False):
    bam_file = os.path.splitext(sam_file)[0] + '.bam'
    if not os.path.exists(bam_file) or force:
        cmd = f'samtools view -h -b -S {sam_file} > {bam_file}'
        print(cmd)
        os.system(cmd)

        # remove sam file
        os.remove(sam_file)

    return bam_file


def sort_bam_file(bam_file):
    out_file = os.path.splitext(bam_file)[0] + '_sort.bam'
    cmd = f"samtools sort -m 1000000000 {bam_file} -T {os.path.dirname(bam_file)} > {out_file};mv {out_file} {bam_file}"
    print(cmd)
    os.system(cmd)
    cmd = f"samtools flagstat {bam_file} > {os.path.splitext(bam_file)[0] + '_bamstat.txt'}"
    print(cmd)
    os.system(cmd)
    return bam_file


def index_bam_file(bam_file):
    cmd = f"samtools index {bam_file}"
    print(cmd)
    os.system(cmd)


def split_unmapped_mapped_reads(bam_file, force):
    unmapped_fastq_file = os.path.splitext(bam_file)[0] + '_unmapped.fastq.gz'
    unmapped_bam_file = os.path.join(os.path.dirname(bam_file), 'tmp_unmapped.bam')
    if not os.path.exists(unmapped_fastq_file) or force:
        # process BAM of unmapped read
        cmd = f"samtools view -b -f 4 {bam_file} > {unmapped_bam_file}"
        print(cmd)
        os.system(cmd)

        # process FASTQ of unmapped read
        cmd = f"samtools fastq {unmapped_bam_file} > {unmapped_fastq_file}"
        print(cmd)
        os.system(cmd)

        # remove unmapped reads BAM file
        os.remove(unmapped_bam_file)

        # process BAM of mapped reads
        out_file = os.path.splitext(bam_file)[0] + '_droped.bam'
        cmd = f"samtools view -b -F 4 {bam_file} > {out_file}"
        print(cmd)
        os.system(cmd)

        # move BAM file
        shutil.move(out_file, bam_file)
    else:
        print(f'\nUnmapped reads fastq file {unmapped_fastq_file} already done\n')
    return bam_file, unmapped_fastq_file


def log_process_output(output, work_dir_path, filename_log):
    try:
        with open(f"{work_dir_path}/{filename_log}", 'w') as log_file:
            log_file.write(output)
    except IOError as e:
        return e


def pre_main(args):
    fasta_file = args.fas_file
    fastq_file1 = args.fastq_file1
    fastq_file2 = args.fastq_file2
    force = args.Force
    threads = int(args.threads)
    main(fasta_file, fastq_file1, fastq_file2, force, threads)


def main(fasta_file, fastq_file1, fastq_file2, force=False, threads=8):

    print("FASTA TO BAM arguments:\n")
    print(f"\t - Fasta File = {fasta_file}")
    print(f"\t - Fastq File 1 = {fastq_file1}")
    print(f"\t - Fastq File 2 = {fastq_file2}")
    print(f"\t - Force = {force}")
    print(f"\t - Threads = {threads}")

    print("*START Align BWA*")
    sam_file = alignment_bwa(fasta_file, fastq_file1, fastq_file2, threads, force)
    print("*END Align BWA*")

    print("*START Convert SAM to BAM*")
    bam_file = convert_sam_to_bam(sam_file, force)
    print("*END Convert SAM to BAM*")

    print("*START Split unmapped and mapped reads*")
    bam_file, unmapped_fastq_file = split_unmapped_mapped_reads(bam_file, force)
    print("*END Split unmapped and mapped reads*")

    print("*START sort BAM*")
    bam_file = sort_bam_file(bam_file)
    print("*END sort BAM*")

    print("*START index BAM*")
    index_bam_file(bam_file)
    print("*END index BAM*")


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='tsv2fasta - Version ' + version())
    parser.add_argument('-f', '--fasFile', dest="fas_file",
                        default='/media/bacteriologie/TX/NGS-caen-enterobacter/fosfo/CNR1558/CNR1558.fasta',
                        help='DNA sequence as multi-fasta file')
    parser.add_argument('-fq1', '--fastq1', dest="fastq_file1",
                        default='/media/bacteriologie/TX/NGS-caen-enterobacter/fosfo/CNR1558/CNR1558_R1.fastq.gz',
                        help='Forward paired end fastq/fastq.gz file')
    parser.add_argument('-fq2', '--fastq2', dest="fastq_file2",
                        default='/media/bacteriologie/TX/NGS-caen-enterobacter/fosfo/CNR1558/CNR1558_R2.fastq.gz',
                        help='Reverse paired end fastq/fastq.gz file')
    parser.add_argument('-F', '--Force', dest="Force", action="store_true", default=False,
                        help="Overwrite the detected files")
    parser.add_argument('-t', '--threads', dest="threads", default='8', help="Number of threads [0]")
    parser.add_argument('-V', '--version', action='version', version='fasta2bam-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
