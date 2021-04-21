#!/usr/bin/python
import logging
import os
import subprocess
import argparse

from seqassembler_lib.seqassembler.fasta2bam import log_process_output


def launch(sample, file1, file2, out_dir):
    logging.info(f'\nAssembly of {sample} in {out_dir} with {file1} and {file2}')
    out_dir = os.path.abspath(out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    logging.info(f'File1: {file1}')
    logging.info(f'File2: {file2}')
    if file1.split('_R1') == file2.split('_R2'):
        logging.info('In process...')
        # set current directory to SKESA work directory
        os.chdir(out_dir)
        # Get version of SKESA
        cmd = 'skesa --version'
        # launch SKESA for version
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
        log = process.decode("utf-8")
        version_skesa = ""
        for n in log.split("\n"):
            logging.info(n)
            if "SKESA" in n:
                version_skesa = n.split(" ")[1]
        logging.info(f"\nVersion SKESA :{version_skesa}\n")

        output_assembly = os.path.join(out_dir, f"{sample}.skesa.fa")
        cmd = f'skesa --reads {file1},{file2} --cores 4 --memory 20 > {output_assembly}'
        logging.info(cmd)

        # launch SKESA
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()

        # make log
        filename_log = "log_skesa_pipeline.txt"
        header = f"Command line executed: {cmd}\n\n\n{process.decode('utf-8')}"
        log_process_output(header, out_dir, filename_log)

        current_dir = os.getcwd()
        # remove useless files
        pivot = 0
        for f in os.listdir(current_dir):
            if os.path.isdir(f):
                continue
            else:
                if f"{sample}.skesa.fa" == f:
                    pivot = 1

        if pivot == 1:
            logging.info(f'Assembly of {sample} done!')
        else:
            logging.error(f'Assembly of {sample} not done! Check error file log : {filename_log}')


def pre_main(arguments):
    file1 = arguments.file1
    file2 = arguments.file2
    sample = arguments.sample
    out_dir = arguments.outdir
    main(file1, file2, sample, out_dir)


def main(file1, file2, sample, out_dir):
    print(f'\nSample: {sample}')
    print(f'Input file names: {file1} {file2}')
    print(f'Output dir: {out_dir}\n')
    launch(sample, file1, file2, out_dir)


def version():
    return "0.0.1"


def run():
    parser = argparse.ArgumentParser(description='launch SKESA assembler - Version ' + version())
    parser.add_argument('-file1', '--readFile_input1', dest='file1', default='file1.fastq',
                        help='Forward fastq or fastq.gz file')
    parser.add_argument('-file2', '--readFile_input2', dest='file2', default='file2.fastq',
                        help='Reverse fastq or fastq.gz file')
    parser.add_argument('-id', '--sampleName', dest='sample', help="ID of sample")
    parser.add_argument('-out', '--outputDir', dest='outdir', help="Name of output directory")
    parser.add_argument('-V', '--version', action='version', version='rgi-' + version(), help="Prints version number")
    return parser.parse_args()


if __name__ == '__main__':
    args = run()
    pre_main(args)
