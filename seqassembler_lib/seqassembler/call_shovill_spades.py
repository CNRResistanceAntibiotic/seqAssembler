#!/usr/bin/python
import logging
import os
import argparse
import subprocess

from seqassembler_lib.seqassembler.fasta2bam import log_process_output


def launch(sample, pe_file1, pe_file2, out_dir):
    # Get version Spades
    cmd = 'spades.py -v'
    # launch spades for version
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    log = process.decode("utf-8")
    logging.info(f"\nVersion Spades :{log.split('SPAdes genome assembler ')[1]}\n")

    # Get version of Shovill
    cmd = 'shovill --version'
    # launch Shovill for version
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    log = process.decode("utf-8")
    version_shovill = ""
    for n in log.split("\n"):
        if "shovill" in n:
            version_shovill = n.split(" ")[1]
    logging.info(f"\nVersion Shovill-spades :{version_shovill}\n")

    cmd = f'shovill --assembler spades --R1 {pe_file1} --R2 {pe_file2} --outdir {out_dir} --ram 20 --force'

    logging.info(f'\nShovill Spades:\n{cmd}\n')
    logging.info('Shovill Spades in process...')

    # launch Shovill Spades
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()

    # make log
    filename_log = "log_shovill_spades_pipeline.txt"
    header = f"Command line executed: {cmd}\n\n\n{process.decode('utf-8')}"
    log_process_output(header, out_dir, filename_log)

    os.remove(os.path.join(out_dir, 'contigs.gfa'))
    os.remove(os.path.join(out_dir, 'spades.fasta'))

    if os.path.exists(os.path.join(out_dir, "contigs.fa")):
        logging.info(f'Assembly of {sample} done!')
    else:
        logging.info(f'Assembly of {sample} not done! Check error file log : {filename_log}')


def pre_main(arguments):
    pe_file1 = arguments.pefile1
    pe_file2 = arguments.pefile2
    s_files = arguments.upfiles
    pacbio = arguments.pacbio
    sanger = arguments.sanger
    tr_contig = arguments.trcontig
    un_contig = arguments.uncontig
    out_dir = arguments.outdir
    plasmid = arguments.plasmid
    cv = arguments.cv
    main("toto", pe_file1, pe_file2, out_dir)


def main(file1, file2, sample, out_dir):
    logging.info(f'\nSample: {sample}')
    logging.info(f'Input file names: {file1} {file2}')
    logging.info(f'Output dir: {out_dir}\n')
    launch(sample, file1, file2, out_dir)


def version():
    return "0.0.1"


def run():
    parser = argparse.ArgumentParser(description='launch shovill spades assembler - Version ' + version())
    parser.add_argument('-1', '--forwardPE', dest='pefile1', action='store', default='sk_s1_pe.fastq',
                        help='Forward fastq or fastq.gz file')
    parser.add_argument('-2', '--reversePE', dest='pefile2', action='store', default='sk_s2_pe.fastq',
                        help='Reverse fastq or fastq.gz file')
    parser.add_argument('-s', '--unpaired', dest='upfiles', action='store', default='sk_s3_up.fastq',
                        help='Unpaired fastq or fastq.gz file')
    parser.add_argument('-pb', '--pacbio', dest='pacbio', action='store', default='', help='PacBio file')
    parser.add_argument('-sg', '--sanger', dest='sanger', action='store', default='', help='Sanger file')
    parser.add_argument('-tc', '--trustedcontig', dest='trcontig', action='store', default='',
                        help='Trusted contig file')
    parser.add_argument('-uc', '--untrustedContig', dest='uncontig', action='store', default='',
                        help='UnTrusted contig file')
    parser.add_argument('-pl', '--plasmid', dest='plasmid', action='store_true', help='Plasmide assembly')
    parser.add_argument('-cv', '--cv', dest='cv', action='store', default='off',
                        help='Coverage cutoff value (a positive float number or \'auto\' or \'off\')'
                             ' [default: \'off\']')
    parser.add_argument('-o', '--outputDir', dest='outdir', action='store', help="Name of output directory")
    parser.add_argument('-V', '--version', action='version', version='rgi-' + version(), help="Prints version number")
    return parser.parse_args()


if __name__ == '__main__':
    args = run()
    pre_main(args)
