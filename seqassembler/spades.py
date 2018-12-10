#!/usr/bin/python
import os
import argparse
import subprocess
from shutil import rmtree


def launch(plasmid, cv, pe_file1, pe_file2, s_files, pacbio, sanger, trcontig, uncontig, out_dir):

    # Get version Spades
    cmd = 'spades.py -v'
    # launch spades for version
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    log = process.decode("utf-8")
    print("\nVersion Spades :{0}\n".format(log.split(" ")[1]))

    if plasmid:
        ass_dir = os.path.join(out_dir, 'plasmidspades')
        i = 0
        while os.path.exists(ass_dir):
            i += 1
            ass_dir = os.path.join(out_dir, 'plasmidspades_{0}'.format(i))

        cmd = 'spades.py --plasmid --careful'
    else:
        ass_dir = os.path.join(out_dir, 'spades')
        i = 0
        while os.path.exists(ass_dir):
            i += 1
            ass_dir = os.path.join(out_dir, 'spades_{0}'.format(i))
        cmd = 'spades.py --careful'

    if pe_file1 != '' and pe_file2 != '':
        cmd = cmd + ' -1 {0} -2 {1}'.format(pe_file1, pe_file2)

    for index, f in enumerate(s_files):
        index += 1
        if f != '':
            cmd = cmd + ' --s{0} {1}'.format(index, f)

    if pacbio != '':
        cmd = cmd + ' --pacbio {0}'.format(pacbio)

    if sanger != '':
        cmd = cmd + ' --sanger {0}'.format(sanger)

    if trcontig != '':
        cmd = cmd + ' --trusted-contigs {0}'.format(trcontig)

    if uncontig != '':
        cmd = cmd + ' --untrusted-contigs {0}'.format(uncontig)

    if cv != 'off':
        cmd = cmd + ' --cov-cutoff {0}'.format(cv)

    cmd = cmd + ' -o {0}'.format(ass_dir)
    print('\nSpades:\n{0}\n'.format(cmd))
    print('Spades in process...')
    subprocess.check_output(cmd, shell=True)

    files_to_remove = ["assembly_graph.gfa", "assembly_graph.fastg", "before_rr.fasta"]
    for file in os.listdir(ass_dir):

        file_path = os.path.join(ass_dir, file)

        if os.path.isfile(file_path):

            if file in files_to_remove:
                os.remove(file_path)
        if os.path.isdir(file_path):
            rmtree(file_path)


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
    main(pe_file1, pe_file2, s_files, pacbio, sanger, tr_contig, un_contig, out_dir, plasmid, cv)


def main(pe_file1="", pe_file2="", s_files="", pacbio="", sanger="", tr_contig="", un_contig="", out_dir="",
         plasmid=False, cv="off"):

    s_files = s_files.split(',')

    print('\nOutput dir: {0}'.format(out_dir))
    print('Input file names: {0} {1} {2} {3} {4} {5} {6}'.format(
        pe_file1, pe_file2, ' '.join(s_files), pacbio, sanger, tr_contig, un_contig))
    # backup_assembly(out_dir, sample)
    launch(plasmid, cv, pe_file1, pe_file2, s_files, pacbio, sanger, tr_contig, un_contig, out_dir)


def version():
    return "0.0.1"


def run():
    parser = argparse.ArgumentParser(description='launch A5 assembler - Version ' + version())
    parser.add_argument('-1', '--forwardPE', dest='pefile1', action='store', default='sk_s1_pe.fastq',
                        help='Forward fastq or fastqz file')
    parser.add_argument('-2', '--reversePE', dest='pefile2', action='store', default='sk_s2_pe.fastq',
                        help='Reverse fastq or fastqz file')
    parser.add_argument('-s', '--unpaired', dest='upfiles', action='store', default='sk_s3_up.fastq',
                        help='Unpaired fastq or fastqz file')
    parser.add_argument('-pb', '--pacbio', dest='pacbio', action='store', default='', help='PacBio file')
    parser.add_argument('-sg', '--sanger', dest='sanger', action='store', default='', help='Sanger file')
    parser.add_argument('-tc', '--trustedcontig', dest='trcontig', action='store', default='',
                        help='Trusted contig file')
    parser.add_argument('-uc', '--untrustedContig', dest='uncontig', action='store', default='',
                        help='UnTrusted contig file')
    parser.add_argument('-pl', '--plasmid', dest='plasmid', action='store_true', help='Plamide assembly')
    parser.add_argument('-cv', '--cv', dest='cv', action='store', default='off',
                        help='Coverage cutoff value (a positive float number or \'auto\' or \'off\')'
                             ' [default: \'off\']')
    parser.add_argument('-o', '--outputDir', dest='outdir', action='store', help="Name of output directory")
    parser.add_argument('-V', '--version', action='version', version='rgi-' + version(), help="Prints version number")
    return parser.parse_args()


if __name__ == '__main__':
    args = run()
    pre_main(args)
