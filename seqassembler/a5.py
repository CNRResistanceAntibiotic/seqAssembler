#!/usr/bin/python
import glob
import os
import subprocess
import argparse
from shutil import rmtree


def backup_assembly(out_dir, sample):
    # BACKUP CONTIG, SCAFFOLD AND STAT FILES
    for ext in ['.assembly_stats.csv', '.contigs.fasta', '.contigs.fastq', '.final.scaffolds.fasta',
                '.final.scaffolds.fastq', 'pe.sort.bam', '.pe.sort.bam.bai']:
        for filename in glob.glob(os.path.join(out_dir, sample + ext)):
            if filename:
                print('Previous assembly file {0} detected'.format(filename))
                print('Backup of file {0} as file {1}_previous'.format(filename, filename))
                os.rename(filename, filename + '_previous')


def launch(sample, file1, file2, out_dir):
    print('\nAssembly of {0} in {1} with {2} and {3}'.format(sample, out_dir, file1, file2))
    out_dir = os.path.abspath(out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print('File1: {0}'.format(file1))
    print('File2: {0}'.format(file2))
    if file1.split('_R1') == file2.split('_R2'):
        print('In process...')

        # set current directory to a5 work directory
        os.chdir(out_dir)

        arguments = " --end=5 {0} {1} {2}".format(file1, file2, sample)

        # launch a5
        subprocess.check_call('a5_pipeline.pl' + arguments, shell=True)

        current_dir = os.getcwd()
        # remove useless files
        for f in os.listdir(current_dir):
            if os.path.isdir(f):
                rmtree(f)
            else:
                if "contigs.fasta" not in f or ".final.scaffolds.fast" not in f:
                    os.remove(os.path.join(current_dir, f))

        print('Assembly of {0} done!'.format(sample))

        inp_stat_file = os.path.join(out_dir, sample + '.assembly_stats.csv')
        if os.path.exists(inp_stat_file):
            out_stat_file = os.path.join(os.path.dirname(out_dir), 'a5_assembly_stats.csv')
            header = ""
            if os.path.exists(out_stat_file):
                IO_type = "a"
            else:
                IO_type = "w"
                header = 'ID\tNombre de contigs\tNombre de scaffolds\tTaille du genome\tScaffold le plus long\tN50\t' \
                         'Nombre de reads\tNombre de reads conserves\t% de reads conservees\tNombre de nucleotides\t' \
                         'Nombre de nucleotides conserves\t% de nucleotides conserves\t' \
                         'Profondeur moyenne\tProfondeur moyenne filtree\tProfondeur mediane\t' \
                         'Profondeur au 10eme percentile\t' \
                         'Nombre de bases >= Q40\tGC %\n'

            with open(out_stat_file, '{0}'.format(IO_type)) as f:
                if header:
                    f.write(header)

                for n, line in enumerate(open(inp_stat_file, 'r')):
                    if n > 0:
                        f.write('\t'.join(line.split('\t')[:-1]) + '\n')


def pre_main(arguments):
    file1 = arguments.file1
    file2 = arguments.file2
    sample = arguments.sample
    out_dir = arguments.outdir
    main(file1, file2, sample, out_dir)


def main(file1, file2, sample, out_dir):
    print('\nSample: %s' % sample)
    print('Input file names: %s %s' % (file1, file2))
    print('Output dir: %s\n' % out_dir)
    backup_assembly(out_dir, sample)
    launch(sample, file1, file2, out_dir)


def version():
    return "0.0.1"


def run():
    parser = argparse.ArgumentParser(description='launch A5 assembler - Version ' + version())
    parser.add_argument('-file1', '--readFile_input1', dest='file1', default='file1.fastq',
                        help='Forward fastq or fastqz file')
    parser.add_argument('-file2', '--readFile_input2', dest='file2', default='file2.fastq',
                        help='Reverse fastq or fastqz file')
    parser.add_argument('-id', '--sampleName', dest='sample', help="ID of sample")
    parser.add_argument('-out', '--outputDir', dest='outdir', help="Name of output directory")
    parser.add_argument('-V', '--version', action='version', version='rgi-' + version(), help="Prints version number")
    return parser.parse_args()


if __name__ == '__main__':
    args = run()
    pre_main(args)
