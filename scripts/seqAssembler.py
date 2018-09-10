#!/usr/bin/python3
# Write by Richard Bonnet
# Date: 19/01/2016
import os
import subprocess
import argparse
from shutil import copy2

from Bio import SeqIO
import time

from Bio.SeqIO.FastaIO import SimpleFastaParser

from seqassembler.scripts import a5, fasta2bam, bam2stats

install_dir = os.path.dirname(os.path.realpath(__file__))


def setup_samples(sample_file):
    # sample_dict = {}
    sample_list = []
    for line in open(sample_file, 'r'):
        try:
            prefix = line.split('\t')[0]
            species = line.strip().split('\t')[1].lower()
            print('  => ID:', prefix, 'Species', species)
            # sample_dict[prefix] = species
            sample_list.append(prefix)
        except IndexError:
            print('File format problem!')
            print('Required format:')
            print('prefix\tgenre species')
            exit()

    return sample_list


def fq_files(sample, fq_dir):
    # SEARCH FASTQ OR FASTQ.GZ
    file_list = []
    for ext in ['.fastq', '.fastq.gz']:
        for dir_path, dir_names, file_names in os.walk(fq_dir):
            file_list = file_list + [os.path.join(dir_path, filename) for filename in file_names
                                     if (filename.endswith(ext) and sample == os.path.basename(filename).split('_')[0])]
    file_list.sort()
    return file_list


def check_trimming(job_dir):
    sk_pe_trim_list = [os.path.join(job_dir, 'sk_s1_pe.fastq.gz'),
                       os.path.join(job_dir, 'sk_s2_pe.fastq.gz'),
                       os.path.join(job_dir, 'sk_s3_up.fastq.gz')]

    sk_se_trim_list = [os.path.join(job_dir, 'sk_s1_se.fastq.gz')]

    tr_pe_trim_list = [os.path.join(job_dir, 'tr_s1_pe.fastq.gz'),
                       os.path.join(job_dir, 'tr_s1_up.fastq.gz'),
                       os.path.join(job_dir, 'tr_s2_pe.fastq.gz'),
                       os.path.join(job_dir, 'tr_s2_up.fastq.gz')]

    tr_se_trim_list = [os.path.join(job_dir, 'tr_s1_se.fastq.gz')]

    for trimming_output_list in [sk_pe_trim_list, tr_pe_trim_list, sk_se_trim_list, tr_se_trim_list]:
        trimming_output_file_found = False
        for trimming_output_file in trimming_output_list:
            if os.path.exists(trimming_output_file):
                trimming_output_file_found = True
                break
        if trimming_output_file_found:
            print('\nTrimming output file:', trimming_output_list)
            print('Trimming already done!\n')
            break
    return trimming_output_file_found


def launch_trimming(job_dir, fq_list, subset, trimmer):
    job_dir = os.path.abspath(job_dir)

    # TRIMMING
    cmd = '{0} -f {1} -r {2} -o {3}'.format(os.path.join(install_dir, 'trimmer.py'), fq_list[0], fq_list[1], job_dir)

    if subset != 'all':
        cmd = cmd + ' -s {0}'.format(subset)
    if trimmer == 'trimmomatic':
        cmd = cmd + ' -tr'
    else:
        cmd = cmd + ' -sk'

    print('\nTrimming launcher:\n{0}\n'.format(cmd))
    out_str = subprocess.check_output(cmd, shell=True)
    print(out_str)


def launch_a5(sample, job_dir, fq_list, force):
    job_dir = os.path.abspath(job_dir)

    # LAUNCH A5
    ass_dir = os.path.join(job_dir, 'a5')
    if not os.path.exists(ass_dir) or force:
        print('\nA5 launcher:')
        a5.main(fq_list[0], fq_list[1], sample, ass_dir)
    else:
        print('\nAssembly a5 already done!\n')


def launch_spades(assembler, sample, job_dir, fastq_dir, force):
    job_dir = os.path.abspath(job_dir)
    fastq_dir = os.path.abspath(fastq_dir)
    # SEARCH FASTQ AND FASTQ.GZ IN TRIM_DIR
    file_list = []
    for ext in ['.fastq', '.fastq.gz']:
        for dir_path, dir_names, file_names in os.walk(job_dir):
            file_list = file_list + [os.path.join(dir_path, filename) for filename in file_names if
                                     filename.endswith(ext)]
    file_list.sort()
    print('\nListing of assembly files:')
    print(file_list)

    # SORT FASTQ AND FASTQ.GZ IN TRIM_DIR
    sk_pe_list = []
    sk_se_list = []
    sk_up_list = []
    tr_pe_list = []
    tr_se_list = []
    tr_up_list = []
    for filename in file_list:
        if os.path.basename(filename).startswith('sk_s') and 'pe.fastq' in os.path.basename(filename):
            sk_pe_list.append(filename)
        elif os.path.basename(filename).startswith('sk_s') and 'se.fastq' in os.path.basename(filename):
            sk_se_list.append(filename)
        elif os.path.basename(filename).startswith('sk_s') and 'up.fastq' in os.path.basename(filename):
            sk_up_list.append(filename)
        elif os.path.basename(filename).startswith('tr_s') and 'pe.fastq' in os.path.basename(filename):
            tr_pe_list.append(filename)
        elif os.path.basename(filename).startswith('tr_s') and 'se.fastq' in os.path.basename(filename):
            tr_se_list.append(filename)
        elif os.path.basename(filename).startswith('tr_s') and 'up.fastq' in os.path.basename(filename):
            tr_up_list.append(filename)

    # SEARCH FASTA IN FASTQ_DIR
    long_reads = []
    for ext in ['.fasta']:
        for dir_path, dir_names, file_names in os.walk(fastq_dir):
            long_reads = long_reads + [os.path.join(dir_path, filename) for filename in file_names if
                                       filename.endswith(ext) and sample in filename]
    if long_reads:
        long_reads = long_reads[0]
    else:
        long_reads = ''

    # LAUNCH SPADES ASSEMBLER
    print('\nSPAdes launcher:', end=' ')
    cmd = os.path.join(install_dir, 'spades.py')
    if assembler == 'plasmidspades':
        cmd = cmd + ' --plasmid'
        print('Plasmid assembly with', end=' ')
    else:
        print('Genomic assembly with', end=' ')

    if not os.path.exists(os.path.join(job_dir, assembler)) or force:

        if sk_pe_list:
            sk_pe_list.sort()
            cmd = cmd + ' -1 {0} -2 {1}'.format(sk_pe_list[0], sk_pe_list[1])
            if sk_up_list:
                cmd = cmd + ' -s {0}'.format(','.join(sk_up_list))
            if long_reads != '':
                cmd = cmd + ' -tc {0}'.format(long_reads)
            cmd = cmd + ' -o {0}'.format(job_dir)
            print('PE sickle files\n{0}'.format(cmd))
            print('SPAdes in process...')
            # os.system(cmd)
            out_str = subprocess.check_output(cmd, shell=True)
            print(out_str)

        elif tr_pe_list:
            tr_pe_list.sort()
            cmd = cmd + ' -1 {0} -2 {1}'.format(tr_pe_list[0], tr_pe_list[1])
            if tr_up_list:
                cmd = cmd + ' -s {0}'.format(','.join(tr_up_list))
            if long_reads != '':
                cmd = cmd + ' -tc {0}'.format(long_reads)
            cmd = cmd + ' -o {0}'.format(job_dir)
            print('PE Trimmomatic files\n{0}'.format(cmd))
            print('SPAdes in process...')
            # os.system(cmd)
            out_str = subprocess.check_output(cmd, shell=True)
            print(out_str)

        elif sk_se_list:
            cmd = cmd + ' -s {0}'.format(sk_se_list[0])
            if long_reads != '':
                cmd = cmd + ' -tc {0}'.format(long_reads)
            cmd = cmd + ' -o {0}'.format(job_dir)
            print('SE sickle files\n{0}'.format(cmd))
            print('SPAdes in process...')
            # os.system(cmd)
            out_str = subprocess.check_output(cmd, shell=True)
            print(out_str)

        elif tr_se_list:
            cmd = cmd + ' -s {0}'.format(tr_se_list[0])
            if long_reads != '':
                cmd = cmd + ' -tc {0}'.format(long_reads)
            cmd = cmd + ' -o {0}'.format(job_dir)
            print('SE Trimmomatic files\n{0}'.format(cmd))
            print('SPAdes in process...')
            # os.system(cmd)
            out_str = subprocess.check_output(cmd, shell=True)
            print(out_str)
    else:
        print('\nAssembly {0} already done!\n'.format(assembler))


def select_assembly(assembler, job_dir, sample, min_size, input_assembler_list):
    source_file = ""
    if assembler == 'spades':
        source_file = os.path.join(job_dir, 'spades', 'scaffolds.fasta')
    elif assembler == 'a5':
        source_file = os.path.join(job_dir, 'a5', sample + '.final.scaffolds.fasta')

    # name of final fasta assembly
    destination_file = os.path.join(job_dir, sample + '.fasta')

    print('\n\n## Assembly with {0} done!  ##'.format(assembler))
    print('Quality check:')
    s_quality_dic = assess_quality(source_file)
    print('Genome size: {0}'.format(s_quality_dic['assembly_len']))
    print('N number: {0}'.format(s_quality_dic['N_number']))
    print('Percentage of N(s): {0}'.format(
        round((100 * s_quality_dic['N_number']) / float(s_quality_dic['assembly_len'])), 5))
    print('Number of scaffold (length >= 500 pb): {0}'.format(s_quality_dic['contig_number']))
    print('N50: {0}'.format(s_quality_dic['N50']))
    # time.sleep(15)

    copy = True

    # check if in case of two or more assembler called wich is the best assembly
    if os.path.exists(destination_file) and len(input_assembler_list)>=2 :
        # get N50
        s_n50 = s_quality_dic['N50']

        # get count of N bases
        s_N = s_quality_dic['N_number']
        d_quality_dic = assess_quality(destination_file)
        d_n50 = d_quality_dic['N50']
        d_n = d_quality_dic['N_number']
        # check N50
        if s_n50 < d_n50:
            copy = False
            print("More N bases than N50")
        elif s_n50 == d_n50:
            if s_N < d_n:
                copy = False
        else:
            # N50 is good
            copy = True
            pass

    # if the newest assembly produce can be copy
    if copy:
        print('\nThe assembly going to be written in {0}'.format(destination_file))
        copy2(source_file, destination_file)
        print('The assembly is written in {0}'.format(destination_file))

        rec_list = []

        print('\nThe assembly contigs are going to be rename in "ctg_/d+" format')
        # get id of sorted by length the contigs
        with open(destination_file, 'r') as f:
            len_and_ids = sorted(((len(seq), title.split(None, 1)[0]) for
                                  title, seq in SimpleFastaParser(f)), reverse=True)
            ids = [id for (length, id) in len_and_ids]
            del len_and_ids  # free this on memory

        with open(destination_file, 'r') as f:
            n = 0
            records = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
            for id in ids:
                rec = records.get(id)
                if len(rec.seq) >= min_size:
                    n += 1
                    rec.id = 'ctg_{0}'.format(n)
                    rec.description = ''
                    rec_list.append(rec)
        SeqIO.write(rec_list, open(destination_file, 'w'), 'fasta')

        print('The assembly contigs have been renamed !')
    else:
        print('The assembly was not selected as the best one by the N50 values')
        exit(1)
    print('##  END  ##\n')

    return destination_file


def assess_quality(filename):
    with open(filename, 'r') as f:
        records = []
        assembly_len = 0
        n_number = 0
        for rec in SeqIO.parse(f, 'fasta'):
            n_number = n_number + str(rec.seq).count('N') + str(rec.seq).count('n')
            assembly_len += len(str(rec.seq))
            records.append(str(rec.seq))
    contig_number = len([x for x in records if len(x) >= 500])
    records.sort(key=len)
    total = 0
    n50 = ""
    for seq in records:
        total += len(seq)
        n50 = len(seq)
        if total >= (assembly_len / 2.0):
            break
    return {'assembly_len': assembly_len, 'N50': n50, 'contig_number': contig_number, 'N_number': n_number}


def launch_fasta2bam(fasta_file, fq_list):
    print('\n#Run the alignment of reads with assembly')

    fasta2bam.main(fasta_file, fq_list[0], fq_list[1], 8, True)

    print('\n#End bam file generated, sorted and indexed')

    return os.path.splitext(fasta_file)[0] + '.bam'


def launch_bam2stats(fasta_file, bam_file, mbam_depth, mbam_size, mbam_basq, mbam_mapq):
    print('\n#Run bam-based filter')

    bam2stats.main(fas_file=fasta_file, bam_file=bam_file, m_depth=mbam_depth, m_size=mbam_size,
                   m_basq=mbam_basq, m_mapq=mbam_mapq, filter_bam=True)

    # copy fasta file
    copy2(os.path.join(os.path.dirname(fasta_file), 'bam_stats', 'assembly_filtered.fas'), fasta_file)

    print('\n#End bam file generated, sorted and indexed')


def launch_plasflow(destination_file, outfile_prep, outfile_plasflow, threshold):
    print('\n\n#Preparation of PlasFlow sequence')

    cmd = 'perl /usr/local/PlasFlow/filter_sequences_by_length.pl -input {0} -output {1} -thresh {2}'.format(
        destination_file, outfile_prep, threshold)
    out_str = subprocess.check_output(cmd, shell=True)
    print(out_str)

    print('\n#End of the preparation of PlasFlow sequence')

    print('\n#Run of PlasFlow')

    cmd = 'python3 /usr/local/PlasFlow/PlasFlow.py --input {0} --output {1}'.format(outfile_prep, outfile_plasflow)
    out_str = subprocess.check_output(cmd, shell=True)
    print(out_str)

    print('\n#End of PlasFlow')

    return outfile_plasflow


def main(args):
    trimmer_list = ['sickle', 'trimmomatic']
    assembler_list = ['a5', 'spades', 'plasmidspades']

    # SETUP FASTQ/FASTGZ DIRECTORY
    fq_dir = args.fqDir
    if not os.path.exists(fq_dir):
        print('\nFastq/Fastq.gz directory was not found: {0}\n'.format(fq_dir))
        exit(1)

    # SETUP THE OUTPUT DIRECTORY
    out_dir = args.outDir
    if out_dir == '':
        out_dir = os.getcwd()
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # SETUP THE PREFIX AND SPECIES OF INPUT FILES AS A LIST FROM SAMPLE INPUT FILE
    sample_file = args.sampleFile
    if sample_file == '':
        sample_file = os.path.join(out_dir, 'sample.csv')
    if not os.path.exists(sample_file):
        print('\nSample file directory was not found: {0}\n'.format(sample_file))
        exit(1)
    sample_list = setup_samples(sample_file)

    # SETUP TRIMMING-fd
    trimmer = args.trimmer

    # SETUP SUBSET
    subset_size = args.subset

    # START THE JOBS
    for sample in sample_list:
        job_dir = os.path.abspath(os.path.join(out_dir, sample))

        if not os.path.exists(job_dir):
            os.makedirs(os.path.abspath(job_dir))

        print('\nStart with {0} in {1}'.format(sample, job_dir))

        # fastq or fastq.gz files
        fq_list = fq_files(sample, fq_dir)

        # SETUP ASSEMBLER
        print('\nLaunch assembly\n')
        input_assembler_list = args.assembler.split(',')
        for assembler in input_assembler_list:
            if assembler not in assembler_list:
                print('\nInvalid assembler name: {0}\n'.format(assembler))
                exit(1)

            elif assembler == 'a5':
                launch_a5(sample, job_dir, fq_list, args.force)

            elif assembler == 'spades' or assembler == 'plasmidspades':
                trimming_output_file_found = check_trimming(job_dir)
                if not trimming_output_file_found or args.force:
                    if trimmer in trimmer_list:
                        launch_trimming(job_dir, fq_list, subset_size, trimmer)
                    else:
                        print('\nTrimmer {0} not found\n'.format(trimmer))
                        exit(1)
                launch_spades(assembler, sample, job_dir, fq_dir, args.force)

            # Make Assembly file with filtering quality
            destination_file = select_assembly(assembler, job_dir, sample, int(args.minSize), input_assembler_list)

            # Bam file
            if args.Bam:
                bam_file = launch_fasta2bam(destination_file, fq_list)
                if not args.nobamFilter:
                    launch_bam2stats(destination_file, bam_file, args.mbamDepth, args.mbamSize, args.mbamBasq,
                                     args.mbamMapq)

            # Make differentiation between chrom vs plasmid
            if args.plasFlow:
                threshold = 0.7
                outfile_prep = os.path.join(out_dir, sample + "filtered_plasflow.fasta")
                outfile_plasflow = os.path.join(out_dir, sample + "plasflow_predictions")
                launch_plasflow(destination_file, outfile_prep, outfile_plasflow, threshold)


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='seqAssembler - Version ' + version())
    parser.add_argument('-fq', '--fqDir', dest="fqDir", help='The directory containing fastq or fastq.gz files')
    parser.add_argument('-o', '--outDir', dest="outDir", default='./', help="The output directory name")
    parser.add_argument('-s', '--sampleFile', dest="sampleFile",
                        help="The sample file with sample names (default: $outDir/sample.csv)")
    parser.add_argument('-tr', '--trimmer', dest="trimmer", default='sickle',
                        help="sickle or trimmomatic or nothing (default: sickle)")
    parser.add_argument('-a', '--assembler', dest="assembler", default="a5",
                        help="Assembler names [a5,spades,plasmidspades] as a comma separated list (default: a5)")
    parser.add_argument('-sb', '--subset', dest="subset", default='all',
                        help="The number of loaded reads by fastq file in trimmomatic-based trimming (default: all)")
    parser.add_argument('-ms', '--minSize', dest="minSize", default=500,
                        help="Min size of scaffolds in assembly (default: 500)")
    parser.add_argument('-b', '--Bam', dest="Bam", action='store_true', default=False, help="Create bam file")
    parser.add_argument('-nbf', '--nobamFilter', dest="nobamFilter", action='store_true', default=False,
                        help="Do not filter from bam file")
    parser.add_argument('-mbd', '--mbamDepth', dest="mbamDepth", default=20,
                        help="Minimum mean depth threshold for the scaffolds (default: 20)")
    parser.add_argument('-mbs', '--mbamSize', dest="mbamSize", default=500,
                        help="Minimum size of scaffolds in bam filtering (default: 500)")
    parser.add_argument('-mbq', '--mbamBasq', dest="mbamBasq", default=20,
                        help="Minimum mean base quality threshold for the scaffolds (default: 20)")
    parser.add_argument('-mmq', '--mbamMapq', dest="mbamMapq", default=30,
                        help="Minimum mapping quality threshold for the scaffolds (default: 30)")
    parser.add_argument('-pf', '--plasFlow', dest="plasFlow", action='store_true', default=False,
                        help="launch the identification of plasmid/chromosome contigs (default: False)")
    parser.add_argument('-F', '--force', dest="force", action='store_true', default=False,
                        help="Force file overwrite (default: False)")

    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
