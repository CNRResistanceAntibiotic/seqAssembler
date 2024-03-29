#!/usr/bin/python3
# Write by Richard Bonnet
# Date: 19/01/2016
import argparse
import os
import subprocess
import sys
from shutil import copy2, move

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

from seqassembler_lib.seqassembler import trimmer, call_a5, call_spades, bam2stats, fasta2bam, call_skesa, \
    call_shovill_skesa, call_shovill_spades, call_shovill_velvet, call_shovill_megahit


def setup_samples(sample_file):
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
    sk_pe_trim_list = [os.path.join(job_dir, 'sk_s1_pe.fastq.gz'), os.path.join(job_dir, 'sk_s2_pe.fastq.gz'),
                       os.path.join(job_dir, 'sk_s3_up.fastq.gz')]
    sk_se_trim_list = [os.path.join(job_dir, 'sk_s1_se.fastq.gz')]
    tr_pe_trim_list = [os.path.join(job_dir, 'tr_s1_pe.fastq.gz'), os.path.join(job_dir, 'tr_s1_up.fastq.gz'),
                       os.path.join(job_dir, 'tr_s2_pe.fastq.gz'), os.path.join(job_dir, 'tr_s2_up.fastq.gz')]
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


def launch_trimming(job_dir, fq_list, subset, trimmer_type):
    job_dir = os.path.abspath(job_dir)
    trimmer_dir = os.path.join(job_dir, f"trimming_{trimmer_type}")
    # create directory of trimming reads
    if not os.path.exists(trimmer_dir):
        os.mkdir(trimmer_dir)

    # TRIMMING
    tr = sk = False
    if trimmer_type == 'trimmomatic':
        tr = True
    else:
        sk = True

    print(f'\nTrimming launcher:\nin_f {fq_list[0]} in_r {fq_list[1]} output_dir {trimmer_dir} subset_size {subset}'
          f' trimmomatic {tr} sickle {sk}\n')
    trimmer.main(in_f=fq_list[0], in_r=fq_list[1], output_dir=trimmer_dir, subset_size=subset, trimmomatic=tr,
                 sickle=sk, sickle_g=True)
    return trimmer_dir


def launch_a5(sample, job_dir, fq_list, force):
    # LAUNCH A5
    ass_dir = os.path.join(job_dir, 'a5')
    if not os.path.exists(ass_dir) or force:
        print('\nA5 launcher:')
        call_a5.main(fq_list[0], fq_list[1], sample, ass_dir)
    else:
        print('\nAssembly a5 already done!\n')


def launch_skesa(sample, job_dir, fq_list, force):
    # LAUNCH SKESA
    ass_dir = os.path.join(job_dir, 'SKESA')
    if not os.path.exists(ass_dir) or force:
        print('\nSKESA launcher:')
        call_skesa.main(fq_list[0], fq_list[1], sample, ass_dir)
    else:
        print('\nAssembly SKESA already done!\n')


def launch_shovill_skesa(sample, job_dir, fq_list, force, temp_dir):
    # LAUNCH shovill SKESA
    ass_dir = os.path.join(job_dir, 'shovill-SKESA')
    if not os.path.exists(ass_dir) or force:
        print('\nShovill-SKESA launcher:')
        call_shovill_skesa.main(fq_list[0], fq_list[1], sample, ass_dir, temp_dir)
    else:
        print('\nAssembly Shovill-SKESA already done!\n')


def launch_shovill_spades(sample, job_dir, fq_list, force, temp_dir):
    # LAUNCH shovill spades
    ass_dir = os.path.join(job_dir, 'shovill-spades')
    if not os.path.exists(ass_dir) or force:
        print('\nShovill-spades launcher:')
        call_shovill_spades.main(fq_list[0], fq_list[1], sample, ass_dir, temp_dir)
    else:
        print('\nAssembly Shovill-spades already done!\n')


def launch_shovill_velvet(sample, job_dir, fq_list, force, temp_dir):
    # LAUNCH shovill velvet
    ass_dir = os.path.join(job_dir, 'shovill-velvet')
    if not os.path.exists(ass_dir) or force:
        print('\nShovill-velvet launcher:')
        call_shovill_velvet.main(fq_list[0], fq_list[1], sample, ass_dir, temp_dir)
    else:
        print('\nAssembly Shovill-velvet already done!\n')


def launch_shovill_megahit(sample, job_dir, fq_list, force, temp_dir):

    # LAUNCH shovill megahit
    ass_dir = os.path.join(job_dir, 'shovill-megahit')
    if not os.path.exists(ass_dir) or force:
        print('\nShovill-megahit launcher:')
        call_shovill_megahit.main(fq_list[0], fq_list[1], sample, ass_dir, temp_dir)
    else:
        print('\nAssembly Shovill-megahit already done!\n')


def launch_spades(assembler, sample, job_dir, fastq_dir, force, trimmer_dir, temp_dir):

    fastq_dir = os.path.abspath(fastq_dir)
    # SEARCH FASTQ AND FASTQ.GZ IN TRIM_DIR
    file_list = []
    for ext in ['.fastq', '.fastq.gz']:
        for dir_path, dir_names, file_names in os.walk(trimmer_dir):
            file_list = file_list + [os.path.join(dir_path, filename) for filename in file_names if
                                     filename.endswith(ext)]
    file_list.sort()
    print('\nListing of assembly files:')
    print(file_list)

    # SORT FASTQ AND FASTQ.GZ IN TRIM_DIR
    sk_pe_list, sk_se_list, sk_up_list, tr_pe_list, tr_se_list, tr_up_list = [], [], [], [], [], []
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
    #
    plasmid = False
    if assembler == 'plasmidspades':
        plasmid = True
        print('Plasmid assembly with', end=' ')
    else:
        print('Genomic assembly with', end=' ')

    if not os.path.exists(os.path.join(job_dir, assembler)) or force:
        s_files = tr_contig = ""
        if long_reads != '':
            tr_contig = long_reads

        if sk_pe_list:
            sk_pe_list.sort()
            if sk_up_list:
                s_files = ','.join(sk_up_list)

            print(f'PE sickle files\n{sk_pe_list[0]} {sk_pe_list[1]}')
            print('SPAdes in process...')
            call_spades.main(pe_file1=sk_pe_list[0], pe_file2=sk_pe_list[1], out_dir=job_dir, plasmid=plasmid,
                             s_files=s_files, tr_contig=tr_contig, temp_dir=temp_dir)

        elif tr_pe_list:
            tr_pe_list.sort()
            if tr_up_list:
                s_files = ','.join(tr_up_list)

            print(f'PE Trimmomatic files\n{tr_pe_list[0]} {tr_pe_list[1]}')
            print('SPAdes in process...')

            call_spades.main(pe_file1=tr_pe_list[0], pe_file2=tr_pe_list[1], out_dir=job_dir, plasmid=plasmid,
                             s_files=s_files, tr_contig=tr_contig, temp_dir=temp_dir)

        elif sk_se_list:
            s_files = sk_se_list[0]
            print(f'SE sickle files\n{sk_se_list[0]}')
            print('SPAdes in process...')
            call_spades.main(out_dir=job_dir, plasmid=plasmid, s_files=s_files, tr_contig=tr_contig, temp_dir=temp_dir)

        elif tr_se_list:
            print(f'SE Trimmomatic files\n{tr_se_list}')
            print('SPAdes in process...')
            call_spades.main(out_dir=job_dir, plasmid=plasmid, s_files=tr_se_list, tr_contig=tr_contig,
                             temp_dir=temp_dir)
    else:
        print(f'\nAssembly {assembler} already done!\n')


def select_assembly(job_dir, sample, min_size, input_assembler_list):
    source_file = final_assembler = ""
    # name of final fasta assembly
    destination_file = os.path.join(job_dir, sample + '.fasta')
    for assembler in input_assembler_list:
        if assembler == 'spades':
            source_file = os.path.join(job_dir, 'spades', 'scaffolds.fasta')
            if os.path.exists(source_file):
                print(f"The assembly file of {assembler} exist.")
            else:
                print(f"The assembly file of {assembler} not exist. This assembler is not keep for further steps.")
                continue

        elif assembler == 'a5':
            source_file = os.path.join(job_dir, 'a5', sample + '.final.scaffolds.fasta')

            if os.path.exists(source_file):
                print(f"The assembly file of {assembler} exist.")
            else:
                print(f"The assembly file of {assembler} not exist. This assembler is not keep for further steps.")
                continue

        elif assembler == 'SKESA':
            source_file = os.path.join(job_dir, 'SKESA', sample + '.skesa.fa')
            if os.path.exists(source_file):
                print(f"The assembly file of {assembler} exist.")
            else:
                print(f"The assembly file of {assembler} not exist. This assembler is not keep for further steps.")
                continue

        elif assembler == 'shovill-SKESA':
            source_file = os.path.join(job_dir, 'shovill-SKESA', 'contigs.fa')
            if os.path.exists(source_file):
                print(f"The assembly file of {assembler} exist.")
            else:
                print(f"The assembly file of {assembler} not exist. This assembler is not keep for further steps.")
                continue

        elif assembler == 'shovill-spades':
            source_file = os.path.join(job_dir, 'shovill-spades', 'contigs.fa')
            if os.path.exists(source_file):
                print(f"The assembly file of {assembler} exist.")
            else:
                print(f"The assembly file of {assembler} not exist. This assembler is not keep for further steps.")
                continue

        elif assembler == 'shovill-megahit':
            source_file = os.path.join(job_dir, 'shovill-megahit', 'contigs.fa')
            if os.path.exists(source_file):
                print(f"The assembly file of {assembler} exist.")
            else:
                print(f"The assembly file of {assembler} not exist. This assembler is not keep for further steps.")
                continue

        elif assembler == 'shovill-velvet':
            source_file = os.path.join(job_dir, 'shovill-velvet', 'contigs.fa')
            if os.path.exists(source_file):
                print(f"The assembly file of {assembler} exist.")
            else:
                print(f"The assembly file of {assembler} not exist. This assembler is not keep for further steps.")
                continue

        print(f'\n\n## Assembly with {assembler} done!  ##')
        print("Order Fasta assembled by Sequence length")

        with open(source_file, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        records.sort(key=lambda r: -len(r))
        source_file_corr = f"{source_file}-corrected"
        with open(source_file_corr, "w") as out_handle:
            SeqIO.write(records, out_handle, "fasta")

        os.remove(source_file)
        move(source_file_corr, source_file)
        print("End order fasta by length")
        print(f"Open {source_file} fasta file for check")
        print('Quality check:')
        s_quality_dic = assess_quality(source_file)
        print(f'Genome size: {s_quality_dic["assembly_len"]}')
        print(f'N number: {s_quality_dic["N_number"]}')
        print(f'Percentage of N(s): {round((100 * s_quality_dic["N_number"]) / float(s_quality_dic["assembly_len"]))}')
        print(f'Number of scaffold (length >= 500 pb): {s_quality_dic["contig_number"]}')
        print(f'N50: {s_quality_dic["N50"]}')

        copy = True

        # check if in case of two or more assembler called wich is the best assembly
        if os.path.exists(destination_file) and len(input_assembler_list) >= 2:
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
                print(f"The assembly with {assembler} have a lower N50 than the previous assembly done")
            elif s_n50 == d_n50:
                if s_N < d_n:
                    copy = False
            else:
                # N50 is good
                copy = True
                print(f"The assembly with {assembler} have a better N50 than the previous assembly done. "
                      f"It will be copy.")
                pass

        # if the newest assembly produce can be copy
        if copy:
            print(f'\nThe assembly going to be written in {destination_file}')
            copy2(source_file, destination_file)
            print(f'The assembly is written in {destination_file}')
            rec_list = []
            print('\nThe assembly contigs are going to be rename in "ctg_/d+" format')
            # get id_contig of sorted by length the contigs
            with open(destination_file, 'r') as f:
                len_and_ids = sorted(((len(seq), title.split(None, 1)[0]) for
                                      title, seq in SimpleFastaParser(f)), reverse=True)
                ids = [id_contig for (length, id_contig) in len_and_ids]
                del len_and_ids  # free this on memory

            with open(destination_file, 'r') as f:
                n = 1
                records = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
                for id_contig in ids:
                    rec = records.get(id_contig)
                    rec.id = f'ctg_{n}'
                    rec.description = ''
                    if len(rec.seq) >= min_size:
                        n += 1
                        rec_list.append(rec)
                    else:
                        print(f"Contig \"{id_contig}\" has been removed on the final assembly (contig length less"
                              f" than {min_size}bp (contig length = {len(rec.seq)}bp))")
            SeqIO.write(rec_list, open(destination_file, 'w'), 'fasta')
            final_assembler = assembler
            print('The assembly contigs have been renamed !')
        else:
            print(f"The assembly with {assembler} was not selected as the best one by the N50 values")
    if not final_assembler:
        print(f"\nNo assembly has been selected !!! The program exit.\n")
        exit()
    else:
        print(f"\nThe best assembly was done by {final_assembler}\n")
    print('##  END  ##\n')
    return destination_file


def assess_quality(filename):
    records = []
    assembly_len = total = n_number = 0
    with open(filename, 'r') as f:
        for rec in SeqIO.parse(f, 'fasta'):
            n_number = n_number + str(rec.seq).count('N') + str(rec.seq).count('n')
            assembly_len += len(str(rec.seq))
            records.append(str(rec.seq))
    contig_number = len([x for x in records if len(x) >= 500])
    records.sort(key=len)
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
    bam2stats.main(fas_file=fasta_file, bam_file=bam_file, m_depth=mbam_depth, m_size=mbam_size, m_basq=mbam_basq,
                   m_mapq=mbam_mapq, filter_bam=True, plt_report=True)
    # copy fasta file
    copy2(os.path.join(os.path.dirname(fasta_file), 'bam_stats', 'assembly_filtered.fas'), fasta_file)
    print('\n#End bam file generated, sorted and indexed')


def launch_plasflow(destination_file, outfile_prep, outfile_plasflow, threshold):
    print('\n\n#Preparation of PlasFlow sequence')
    cmd = f'perl /usr/local/PlasFlow/filter_sequences_by_length.pl -input {destination_file} -output {outfile_prep}' \
          f' -thresh {threshold}'
    out_str = subprocess.check_output(cmd, shell=True)
    print(out_str)
    print('\n#End of the preparation of PlasFlow sequence')
    #########################
    print('\n#Run of PlasFlow')
    cmd = f'python3 /usr/local/PlasFlow/PlasFlow.py --input {outfile_prep} --output {outfile_plasflow}'
    out_str = subprocess.check_output(cmd, shell=True)
    print(out_str)
    print('\n#End of PlasFlow')
    return outfile_plasflow


def main(args):
    trimmer_list = ['sickle', 'trimmomatic']
    assembler_list = ['a5', 'spades', 'plasmidspades', 'SKESA', 'shovill-spades', 'shovill-SKESA', 'shovill-velvet',
                      'shovill-megahit']

    # SETUP FASTQ/FASTGZ DIRECTORY
    fq_dir = args.fqDir
    if not os.path.exists(fq_dir):
        print(f'\nFastq/Fastq.gz directory was not found: {fq_dir}\n')
        exit(1)

    # SETUP THE OUTPUT DIRECTORY
    out_dir = args.outDir
    if out_dir == '':
        out_dir = os.getcwd()
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    temp_dir = args.tempDir
    if temp_dir:
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
    else:
        temp_dir = "/tmp"

    # SETUP THE PREFIX AND SPECIES OF INPUT FILES AS A LIST FROM SAMPLE INPUT FILE
    sample_file = args.sampleFile
    if sample_file == '':
        sample_file = os.path.join(out_dir, 'sample.csv')
    if not os.path.exists(sample_file):
        print(f'\nSample file directory was not found: {sample_file}\n')
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
        print(f'\nStart with {sample} in {job_dir}')
        # fastq or fastq.gz files
        fq_list = fq_files(sample, fq_dir)
        print(f'\nThe fastq files treated are (in PE) : {fq_list[0]} and {fq_list[1]}')
        # SETUP ASSEMBLER
        print('\nLaunch assembly\n')
        input_assembler_list = args.assembler.split(',')
        destination_file = ""
        for assembler in input_assembler_list:
            if assembler not in assembler_list:
                print(f'\nInvalid assembler name: {assembler}\n')
                exit(1)
            elif assembler == 'a5':
                launch_a5(sample, job_dir, fq_list, args.force)
            elif assembler == 'spades' or assembler == 'plasmidspades':
                trimming_output_file_found = check_trimming(job_dir)
                trimmer_dir = ""
                if not trimming_output_file_found or args.force:
                    if trimmer in trimmer_list:
                        trimmer_dir = launch_trimming(job_dir, fq_list, subset_size, trimmer)
                    else:
                        print(f'\nTrimmer {trimmer} not found\n')
                        exit(1)
                launch_spades(assembler, sample, job_dir, fq_dir, args.force, trimmer_dir, temp_dir)
            elif assembler == 'SKESA':
                launch_skesa(sample, job_dir, fq_list, args.force)
            elif assembler == 'shovill-spades':
                launch_shovill_spades(sample, job_dir, fq_list, args.force, temp_dir)
            elif assembler == 'shovill-SKESA':
                launch_shovill_skesa(sample, job_dir, fq_list, args.force, temp_dir)
            elif assembler == 'shovill-velvet':
                launch_shovill_velvet(sample, job_dir, fq_list, args.force, temp_dir)
            elif assembler == 'shovill-megahit':
                launch_shovill_megahit(sample, job_dir, fq_list, args.force, temp_dir)
        # Make Assembly file with filtering quality
        destination_file = select_assembly(job_dir, sample, int(args.minSize), input_assembler_list)

        # Bam file
        if args.Bam:
            bam_file = launch_fasta2bam(destination_file, fq_list)
            if not args.nobamFilter:
                launch_bam2stats(destination_file, bam_file, args.mbamDepth, args.mbamSize, args.mbamBasq,
                                 args.mbamMapq)

        # Make differentiation between chrom vs plasmid
        if args.plasFlow:
            threshold = 0.7
            outfile_prep = os.path.join(out_dir, sample, "filtered_plasflow.fasta")
            outfile_plasflow = os.path.join(out_dir, sample, "plasflow_predictions")
            launch_plasflow(destination_file, outfile_prep, outfile_plasflow, threshold)


def version():
    return "1.0"


def run():
    global usage

    usage = "seqAssembler.py [-fq fastq reads directory] [-o output directory] [-s sample file] "

    parser = argparse.ArgumentParser(
        prog='seqAssembler',
        usage=usage,
        description='SeqAssembler: pipeline CNR Resistance for Assembling - Version ' + version(),
    )

    parser.add_argument('-fq', '--fqDir', dest="fqDir", help='The directory containing fastq or fastq.gz files')
    parser.add_argument('-o', '--outDir', dest="outDir", default='./', help="The output directory name")
    parser.add_argument('-s', '--sampleFile', dest="sampleFile",
                        help="The sample file with sample names (default: $outDir/sample.csv)")
    parser.add_argument('-tr', '--trimmer', dest="trimmer", default='sickle',
                        help="sickle or trimmomatic or nothing (default: sickle)")
    parser.add_argument('-a', '--assembler', dest="assembler", default="a5",
                        help="Assembler names [a5,spades,plasmidspades, SKESA, 'shovill-spades', 'shovill-SKESA',"
                             " 'shovill-velvet', 'shovill-megahit'] as a comma separated list (default: a5)")
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
    parser.add_argument('-tmp', '--tempDir', dest="tempDir",
                        help="The tmp directory. (Not set by default)")
    parser.add_argument('-F', '--force', dest="force", action='store_true', default=False,
                        help="Force file overwrite (default: False)")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
