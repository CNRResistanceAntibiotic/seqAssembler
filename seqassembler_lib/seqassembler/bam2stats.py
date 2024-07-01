# !/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import logging
import os
from collections import OrderedDict

import matplotlib.pyplot as plt
import pandas as pd
import pysam
from Bio import SeqIO
from statistics import mean


def read_fasta_file(fas_file):
    with open(fas_file) as in_f:
        contigs = list(SeqIO.parse(in_f, 'fasta'))
    return contigs


def extract_bam_stats(bam_file, fas_file, out_dir, ext_report, plt_report, force=False):

    # extract contig list
    contigs = read_fasta_file(fas_file)

    # Open bam file
    bam = pysam.AlignmentFile(bam_file)
    bam_header_SQ = bam.header['SQ']

    tailles_contigs = {}
    for contig in bam_header_SQ:
        nom_contig = contig['SN']
        taille_contig = contig['LN']
        tailles_contigs[nom_contig] = taille_contig

    # Extract data for each contigs
    out_file = os.path.join(out_dir, 'assembly_stats.tsv')
    if not os.path.exists(out_file) or force:
        df = pd.DataFrame()
        for ctg in contigs:
            logging.info(f'Bam data for {ctg.id} ({len(ctg)}-bp) in process...')
            dt = {'Depth': round(mean(bam.count_coverage(ctg.id)[0]), 2), 'ctg': [ctg.id],
                  "Size": bam.get_reference_length(ctg.id)}
            total_mapq = 0
            count = 0
            for read in bam.fetch(ctg.id):
                # Ajouter la valeur MAPQ au total et incrémenter le compteur
                total_mapq += read.mapping_quality
                count += 1
            if count > 0:
                mapq_moyenne = total_mapq / count
            else:
                mapq_moyenne = 0
            dt['Mapq'] = round(mapq_moyenne, 2)
            print(dt)
            df = pd.concat([df, pd.DataFrame.from_dict(dt)])

        """
        for ctg in contigs:
            logging.info(f'Bam data for {ctg.id} ({len(ctg)}-bp) in process...')

            
            # Extract depth stats
            data = pysamstats.load_variation(bam, truncate=True, pad=True, max_depth=400, fafile=fas_file,
                                             chrom=ctg.id, start=1, end=len(ctg.seq))
            positions = data.pos
            dt = {'Depth': data.reads_all, 'Match_depth': data.matches, 'Seq': data.ref,
                  'ctg': [ctg.id] * len(positions)}
            
            
            if ext_report:
                # Depth contig report as tsv file
                cv_dic = OrderedDict([('Total_reads_depth', data.reads_all), ('Paired_reads_depth', data.reads_pp),
                                     ('Match_depth', data.matches), ('Mismatch_depth', data.mismatches),
                                     ('A_depth', data['A']),
                                     ('T_depth', data['T']), ('C_depth', data['C']), ('G_depth', data['G']),
                                     ('N_depth', data['N']),
                                     ('Deletion_depth', data.deletions), ('Insertion_depth', data.insertions)])
                df_1 = pd.DataFrame(cv_dic, index=positions)
                df_1.index.name = 'Positions'
                out_prefix = os.path.join(out_dir, '{0}_bam_stats'.format(ctg.id))
                df_1.to_csv(out_prefix + '.tsv', sep='\t', index=True)

                if plt_report:
                    # Depth contig report as png plot
                    subplots = {'reads_depth': {'data': ['Total_reads_depth', 'Paired_reads_depth'],
                                                'param': {'loc': 221, 'ylim': ''}},
                                'match_depth': {'data': ['Match_depth', 'Mismatch_depth'],
                                                'param': {'loc': 222, 'ylim': ''}},
                                'ambiguous_depth': {'data': ['N_depth'], 'param': {'loc': 223, 'ylim': ''}},
                                'indel_depth': {'data': ['Deletion_depth', 'Insertion_depth'],
                                                'param': {'loc': 224, 'ylim': ''}}}
                    plt.figure(figsize=(20, 10))
                    for plot in subplots:
                        for key in subplots[plot]['data']:
                            loc = subplots[plot]['param']['loc']
                            plt.subplot(loc)
                            values = cv_dic[key]
                            plt.plot(positions, values, linewidth=0.5, linestyle='-', label=key, alpha=0.7)
                        plt.legend(fontsize=12, loc=1)
                        plt.tick_params(axis='both', which='major', labelsize=12)
                        plt.tick_params(axis='both', which='minor', labelsize=12)
                        plt.ylabel('Depth', fontsize=12)
                        plt.xlabel('Positions', fontsize=12)
                    plt.savefig(out_prefix + '.png')
                    plt.close()

            # Extract base quality stats
            data = pysamstats.load_baseq_ext(bam, truncate=True, pad=True, max_depth=400, fafile=fas_file,
                                             chrom=ctg.id, start=1, end=len(ctg.seq))
            dt['Basq'] = data.rms_baseq
            dt['Match_basq'] = data.rms_baseq_matches

            if ext_report:
                # Base quality report as tsv file
                bqDic = OrderedDict(
                    [('RMS_base_quality', data.rms_baseq), ('RMS_match_base_quality', data.rms_baseq_matches),
                     ('RMS_mismatch_base_quality', data.rms_baseq_mismatches)])
                df_2 = pd.DataFrame(bqDic, index=positions)
                df_2.index.name = 'Positions'
                out_prefix = os.path.join(out_dir, f'{ctg.id}_assembly_qual')
                df_2.to_csv(out_prefix + '.tsv', sep='\t', index=True)

                if plt_report:
                    # Base quality report as png plot
                    plt.figure(figsize=(20, 5))
                    for n, key in enumerate(bqDic):
                        values = bqDic[key]
                        plt.subplot(1, 3, n + 1)
                        plt.plot(positions, values, linewidth=0.5, linestyle='-', label=key)
                        plt.legend(fontsize=12, loc=1)
                        # plt.ylim(0, 500)
                        plt.tick_params(axis='both', which='major', labelsize=12)
                        plt.tick_params(axis='both', which='minor', labelsize=12)
                        plt.ylabel('Phred quality', fontsize=12)
                        plt.xlabel('Positions', fontsize=12)
                    plt.savefig(out_prefix + '.png')
                    plt.close()
            
            # Extract mapping quality data
            data = pysamstats.load_mapq(bam, truncate=True, pad=True, max_depth=400, fafile=fas_file,
                                        chrom=ctg.id, start=1, end=len(ctg.seq))
            dt['Mapq'] = data.rms_mapq
            df = pd.concat([df, pd.DataFrame.from_dict(dt)])

            if ext_report:
                # Mapping quality report as tsv file
                mq_dic = OrderedDict([('RMS_mapq', data.rms_mapq), ('MAX_mapq', data.max_mapq),
                                     ('Nbr_mapqO', data.reads_mapq0)])
                df_3 = pd.DataFrame(mq_dic, index=positions)
                df_3.index.name = 'Positions'
                out_prefix = os.path.join(out_dir, f'{ctg.id}_assembly_mapq')
                df_3.to_csv(out_prefix + '.tsv', sep='\t', index=True)

                if plt_report:
                    # Mapping quality report as png plot
                    plt.figure(figsize=(20, 20))
                    for n, key in enumerate(mq_dic):
                        values = mq_dic[key]
                        plt.subplot(3, 1, n + 1)
                        plt.plot(positions, values, linewidth=0.5, linestyle='-', label=key)
                        plt.legend(fontsize=12, loc=1)
                        # plt.ylim(0, 500)
                        plt.tick_params(axis='both', which='major', labelsize=12)
                        plt.tick_params(axis='both', which='minor', labelsize=12)
                        if n < 2:
                            plt.ylabel('Mapping quality', fontsize=12)
                        else:
                            plt.ylabel('Number of reads with mapping quality = 0', fontsize=12)
                        plt.xlabel('Positions', fontsize=12)
                    plt.savefig(out_prefix + '.png')
                    plt.close()
        """
        logging.info(f'\nWrite the main results in {out_file}:')
        # df.boxplot(by='ctg', column=['Depth', 'Match_depth', 'Basq', 'Match_basq', 'Mapq'])
        # plt.savefig(os.path.splitext(outfile)[0]+'.png')
        results = []
        for ctg in [x.id for x in contigs] + ['overall']:
            logging.info(f'{ctg} in process...')
            res_dic = OrderedDict([('ID', ctg)])
            if ctg == 'overall':
                values = df
                data = values.describe(percentiles=[0.10, 0.90])
            else:
                values = df[df['ctg'] == ctg]
                data = values.describe(percentiles=[0.10, 0.50, 0.90])

            for i in ['Depth', 'Mapq']:
                N20 = round(100 * values[values[i] >= 20].index.size / float(values.index.size), 2)
                res_dic[f'Perc_{i}_>=20'] = N20
                N30 = round(100 * values[values[i] >= 30].index.size / float(values.index.size), 2)
                res_dic[f'Perc_{i}_>=30'] = N30
                for j in ["mean", "std", "max", "min"]:
                    key = f'{i}_{j}'
                    value = data.loc[j, i].round(2)
                    res_dic[key] = value
                for j in ["10%", "50%", "90%"]:
                    key = f'{i}_{j}_percentile'
                    value = data.loc[j, i].round(2)
                    res_dic[key] = value
            results.append(res_dic)
        df = pd.DataFrame(results)
        df.to_csv(out_file, sep='\t', index=False)
    else:
        logging.info('\nThe main output file already done!\n')

    return out_file, contigs


def filter_contigs(result_file, contigs, m_size, m_basq, m_mapq, m_depth, rename):
    # Load the main results
    results = []
    with open(result_file) as in_f:
        header = ""
        for n, line in enumerate(in_f):
            if n == 0:
                header = line.strip().split('\t')
            else:
                line = line.strip().split('\t')
                results.append(dict(zip(header, line)))

    # Extract outfiltered contigs
    print('\nStart filtering:')
    del_IDs = []
    for data in results:
        ID = data['ID']
        print(f'Filtering {ID}')
        if float(data['Depth_mean']) <= m_depth:
            del_IDs.append(ID)
        if float(data['Mapq_mean']) <= m_mapq:
            del_IDs.append(ID)
        """
        if float(data['Basq_mean']) <= m_basq:
            del_IDs.append(ID)
        """
        if float(data['Size']) <= m_size:
            del_IDs.append(ID)

    del_IDs = list(set(del_IDs))
    for ID in del_IDs + ['overall']:
        for n, data in enumerate(results):
            if ID == data['ID']:
                del results[n]
                break

    records = []
    n = 0
    for ctg in contigs:
        if ctg.id not in del_IDs:
            n += 1
            if rename:
                ctg.id = f'ctg_{n}'
                ctg.description = ''
            records.append(ctg)
    logging.info(f'\n{len(del_IDs)} deleted contigs')
    logging.info(f'{len(records)} remaining contigs\n')

    # Write the filtered contigs
    out_file = os.path.join(os.path.dirname(result_file), 'assembly_filtered.fas')
    with open(out_file, 'w') as out_f:
        SeqIO.write(records, out_f, 'fasta')

    # Assembly stats after filtering:
    ds = pd.DataFrame(results)
    nw = {}

    for item in header:
        if item == 'ID':
            nw[item] = 'overall'
        elif 'max' in item:
            nw[item] = ds[item].astype(float).max()
        elif 'min' in item:
            nw[item] = ds[item].astype(float).min()
        elif item == 'Size' or item == 'Nbr_ambiguous':
            nw[item] = ds[item].astype(float).sum()
        else:
            d_item = ds[item].astype(float) * ds['Size'].astype(int)
            nw[item] = round(ds[item].astype(float).sum() / len(records), 2)

    results.append(nw)

    df = pd.DataFrame(results)
    df = df[header]
    out_file = os.path.splitext(result_file)[0] + '_filtered.tsv'
    df.to_csv(out_file, sep='\t', index=False)


def pre_main(args):
    fas_file = args.fasFile
    bam_file = args.bamFile
    ext_report = args.extReport
    plt_report = args.pltReport
    force = args.force
    filter_bam = args.filter
    out_dir = args.outDir
    m_size = args.mSize
    m_basq = args.mBasq
    m_mapq = args.mMapq
    m_depth = args.mDepth
    rename = args.rename

    main(fas_file, bam_file, ext_report, plt_report, force, filter_bam, out_dir, m_size, m_basq, m_mapq, m_depth,
         rename)


def main(fas_file="", bam_file="", ext_report=False, plt_report=False, force=True, filter_bam=True, out_dir="",
         m_size=500, m_basq=20, m_mapq=30, m_depth=20, rename=True):

    if out_dir == '':
        out_dir = os.path.dirname(bam_file)

    out_dir_path = os.path.join(out_dir, 'bam_stats')
    if not os.path.exists(out_dir_path):
        os.mkdir(out_dir_path)

    result_file, contigs = extract_bam_stats(bam_file, fas_file, out_dir_path, ext_report, plt_report, force)
    if filter_bam:
        filter_contigs(result_file, contigs, m_size, m_basq, m_mapq, m_depth, rename)


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='bam2stats - Version ' + version())
    parser.add_argument('-f', '--fasFile', dest="fasFile",
                        default='/media/bacteriologie/TX/NGS-caen-enterobacter/test/CNR1717/CNR1717.fasta',
                        help='Reference fasta file')
    parser.add_argument('-b', '--bamFile', dest="bamFile",
                        default='/media/bacteriologie/TX/NGS-caen-enterobacter/test/CNR1717/CNR1717.bam',
                        help='Bam file')
    parser.add_argument('-o', '--outDir', dest="outDir", default='',
                        help="Output directory name (default: bam file directory)")
    parser.add_argument('-ext', '--extReport', dest="extReport", action='store_true', default=False,
                        help="Make an extended report (default: False)")
    parser.add_argument('-plt', '--pltReport', dest="pltReport", action='store_true', default=False,
                        help="Make plots (default: False)")
    parser.add_argument('-nfl', '--no_filter', dest="filter", action='store_false', default=True,
                        help="Do not filter fasta file (default: False)")
    parser.add_argument('-nrn', '--no_rename', dest="rename", action='store_false', default=True,
                        help="Do not rename contigs during filtering (default: False)")
    parser.add_argument('-md', '--mDepth', dest="mDepth", default=20,
                        help="Mean depth threshold for the scaffolds (default: 20)")
    parser.add_argument('-ms', '--mSize', dest="mSize", default=500, help="Minimum size of scaffolds (default: 500)")
    parser.add_argument('-mbq', '--mBasq', dest="mBasq", default=20,
                        help="Mean base quality threshold for the scaffolds (default: 20)")
    parser.add_argument('-mmq', '--mMapq', dest="mMapq", default=30,
                        help="Mean mapping quality threshold for the scaffolds (default: 30)")
    parser.add_argument('-F', '--force', dest="force", action='store_true', default=True,
                        help="Force file overwrite (default: False)")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
