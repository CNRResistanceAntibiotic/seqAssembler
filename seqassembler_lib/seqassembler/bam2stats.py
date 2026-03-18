# !/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import os
from collections import OrderedDict, defaultdict

import matplotlib.pyplot as plt
import pandas as pd
import pysam
from Bio import SeqIO
from statistics import mean
import numpy as np


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
            print(f'Bam data for {ctg.id} ({len(ctg)}-bp) in process...')
            coverage_position_list = []
            ctg_list = []
            size_list = []
            x = 0
            for base_coverage_list in bam.count_coverage(ctg.id):
                for i, depth in enumerate(base_coverage_list):
                    if x != 0:
                        tmp_depth = coverage_position_list[i]
                        coverage_position_list[i] = tmp_depth + depth
                    else:
                        coverage_position_list.append(depth)
                        ctg_list.append(ctg.id)
                        if size_list:
                            size_list.append(0)
                        else:
                            size_list.append(bam.get_reference_length(ctg.id))
                x += 1

            dt = {'Depth': coverage_position_list, 'ctg': ctg_list, "Size": size_list}
            count = 0
            mapq_by_position = defaultdict(list)
            for read in bam.fetch(ctg.id):
                if not read.is_unmapped:
                    start = read.reference_start
                    end = read.reference_end
                    # Parcourir chaque position couverte par la lecture
                    for pos in range(start, end):
                        mapq_by_position[pos].append(read.mapping_quality)
                count += 1
            mapq_list = []
            for pos in range(0, len(dt["Depth"])):
                if pos in mapq_by_position:
                    qmap_l = mapq_by_position[pos]
                    mapq_list.append(round(mean(qmap_l), 2))
                else:
                    mapq_list.append(0)
            dt['Mapq'] = mapq_list
            df = pd.concat([df, pd.DataFrame.from_dict(dt)])

        print(f'\nWrite the main results in {out_file}:')
        # df.boxplot(by='ctg', column=['Depth', 'Match_depth', 'Basq', 'Match_basq', 'Mapq'])
        # plt.savefig(os.path.splitext(outfile)[0]+'.png')
        results = []
        for ctg in [x.id for x in contigs] + ['overall']:
            print(f'{ctg} in process...')
            res_dic = OrderedDict([('ID', ctg)])
            if ctg == 'overall':
                values = df
                res_dic['Size'] = values['Size'].sum()
            else:
                values = df[df['ctg'] == ctg]
                res_dic['Size'] = df[df['ctg'] == ctg]["Size"].values[0]
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
        df_result = pd.DataFrame(results)
        df_result.replace(np.nan, 0, inplace=True)
        df_result.to_csv(out_file, sep='\t', index=False)
    else:
        print('\nThe main output file already done!\n')

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
    print(f'\n{len(del_IDs)} deleted contigs')
    print(f'{len(records)} remaining contigs\n')

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
