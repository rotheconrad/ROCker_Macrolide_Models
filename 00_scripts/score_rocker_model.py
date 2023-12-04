#!/usr/bin/env python

'''Score the rocker model against static filters

Takes inputs:

1) original mock metagenome fasta
2) best hit filtered hmmsearch result (besthit_filter_hmm.py)
3) best hit filtered blastx alignment (besthit_filter_blast.py)
4) best hit filtered passing rocker alignments
5) best hit filtered failing rocker alignments

* expects inputs to be best hit filtered.
* newest rocker should have built in best hit filter.

Outputs:

1) tsv data table
2) pdf bar plot

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: April 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def read_fasta(fp):
    ''' parses a fasta file into name, seq '''
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def parse_mm(mm, outdir):
    ''' Read fasta and counts positive labeled reads. returns count '''

    labels = {}
    pos_reads = {}
    neg_reads = {}
    pfo = open(f'{outdir}/all_reads-pos.fasta', 'w')
    nfo = open(f'{outdir}/all_reads-neg.fasta', 'w')

    with open(mm, 'r') as file:
        for name, seq in read_fasta(file):
            name = name[1:]
            label = name.split(';')[-1]
            labels[label] = ''
            if label == 'Positive':
                entry = f'>{name}\n{seq}\n'
                pos_reads[name] = entry
                pfo.write(entry)

            elif label == 'Non_Target' or label == 'Negative':
                entry = f'>{name}\n{seq}\n'
                neg_reads[name] = entry
                nfo.write(entry)

            else:
                print(name, label)

    print('Labels in this dataset:', list(labels.keys()))

    pfo.close()
    nfo.close()

    return pos_reads, neg_reads


def score_results(TP, FP, TN, FN):

    P = TP + FP
    N = TN + FN
    Total = P + N
    FNR = FN/(TP + FN)
    FPR = FP/P if P > 0 else 0
    Sensitivity = 1 - FNR
    Specificity = 1 - FPR
    Accuracy = (TP + TN) / (P + N)
    Precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    Recall = TP / (TP + FN)
    if (Precision + Recall) > 0:
        F1 = 2 * (Precision * Recall) / (Precision + Recall)
    else: F1 = 0

    results = [
                Total, P, N, TP, FP, TN, FN, FNR, FPR,
                Sensitivity, Specificity, Accuracy,
                Precision, Recall, F1
                ]

    return results


def score_blastx(bx, pos_reads, neg_reads, outdir):

    pos_reads = pos_reads.copy()
    neg_reads = neg_reads.copy()

    customA = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    customB = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    e30 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    e20 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    e10 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    e03 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}

    pos_total = len(pos_reads)
    pos_count = 0

    fasta = {
                'cA-TP': [], 'cA-FP': [], 'cA-FN': [], 'cA-TN': [],
                'cB-TP': [], 'cB-FP': [], 'cB-FN': [], 'cB-TN': [],
                'e30-TP': [], 'e30-FP': [], 'e30-FN': [], 'e30-TN': [],
                'e20-TP': [], 'e20-FP': [], 'e20-FN': [], 'e20-TN': [],
                'e10-TP': [], 'e10-FP': [], 'e10-FN': [], 'e10-TN': [],
                'e03-TP': [], 'e03-FP': [], 'e03-FN': [], 'e03-TN': [],
                }

    with open(bx, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]
            evalue = float(X[10]) # evalue
            pid = float(X[2]) # percent sequence identity
            aLen = int(X[3]) # read alignment length
            qLen = int(X[12]) / 3 # full length of read divided by 3 amino acids
            pMatch = aLen / qLen # percent match length of read length

            # count positive reads in blastX alignment
            if label == 'Positive':
                pos_count += 1
                pf = pos_reads.pop(query)
            else:
                nf = neg_reads.pop(query)
            ###############################
            # count customA filter ########
            if pMatch >= 0.7 and pid >= 95:
                if label == 'Positive':
                    customA['TP'] += 1
                    fasta['cA-TP'].append(pf)

                else:
                    customA['FP'] += 1
                    fasta['cA-FP'].append(nf)

            else:
                if label == 'Positive':
                    customA['FN'] += 1
                    fasta['cA-FN'].append(pf)

                else:
                    customA['TN'] += 1
                    fasta['cA-TN'].append(nf)
            ###############################
            # count customB filter ########
            if pMatch >= 0.5 and pid >= 90:
                if label == 'Positive':
                    customB['TP'] += 1
                    fasta['cB-TP'].append(pf)

                else:
                    customB['FP'] += 1
                    fasta['cB-FP'].append(nf)

            else:
                if label == 'Positive':
                    customB['FN'] += 1
                    fasta['cB-FN'].append(pf)

                else:
                    customB['TN'] += 1
                    fasta['cB-TN'].append(nf)
            ###############################
            # count e30 filter ############
            if evalue <= 1e-30:
                if label == 'Positive':
                    e30['TP'] += 1
                    fasta['e30-TP'].append(pf)

                else:
                    e30['FP'] += 1
                    fasta['e30-FP'].append(nf)

            else:
                if label == 'Positive':
                    e30['FN'] += 1
                    fasta['e30-FN'].append(pf)

                else:
                    e30['TN'] += 1
                    fasta['e30-TN'].append(nf)
            ###############################
            # count e20 filter ############
            if evalue <= 1e-20:
                if label == 'Positive':
                    e20['TP'] += 1
                    fasta['e20-TP'].append(pf)

                else:
                    e20['FP'] += 1
                    fasta['e20-FP'].append(nf)

            else:
                if label == 'Positive':
                    e20['FN'] += 1
                    fasta['e20-FN'].append(pf)

                else:
                    e20['TN'] += 1
                    fasta['e20-TN'].append(nf)
            ###############################=
            # count e10 filter ############
            if evalue <= 1e-10:
                if label == 'Positive':
                    e10['TP'] += 1
                    fasta['e10-TP'].append(pf)

                else:
                    e10['FP'] += 1
                    fasta['e10-FP'].append(nf)

            else:
                if label == 'Positive':
                    e10['FN'] += 1
                    fasta['e10-FN'].append(pf)

                else:
                    e10['TN'] += 1
                    fasta['e10-TN'].append(nf)
            ###############################
            # count e03 filter ############
            if evalue <= 1e-03:
                if label == 'Positive':
                    e03['TP'] += 1
                    fasta['e03-TP'].append(pf)

                else:
                    e03['FP'] += 1
                    fasta['e03-FP'].append(nf)

            else:
                if label == 'Positive':
                    e03['FN'] += 1
                    fasta['e03-FN'].append(pf)

                else:
                    e03['TN'] += 1
                    fasta['e03-TN'].append(nf)
            ###############################

    # group the data to return for ease
    data = [customA, customB, e30, e20, e10, e03]

    # check all positive reads are counted
    print('\n\nPositive read count from metagenome:', pos_total)
    print('\nPositive read count from Diamond/BlastX alignment:', pos_count)
    pos_diff = pos_total - pos_count # positive reads not aligned
    print('Positive reads not aligned by BlastX:', pos_diff)

    # write the fasta files for each filter and category:
    for val, fastas in fasta.items():
        with open(f'{outdir}/{val}.fasta', 'w') as o:
            for entry in fastas:
                o.write(entry)

    # write fasta files not aligned by blastx
    with open(f'{outdir}/blastx_failed-pos.fasta', 'w') as o:
        for val, fasta in pos_reads.items():
            o.write(fasta)

    with open(f'{outdir}/blastx_failed-neg.fasta', 'w') as o:
        for val, fasta in neg_reads.items():
            o.write(fasta)

    return data, pos_diff


def score_rocker(rp, rf, pos_reads, neg_reads, outdir):

    pos_reads = pos_reads.copy()
    neg_reads = neg_reads.copy()

    roc = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}

    passing_reads = {}

    pos_total = len(pos_reads)
    pos_count = 0

    fasta = {'roc-TP': [], 'roc-FP': [], 'roc-FN': [], 'roc-TN': []}

    # rocker passing TP, FP
    with open(rp, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]
            evalue = float(X[10]) # evalue
            passing_reads[query] = ''

            if label == 'Positive':
                roc['TP'] += 1
                pos_count += 1
                pf = pos_reads.pop(query)
                fasta['roc-TP'].append(pf)
            else:
                roc['FP'] += 1
                nf = neg_reads.pop(query)
                fasta['roc-FP'].append(nf)

    # rocker failing FN, TN
    with open(rf, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]
            evalue = float(X[10]) # evalue

            # don't count the same read twice
            if query in passing_reads: continue

            if label == 'Positive':
                roc['FN'] += 1
                pos_count += 1
                pf = pos_reads.pop(query)
                fasta['roc-FN'].append(pf)
            else:
                roc['TN'] += 1
                nf = neg_reads.pop(query)
                fasta['roc-TN'].append(nf)

    # check all positive reads are counted
    print('\nPositive read count from ROCkOut:', pos_count)
    pos_diff = pos_total - pos_count
    print('Positive reads not reported by ROCkOut:', pos_diff)

    # write the fasta files for each filter and category:
    for val, fastas in fasta.items():
        with open(f'{outdir}/{val}.fasta', 'w') as o:
            for entry in fastas:
                o.write(entry)

    # write fasta files not aligned by blastx
    with open(f'{outdir}/rocker_failed-pos.fasta', 'w') as o:
        for val, fasta in pos_reads.items():
            o.write(fasta)

    with open(f'{outdir}/rocker_failed-neg.fasta', 'w') as o:
        for val, fasta in neg_reads.items():
            o.write(fasta)

    return roc


def score_hmm(hm, pos_reads, neg_reads, outdir):

    pos_reads = pos_reads.copy()
    neg_reads = neg_reads.copy()

    hDf = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    h30 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    h20 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    h10 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    h03 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}

    pos_total = len(pos_reads)
    pos_count = 0

    fasta = {
                'Hdf-TP': [], 'Hdf-FP': [], 'Hdf-FN': [], 'Hdf-TN': [],
                'He30-TP': [], 'He30-FP': [], 'He30-FN': [], 'He30-TN': [],
                'He20-TP': [], 'He20-FP': [], 'He20-FN': [], 'He20-TN': [],
                'He10-TP': [], 'He10-FP': [], 'He10-FN': [], 'He10-TN': [],
                'He03-TP': [], 'He03-FP': [], 'He03-FN': [], 'He03-TN': [],
                }

    with open(hm, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]
            evalue = float(X[1]) # evalue

            # Count positive reads mapped by hmmsearch
            if label == 'Positive':
                pos_count += 1
                pf = pos_reads.pop(query)
            else:
                nf = neg_reads.pop(query)
            ###############################
            # count hmm default ###########
            if evalue <= 0.01:
                if label == 'Positive':
                    hDf['TP'] += 1
                    fasta['Hdf-TP'].append(pf)
                else:
                    hDf['FP'] += 1
                    fasta['Hdf-FP'].append(nf)
            else:
                if label == 'Positive':
                    hDf['FN'] += 1
                    fasta['Hdf-FN'].append(pf)
                else:
                    hDf['TN'] += 1
                    fasta['Hdf-TN'].append(nf)
            ###############################
            # count h30 filter ############
            if evalue <= 1e-30:
                if label == 'Positive':
                    h30['TP'] += 1
                    fasta['He30-TP'].append(pf)
                else:
                    h30['FP'] += 1
                    fasta['He30-FP'].append(nf)
            else:
                if label == 'Positive':
                    h30['FN'] += 1
                    fasta['He30-FN'].append(pf)
                else:
                    h30['TN'] += 1
                    fasta['He30-TN'].append(nf)
            ###############################
            # count h20 filter ############
            if evalue <= 1e-20:
                if label == 'Positive':
                    h20['TP'] += 1
                    fasta['He20-TP'].append(pf)
                else:
                    h20['FP'] += 1
                    fasta['He20-FP'].append(nf)
            else:
                if label == 'Positive':
                    h20['FN'] += 1
                    fasta['He20-FN'].append(pf)
                else:
                    h20['TN'] += 1
                    fasta['He20-TN'].append(nf)
            ###############################
            # count h10 filter ############
            if evalue <= 1e-10:
                if label == 'Positive':
                    h10['TP'] += 1
                    fasta['He10-TP'].append(pf)
                else:
                    h10['FP'] += 1
                    fasta['He10-FP'].append(nf)
            else:
                if label == 'Positive':
                    h10['FN'] += 1
                    fasta['He10-FN'].append(pf)
                else:
                    h10['TN'] += 1
                    fasta['He10-TN'].append(nf)
            ###############################
            # count h03 filter ############
            if evalue <= 1e-03:
                if label == 'Positive':
                    h03['TP'] += 1
                    fasta['He03-TP'].append(pf)
                else:
                    h03['FP'] += 1
                    fasta['He03-FP'].append(nf)
            else:
                if label == 'Positive':
                    h03['FN'] += 1
                    fasta['He03-FN'].append(pf)
                else:
                    h03['TN'] += 1
                    fasta['He03-TN'].append(nf)
            ###############################

    # group the data for ease
    data = [h30, h20, h10, h03, hDf]

    # check all positive reads are counted
    print('\nPositive read count from hmmsearch:', pos_count)
    pos_diff = pos_total - pos_count # positive reads not aligned
    print('Positive reads not aligned by hmmsearch:', pos_diff)

    # write the fasta files for each filter and category:
    for val, fastas in fasta.items():
        with open(f'{outdir}/{val}.fasta', 'w') as o:
            for entry in fastas:
                o.write(entry)

    # write fasta files not aligned by blastx
    with open(f'{outdir}/hmm_failed-pos.fasta', 'w') as o:
        for val, fasta in pos_reads.items():
            o.write(fasta)

    with open(f'{outdir}/hmm_failed-neg.fasta', 'w') as o:
        for val, fasta in neg_reads.items():
            o.write(fasta)

    return data, pos_diff


def build_bar_plots(df, out):

    fs = 12 # set font size
    metrics = ['F1', 'FPR', 'FNR']
    labels = df.columns.to_list()

    fig, axes = plt.subplots(len(metrics), 1, figsize=(4.25, 10), sharex=True)

    for i, met in enumerate(metrics):
        ax = axes[i]
        data = df.T[met].to_list()
        _ = ax.bar(labels, data, color='#636363', width=0.5, alpha=0.75)
        ax.set_ylabel(met, fontsize=fs)
        ax.set_ylim([0, 1])
        ax.vlines(x=0.5, ymin=0, ymax=1, color='#000000', ls='--', lw=1.5,)
        ax.vlines(x=6.5, ymin=0, ymax=1, color='#000000', ls='--', lw=1.5,)
        ax.yaxis.grid(
            which="major", color='#d9d9d9', linestyle='--',
            linewidth=1, zorder=1
            )
        ax.set_axisbelow(True)

    axes[0].text(3.5, 1.02, 'BLASTx', fontweight='black', ha='center')
    axes[0].text(9.2, 1.02, 'HMM', fontweight='black', ha='center')

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

    fig.set_tight_layout(True)
    plt.savefig(f'{out}_bar_plot.pdf')
    plt.close()

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-mm', '--mock_metagenome_fasta',
        help='Original mock metagenome fasta file.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-bx', '--blastx_alignments_file',
        help='Best hit filtered blastx alignment file.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-hm', '--filtered_hmmer_file',
        help='best hit filtered hmmsearch results.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-rp', '--rocker_passing_alignments',
        help='BlastX alignments passing the rocker filter.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-rf', '--rocker_failing_alignments',
        help='BlastX alignments failing the rocker filter.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Prefix to use for output files.',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # define input params
    mm = args['mock_metagenome_fasta'] # original mock metagenome fasta
    bx = args['blastx_alignments_file'] # best hit filtered blastx file
    hm = args['filtered_hmmer_file'] # best hit filtered hmmsearch file
    rp = args['rocker_passing_alignments'] # rocker passing alignments
    rf = args['rocker_failing_alignments'] # rocker failing alignments
    out = args['output_prefix']  # output file prefix

    # Do what you came here to do:
    print('\n\nRunning Script...\n\n')

    # create directory for fasta file outputs
    outdir = '/'.join(out.split('/')[:-1]) + '/fastas'
    Path(outdir).mkdir(parents=True, exist_ok=True)
    # Get all reads labeled as positive
    pos_reads, neg_reads = parse_mm(mm, outdir)

    ####
    # PHASE 1: Count TP, FP, TN, FN
    ###############################

    # BLASTx
    bx_data, bx_pos_fail = score_blastx(bx, pos_reads, neg_reads, outdir)
    bx_results = []
    for d in bx_data:
        r = score_results(d['TP'], d['FP'], d['TN'], d['FN'])
        bx_results.append(r)

    # ROCker
    # Count TP, FP, TN, FN
    roc_data = score_rocker(rp, rf, pos_reads, neg_reads, outdir)
    roc_results = score_results(
                            roc_data['TP'], roc_data['FP'],
                            roc_data['TN'], roc_data['FN']
                            )

    # HMM
    # Count TP, FP, TN, FN
    hm_data, hm_pos_fail = score_hmm(hm, pos_reads, neg_reads, outdir)
    hm_results = []
    # hmm models tend to have more unmapped reads than blastx algos
    # we count the unmapped pos reads from hmm model as FN
    # subtract unmapped bx reads from hm reads
    pos_fail = hm_pos_fail - bx_pos_fail
    for d in hm_data:
        r = score_results(d['TP'], d['FP'], d['TN'], d['FN']+pos_fail)
        hm_results.append(r)

    ####
    # PHASE 02: output table and pdf
    ################################

    scores = {
                'ROCker': roc_results, 'Custom-A': bx_results[0],
                'Custom-B': bx_results[1], 'evalue 1e-30': bx_results[2],
                'evalue 1e-20': bx_results[3], 'evalue 1e-10': bx_results[4],
                'evalue 1e-3': bx_results[5], 'HMM 1e-30': hm_results[0],
                'HMM 1e-20': hm_results[1], 'HMM 1e-10': hm_results[2],
                'HMM 1e-3': hm_results[3], 'HMM default': hm_results[4],
                }
    rows = [
            'Total', 'P', 'N', 'TP', 'FP', 'TN', 'FN', 'FNR', 'FPR',
            'Sensitivity', 'Specificity', 'Accuracy',
            'Precision', 'Recall', 'F1'
            ]

    df = pd.DataFrame(scores, index=rows).round(2)
    df.to_csv(f'{out}_score_table.tsv', sep='\t')

    print('\n\n', df)

    _ = build_bar_plots(df, out)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')

if __name__ == "__main__":
    main()
