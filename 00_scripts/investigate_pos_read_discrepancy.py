#!/usr/bin/env python

'''Investigate the positive read discrepancy

rocker align aligns fewer reads than rocker labels as positive during
the read simulation step.

Possibly because a read is labeled as positive if only 1 bp overlaps
with the target gene region.

The purpose of this script is to verify this hypothesis.

It plots histograms for the simulated reads labeled as positive
that do not align by blastx or hmm or don't pass the rocker filter.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: July 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import pandas as pd
import matplotlib.pyplot as plt

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


def parse_blastx(bx):

    aligned_pos = {}

    with open(bx, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]

            # count positive reads in blastX alignment
            if label == 'Positive':
                aligned_pos[query] = ''

    return aligned_pos


def parse_rocker(rp, rf):

    roc_pos = {}

    # rocker passing TP, FP
    with open(rp, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]

            if label == 'Positive':
                roc_pos[query] = ''

    # rocker failing FN, TN
    with open(rf, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]

            if label == 'Positive':
                roc_pos[query] = ''

    return roc_pos


def parse_hmm(hm):

    hmm_pos = {}

    with open(hm, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]
            if label == 'Positive':
                hmm_pos[query] = ''

    return hmm_pos


def plot_hists(data1, data2, outfile):

    fs = 12 # set font size
    # Set the colors
    bar_color = '#2171b5'
    gridM = '#bdbdbd'
    gridm = '#d9d9d9'
    alpha = 0.6

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(7, 7), sharex=True)

    # Plot the data
    ax1.hist(data1['bp'], bins=30, rwidth=0.9, color=bar_color, alpha=alpha)
    ax2.hist(data1['pct'], bins=30, rwidth=0.9, color=bar_color, alpha=alpha)
    ax3.hist(data2['bp'], bins=30, rwidth=0.9, color=bar_color, alpha=alpha)
    ax4.hist(data2['pct'], bins=30, rwidth=0.9, color=bar_color, alpha=alpha)

    # Plot labels
    ax1.set_xlabel('', fontsize=fs)
    ax1.set_ylabel('Count', fontsize=fs)
    ax2.set_xlabel('', fontsize=fs)
    ax2.set_ylabel('', fontsize=fs)
    ax3.set_xlabel('Gene overlap of read (bp)', fontsize=fs)
    ax3.set_ylabel('Count', fontsize=fs)
    ax4.set_xlabel('Gene overlap of read (%)', fontsize=fs)
    ax4.set_ylabel('', fontsize=fs)

    ax1.text(
            0.02, 1.02, f"Positive aligned reads: {len(data1['names'])}",
            fontsize=12, color='r',
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax1.transAxes
            )
    ax3.text(
            0.02, 1.02, f"Positive un-aligned reads: {len(data2['names'])}",
            fontsize=12, color='r',
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax3.transAxes
            )

    for ax in [ax1, ax2, ax3, ax4]:
        # Set plot/grid style
        ax.minorticks_on()
        ax.tick_params(
            which='minor', axis='both', left=False, bottom=False
            )
        ax.tick_params(
                    which='major', axis='both',
                    left=True, bottom=True,
                    size=6, width=2, tickdir='inout',
                    labelsize=12, zorder=10
                    )
        ax.yaxis.grid(
            which="minor", color=gridm, linestyle='--',
            linewidth=1, alpha=0.6, zorder=1
            )
        ax.yaxis.grid(
            which="major", color=gridM, linestyle='--',
            linewidth=1.5, alpha=0.4, zorder=1
            )
        ax.set_axisbelow(True)
        for spine in ax.spines.values(): spine.set_linewidth(2)

    fig.set_tight_layout(True)
    plt.savefig(f'{outfile}.pdf')
    plt.close()

    with open(f'{outfile}.tsv', 'w') as ofile:
        for name, bp, pct in zip(data1['names'], data1['bp'], data1['pct']):
            ofile.write(f'{name}\t{bp}\t{pct}\tAligned\n')
        for name, bp, pct in zip(data2['names'], data2['bp'], data2['pct']):
            ofile.write(f'{name}\t{bp}\t{pct}\tNot Aligned\n')

    return True


def parse_mm(mm, aligned_pos, roc_pos, hmm_pos, out):

    pos_count = 0
    not_aligned = {'bp': [], 'pct': [], 'names': []}
    aligned = {'bp': [], 'pct': [], 'names': []}
    not_rocked = {'bp': [], 'pct': [], 'names': []}
    rocked = {'bp': [], 'pct': [], 'names': []}
    not_hmmed = {'bp': [], 'pct': [], 'names': []}
    hmmed = {'bp': [], 'pct': [], 'names': []}

    with open(mm, 'r') as file:
        for name, seq in read_fasta(file):
            name = name[1:]
            X = name.split(';')
            label = X[-1]
            overlapping_bp = int(X[3])
            pct_overlap = float(X[4])
            
            if label == 'Positive':
                pos_count += 1

                if name in aligned_pos:
                    aligned['bp'].append(overlapping_bp)
                    aligned['pct'].append(pct_overlap)
                    aligned['names'].append(name)

                elif name not in aligned_pos:
                    not_aligned['bp'].append(overlapping_bp)
                    not_aligned['pct'].append(pct_overlap)
                    not_aligned['names'].append(name)

                if name in roc_pos:
                    rocked['bp'].append(overlapping_bp)
                    rocked['pct'].append(pct_overlap)
                    rocked['names'].append(name)

                elif name not in roc_pos:
                    not_rocked['bp'].append(overlapping_bp)
                    not_rocked['pct'].append(pct_overlap)
                    not_rocked['names'].append(name)

                if name in hmm_pos:
                    hmmed['bp'].append(overlapping_bp)
                    hmmed['pct'].append(pct_overlap)
                    hmmed['names'].append(name)

                elif name not in hmm_pos:
                    not_hmmed['bp'].append(overlapping_bp)
                    not_hmmed['pct'].append(pct_overlap)
                    not_hmmed['names'].append(name)

    # build hist plots
    _ = plot_hists(aligned, not_aligned, f'{out}_aligned')
    _ = plot_hists(rocked, not_rocked, f'{out}_rocked')
    _ = plot_hists(hmmed, not_hmmed, f'{out}_hmmed')

    print(f'\n\nTotal Positives reads created: {pos_count}')
    print(f"Total reads not aligned: {len(not_aligned['bp'])}")
    print(f"Total reads not in rocker filter output: {len(not_rocked['bp'])}")
    print(f"Total reads not in hmm search results: {len(not_hmmed['bp'])}")

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

    # parse blastx alignment file for positive reads
    aligned_pos = parse_blastx(bx)

    # parse rocker filter files for positive reads
    roc_pos = parse_rocker(rp, rf)

    # parse hmm file for positive reads
    hmm_pos = parse_hmm(hm)

    # parse mm file and output read names / metrics for missing reads
    _ = parse_mm(mm, aligned_pos, roc_pos, hmm_pos, out)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')

if __name__ == "__main__":
    main()
