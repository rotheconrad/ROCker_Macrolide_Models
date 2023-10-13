#!/usr/bin/env python

''' Visualize read alignments and scores based on filter labels

passing, failing, TP, FP, TN, FN.

Reads blastn and blastx files mapping the same reads to the genes with
sequences in nucleotides (blastn) and in amino acids (blastx).

Reads a directory of fasta files from score_rocker_model.py where the 
filenames indicate the filtering result of its contained reads.

Creates plots to visualize the filtering labels and the alignment stats.

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

import argparse, glob
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns; sns.set(style="white", color_codes=True)
from scipy.stats import pearsonr as corr
from pathlib import Path


class TickRedrawer(matplotlib.artist.Artist):
    #https://stackoverflow.com/questions/19677963/
    #matplotlib-keep-grid-lines-behind-the-graph-but-the-y-and-x-axis-above
    """Artist to redraw ticks."""

    __name__ = "ticks"

    zorder = 10

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer: matplotlib.backend_bases.RendererBase) -> None:
        """Draw the ticks."""
        if not self.get_visible():
            self.stale = False
            return

        renderer.open_group(self.__name__, gid=self.get_gid())

        for axis in (self.axes.xaxis, self.axes.yaxis):
            loc_min, loc_max = axis.get_view_interval()

            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if tick.get_visible() and loc_min <= tick.get_loc() <= loc_max:
                    for artist in (tick.tick1line, tick.tick2line):
                        artist.draw(renderer)

        renderer.close_group(self.__name__)
        self.stale = False


def parse_blast(blast, blastx=False):

    data = {}

    # read input fasta and split
    with open(blast, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            name = X[0]
            pid = float(X[2])
            alen = int(X[3])
            sstart = int(X[8])
            ssend = int(X[9])
            midp = sstart + ssend / 2
            evalue = float(X[10])
            bitscore = float(X[11])
            qlen = int(X[12])
            # blastx amino acid correction
            if blastx == True:
                aqlen = qlen / 3
                afrac = (alen / aqlen) * 100
            else:
                afrac = (alen / qlen) * 100

            data[name] = [pid, evalue, bitscore, afrac, midp, qlen]

    return data


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def compile_data(fasta_dir, blastn, blastx, dmndx, outdir):

    # get list of fasta files
    fasta_files = [f for f in glob.glob(f'{fasta_dir}/*.fasta')]
    # create dictionaries to store data
    data = {
            'Filter Type': [], 'Filter Result': [], 'BLAST': [], 'Label': [],
            'Mid Point': [], 'Sequence Identity': [], 'E-value': [],
            'Bitscore': [], 'Aligned Fraction': [], 'Sim Overlap': [],
            'Read Length': [], 'Read Name': []
            }

    # parse fasta files and store read names, filter type, and filter result.
    for fasta in fasta_files:
        filter_type = fasta.split('/')[-1].split('-')[0]
        filter_result = fasta.split('/')[-1].split('-')[1].split('.')[0]

        with open(fasta, 'r') as file:
            for name, seq in read_fasta(file):
                name = name[1:] # remove '>' from fasta defline (read name)
                label = name.split(';')[-1]
                # percent overlap of simulated read with target gene (POSRTG)
                POSRTG = float(name.split(';')[4])
                # enter blastn stats
                bns = blastn.get(name, [-1, -1, -1, -1, -1, -1])
                pid, evalue, bitscore = bns[0], bns[1], bns[2]
                afrac, midp, qlen = bns[3], bns[4], bns[5]
                data['Filter Type'].append(filter_type)
                data['Filter Result'].append(filter_result)
                data['Sequence Identity'].append(pid)
                data['E-value'].append(evalue)
                data['Bitscore'].append(bitscore)
                data['Aligned Fraction'].append(afrac)
                data['BLAST'].append('BLASTn')
                data['Read Name'].append(name)
                data['Mid Point'].append(midp)
                data['Label'].append(label)
                data['Sim Overlap'].append(POSRTG)
                data['Read Length'].append(qlen)
                # enter blastx stats
                bxs = blastx.get(name, [-1, -1, -1, -1, -1, -1])
                pid, evalue, bitscore = bxs[0], bxs[1], bxs[2]
                afrac, midp, qlen = bxs[3], bxs[4], bxs[5]
                data['Filter Type'].append(filter_type)
                data['Filter Result'].append(filter_result)
                data['Sequence Identity'].append(pid)
                data['E-value'].append(evalue)
                data['Bitscore'].append(bitscore)
                data['Aligned Fraction'].append(afrac)
                data['BLAST'].append('BLASTx')
                data['Read Name'].append(name)
                data['Mid Point'].append(midp)
                data['Label'].append(label)
                data['Sim Overlap'].append(POSRTG)
                data['Read Length'].append(qlen)
                # enter diamond blastx stats
                dxs = dmndx.get(name, [-1, -1, -1, -1, -1, -1])
                pid, evalue, bitscore = dxs[0], dxs[1], dxs[2]
                afrac, midp, qlen = dxs[3], dxs[4], dxs[5]
                data['Filter Type'].append(filter_type)
                data['Filter Result'].append(filter_result)
                data['Sequence Identity'].append(pid)
                data['E-value'].append(evalue)
                data['Bitscore'].append(bitscore)
                data['Aligned Fraction'].append(afrac)
                data['BLAST'].append('DMNDx')
                data['Read Name'].append(name)
                data['Mid Point'].append(midp)
                data['Label'].append(label)
                data['Sim Overlap'].append(POSRTG)
                data['Read Length'].append(qlen)

    df = pd.DataFrame(data)
    df.to_csv(f'{outdir}/Compiled_Data.tsv', sep='\t', index=False)

    return df


def all_reads_bx_bn_plots(df, outdir):
    ''' compare blastn result with blastx result alignment position,
        pid, evalue, and bitscore '''

    subdir = f'{outdir}/01_All_Reads_bnbxdx'
    _ = Path(subdir).mkdir(parents=True, exist_ok=True)

    # Color dicts
    C = {'Positive': '#3182bd', 'Negative': '#de2d26', 'Non_Target': '#636363'}
    # markers
    M = '.'

    # organize the data
    dftemp = df[(df['Filter Type'] == 'all_reads')]
    dfbn = dftemp[dftemp['BLAST'] == 'BLASTn']
    dfbx = dftemp[dftemp['BLAST'] == 'BLASTx']
    dfdx = dftemp[dftemp['BLAST'] == 'DMNDx']

    ####################################################################
    ### MID POINT PLOT
    ####################################################################
    # all BLASTx vs BLASTn
    ptitle = 'BLASTx vs. BLASTn (Mid Point)'
    xlab = 'BLASTn (base pairs)'
    ylab = 'BLASTx (amino acids)'
    labs = dfbn['Label'].tolist()
    xs = dfbn['Mid Point'].tolist()
    ys = dfbx['Mid Point'].tolist()
    if len(xs) >= 1:
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/01a_BnBx_Midp.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # all DMNDx vs BLASTn
    ptitle = 'DMNDx vs. BLASTn (Mid Point)'
    xlab = 'BLASTn (base pairs)'
    ylab = 'DMNDx (amino acids)'
    labs = dfbn['Label'].tolist()
    xs = dfbn['Mid Point'].tolist()
    ys = dfdx['Mid Point'].tolist()
    if len(xs) >= 1:
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/01b_BnDx_Midp.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # all BLASTx vs DMNDx
    ptitle = 'BLASTx vs. DMNDx (Mid Point)'
    xlab = 'DMNDx (amino acids)'
    ylab = 'BLASTx (amino acids)'
    labs = dfbn['Label'].tolist()
    xs = dfdx['Mid Point'].tolist()
    ys = dfbx['Mid Point'].tolist()
    if len(xs) >= 1:
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/01c_DxBx_Midp.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    ####################################################################
    ### PID PLOT
    ####################################################################
    # all BLASTx vs BLASTn
    ptitle = 'BLASTx vs. BLASTn (Sequence Identity)'
    xlab = 'BLASTn (%)'
    ylab = 'BLASTx (%)'
    labs = dfbn['Label'].tolist()
    xs = dfbn['Sequence Identity'].tolist()
    ys = dfbx['Sequence Identity'].tolist()
    if len(xs) >= 1:
        xs, ys, labs = remove_no_alignments(xs, ys, labs)
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/02a_BnBx_pid.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # all DMNDx vs BLASTn
    ptitle = 'DMNDx vs. BLASTn (Sequence Identity)'
    xlab = 'BLASTn (%)'
    ylab = 'DMNDx (%)'
    labs = dfbn['Label'].tolist()
    xs = dfbn['Sequence Identity'].tolist()
    ys = dfdx['Sequence Identity'].tolist()
    if len(xs) >= 1:
        xs, ys, labs = remove_no_alignments(xs, ys, labs)
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/02b_BnDx_pid.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # all BLASTx vs DMNDx
    ptitle = 'BLASTx vs DMNDx (Sequence Identity)'
    xlab = 'DMNDx (%)'
    ylab = 'BLASTx (%)'
    labs = dfbn['Label'].tolist()
    xs = dfdx['Sequence Identity'].tolist()
    ys = dfbx['Sequence Identity'].tolist()
    if len(xs) >= 1:
        xs, ys, labs = remove_no_alignments(xs, ys, labs)
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/02c_DxBx_pid.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    ####################################################################
    ### EVALUE PLOT
    ####################################################################
    # all BLASTx vs BLASTn
    ptitle = 'BLASTx vs. BLASTn (E-value)'
    xlab = 'BLASTn'
    ylab = 'BLASTx'
    labs = dfbn['Label'].tolist()
    xs = dfbn['E-value'].tolist()
    ys = dfbx['E-value'].tolist()
    if len(xs) >= 1:
        xs, ys, labs = remove_no_alignments(xs, ys, labs)
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/03a_BnBx_ev.pdf'
        _ = scatter_plot(
                ptitle, xlab, ylab, xs, ys, C, cmap, outfile, log='both'
                )

    # all DMNDx vs BLASTn
    ptitle = 'DMNDx vs. BLASTn (E-value)'
    xlab = 'BLASTn'
    ylab = 'DMNDx'
    labs = dfbn['Label'].tolist()
    xs = dfbn['E-value'].tolist()
    ys = dfdx['E-value'].tolist()
    if len(xs) >= 1:
        xs, ys, labs = remove_no_alignments(xs, ys, labs)
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/03b_BnDx_ev.pdf'
        _ = scatter_plot(
                ptitle, xlab, ylab, xs, ys, C, cmap, outfile, log='both'
                )

    # all BLASTx vs DMNDx
    ptitle = 'BLASTx vs DMNDx (E-value)'
    xlab = 'DMNDx'
    ylab = 'BLASTx'
    labs = dfbn['Label'].tolist()
    xs = dfdx['E-value'].tolist()
    ys = dfbx['E-value'].tolist()
    if len(xs) >= 1:
        xs, ys, labs = remove_no_alignments(xs, ys, labs)
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/03c_DxBx_ev.pdf'
        _ = scatter_plot(
                ptitle, xlab, ylab, xs, ys, C, cmap, outfile, log='both'
                )

    ####################################################################
    ### BITSCORE PLOT
    ####################################################################
    # all BLASTx vs BLASTn
    ptitle = 'Bitscore'
    xlab = 'BLASTn'
    ylab = 'BLASTx'
    labs = dfbn['Label'].tolist()
    xs = dfbn['Bitscore'].tolist()
    ys = dfbx['Bitscore'].tolist()
    if len(xs) >= 1:
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/04a_BnBx_bts.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # all DMNDx vs BLASTn
    ptitle = 'DMNDx vs BLASTn (Bitscore)'
    xlab = 'BLASTn'
    ylab = 'DMNDx'
    labs = dfbn['Label'].tolist()
    xs = dfbn['Bitscore'].tolist()
    ys = dfdx['Bitscore'].tolist()
    if len(xs) >= 1:
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/04b_BnDx_bts.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # all BLASTx vs DMNDx
    ptitle = 'BLASTx vs DMNDx (Bitscore)'
    xlab = 'DMNDx'
    ylab = 'BLASTx'
    labs = dfbn['Label'].tolist()
    xs = dfdx['Bitscore'].tolist()
    ys = dfbx['Bitscore'].tolist()
    if len(xs) >= 1:
        cmap = [ C[i] for i in labs ]
        outfile = f'{subdir}/04c_DxBx_bts.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    return True


def alignment_failed_plots(df, outdir, ftype):
    ''' plot positive vs pid of positive reads that failed hmm '''

    switch = {
            'hmm': ['02_HMM_Fails', 'hmm_failed', 'HMM'],
            'DorBx': ['03_DorBx_Fails', 'blastx_failed', 'DorBx']
            }
    fdir = switch[ftype][0]
    subdir = f'{outdir}/{fdir}'
    _ = Path(subdir).mkdir(parents=True, exist_ok=True)

    # Color dicts
    C = {'BLASTn': '#3182bd', 'BLASTx': '#de2d26', 'DMNDx': '#7fbc41'}
    # hlines will be a dict of the lines to draw.
    hlines = {'e-30': 1e-30, 'e-20': 1e-20, 'e-10': 1e-10, 'e-03': 1e-03}

    # organize the data
    fselect = switch[ftype][1]
    dftemp = df[(df['Filter Type'] == fselect) & (df['Filter Result'] =='pos')]

    # Sequence Identity
    ftitle = switch[ftype][2]
    ptitle = f'Positive read failed {ftitle} (Sequence Identity)'
    xlab = 'Sequence position (base pairs or amino acids)'
    ylab = 'Sequence identity (%)'
    labs = dftemp['BLAST'].tolist()
    xs = dftemp['Mid Point'].tolist()
    ys = dftemp['Sequence Identity'].tolist()
    cmap = [ C[i] for i in labs ]
    if len(xs) >= 1:
        outfile = f'{subdir}/{fselect}_pos_pid.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # E-Value
    ptitle = f'Positive read failed {ftitle} (E-value)'
    xlab = 'Sequence position (base pairs or amino acids)'
    ylab = 'E-value'
    labs = dftemp['BLAST'].tolist()
    xs = dftemp['Mid Point'].tolist()
    ys = dftemp['E-value'].tolist()
    cmap = [ C[i] for i in labs ]
    if len(xs) >= 1:
        outfile = f'{subdir}/{fselect}_pos_evalue.pdf'
        _ = scatter_plot(
                ptitle, xlab, ylab, xs, ys, C, cmap, outfile,
                log='y', hlines=hlines
                )

    # Bitscore
    ptitle = f'Positive read failed {ftitle} (Bitscore)'
    xlab = 'Sequence position (base pairs or amino acids)'
    ylab = 'Bitscore'
    labs = dftemp['BLAST'].tolist()
    xs = dftemp['Mid Point'].tolist()
    ys = dftemp['Bitscore'].tolist()
    cmap = [ C[i] for i in labs ]
    if len(xs) >= 1:
        outfile = f'{subdir}/{fselect}_pos_bitscore.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # Percent overlap of simulated read with target gene
    ptitle = f'Positive read failed {ftitle} (POSRTG)'
    xlab = 'Sequence position (base pairs or amino acids)'
    ylab = 'Percent overlap of simulated read with target gene (%)'
    labs = dftemp['BLAST'].tolist()
    xs = dftemp['Mid Point'].tolist()
    ys = dftemp['Sim Overlap'].tolist()
    cmap = [ C[i] for i in labs ]
    if len(xs) >= 1:
        outfile = f'{subdir}/{fselect}_pos_POSRTG.pdf'
        _ = scatter_plot(ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    return True


def filter_result_vs_filter_type(df, outdir):
    ''' plot positive vs pid of positive reads that failed hmm '''

    subdir = f'{outdir}/04_Filter_Type_Results'
    _ = Path(subdir).mkdir(parents=True, exist_ok=True)

    # Color dicts
    C = {'TP': '#3182bd', 'FP': '#de2d26', 'TN': '#7fbc41', 'FN': '#c51b7d'}
    #C = {'Positive': '#3182bd', 'Negative': '#de2d26', 'Non_Target': '#636363'}
    # hlines will be a dict of the lines to draw.
    hlines = {'e-30': 1e-30, 'e-20': 1e-20, 'e-10': 1e-10, 'e-03': 1e-03}
    # filter name values to plot
    filters = ['roc', 'cA', 'cB', 'e30', 'e20', 'e10', 'e03']
    markers = {'Positive': '.', 'Negative': 'x', 'Non_Target': '_'}
    results = df['Filter Result'].unique()
    blasts = df['BLAST'].unique()

    for fil in filters:
        for blast in blasts:
            # organize the data
            dftemp = df[(df['Filter Type'] == fil) & (df['BLAST'] == blast)]
            labs = dftemp['Filter Result'].tolist()
            xs = dftemp['Mid Point'].tolist()
            cmap = [ C[i] for i in labs ]
            xlabdict = {
                'DMNDx': 'Sequence position (amino acids)',
                'BLASTx': 'Sequence position (amino acids)',
                'BLASTn': 'Sequence position (base pairs)'
                }
            xlabel = xlabdict[blast]

            # E-value
            ptitle = f'{blast}-{fil} (E-value)'
            xlab = xlabel
            ylab = 'E-value'
            ys = dftemp['E-value'].tolist()
            outfile = f'{subdir}/{blast}-{fil}_ev.pdf'
            _ = scatter_plot(
                ptitle, xlab, ylab, xs, ys, C, cmap, outfile,
                hlines=hlines, log='y'
                )
            # Sequence Identity
            ptitle = f'{blast}-{fil} (Sequence Identity)'
            xlab = xlabel
            ylab = 'Sequence identity (%)'
            ys = dftemp['Sequence Identity'].tolist()
            if len(xs) >= 1:
                outfile = f'{subdir}/{blast}-{fil}_pid.pdf'
                _ = scatter_plot(
                    ptitle, xlab, ylab, xs, ys, C, cmap, outfile)
            # Bitscore
            ptitle = f'{blast}-{fil} (Bitscore)'
            xlab = xlabel
            ylab = 'Bitscore'
            ys = dftemp['Bitscore'].tolist()
            if len(xs) >= 1:
                outfile = f'{subdir}/{blast}-{fil}_bts.pdf'
                _ = scatter_plot(
                    ptitle, xlab, ylab, xs, ys, C, cmap, outfile)
            # Aligned Fraction
            ptitle = f'{blast}-{fil} (Aligned Fraction)'
            xlab = xlabel
            ylab = 'Aligned fraction (%)'
            ys = dftemp['Aligned Fraction'].tolist()
            if len(xs) >= 1:
                outfile = f'{subdir}/{blast}-{fil}_alfr.pdf'
                _ = scatter_plot(
                    ptitle, xlab, ylab, xs, ys, C, cmap, outfile)
                
    return True


def remove_no_alignments(xs, ys, labs):

    # remove any positions in the list where xs or ys has a -1
    # -1 is for no alignment
    xt, yt, lt = [], [], []

    for x, y, l in zip(xs, ys, labs):
        if x != -1 and y != -1:
            xt.append(x)
            yt.append(y)
            lt.append(l)

    return xt, yt, lt


def scatter_plot(
    ptitle, xlab, ylab, xs, ys, C, cmap, outfile, log=None, hlines=None
    ):

    # number of data points
    n = len(xs)
    # Set Colors and markers
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'
    color = '#252525'
    marker = 'o'
    p = 50 # point size
    a = 0.25 # alpha

    # build plot
    g = sns.JointGrid(x=xs, y=ys, palette=cmap)
    # x margin kde plot
    sns.kdeplot(
            x=xs,
            ax=g.ax_marg_x,
            legend=False,
            color=color
            )
    # y margin kde plot
    sns.kdeplot(
            y=ys,
            ax=g.ax_marg_y,
            legend=False,
            color=color
            )
    # main scatter plot
    rast = True if n >= 1000000 else False
    g.ax_joint.scatter(
            xs,
            ys,
            marker=marker,
            s=p,
            alpha=a,
            rasterized=rast,
            c=cmap,
            edgecolors='none'
            )

    # plot title, labels, text
    g.ax_marg_x.set_title(ptitle, fontsize=18, y=1.02)

    g.ax_joint.set_xlabel(
        xlab,
        fontsize=12, y=-0.02
        )
    g.ax_joint.set_ylabel(
        ylab,
        fontsize=12, x=-0.02
        )
    # Gather correlation
    if n >= 2:
        pcorr = corr(xs, ys)
        stats_line = (
            f"Pearson r: {round(pcorr[0], 2)}; "
            f"p value: {round(pcorr[1], 2)}"
            )
        g.ax_joint.text(
            0.05, 0.99, stats_line,
            fontsize=10, color=second_color,
            verticalalignment='top', horizontalalignment='left',
            transform=g.ax_joint.transAxes
            )

    # log axis?
    if log == 'both':
        g.ax_joint.set_yscale('log')
        g.ax_joint.invert_yaxis()
        g.ax_joint.set_xscale('log')
    elif log == 'y':
        g.ax_joint.set_yscale('log')
        g.ax_joint.invert_yaxis()

    # hlines?
    if hlines:
        for k, v in hlines.items():
            g.ax_joint.axhline(
                y=v, xmin=0, xmax=1,
                color=vline_color, linewidth=2, linestyle='--'
                )
            g.ax_joint.text(
                g.ax_joint.get_xlim()[1]-5, v, k, fontsize=10, color='red',
                verticalalignment='top', horizontalalignment='right',
                #transform=g.ax_joint.transAxes
                )

    # set the axis parameters / style
    g.ax_joint.tick_params(
        axis='both', which='major', direction='inout', color='k',
        width=2, length=6, bottom=True, left=True, zorder=3
        )

    # set grid style
    g.ax_joint.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    g.ax_joint.xaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    g.ax_joint.set_axisbelow(True)
    g.ax_joint.add_artist(TickRedrawer())

    # build the legend
    legend_labels = []
    legend_elements = []

    nl = '(-1) No Alignment'
    n = Line2D(
            [0], [0], color='w', label=nl, marker=marker,
            markersize=12, markerfacecolor='w'
            )
    legend_labels.append(nl)
    legend_elements.append(n)

    for k,v in C.items():
        n = Line2D(
            [0], [0], color='w', label=k, marker=marker,
            markersize=12, markerfacecolor=v
            )
        legend_labels.append(k)
        legend_elements.append(n)

    g.ax_joint.legend(
        handles=legend_elements,
        fontsize=8,
        fancybox=True,
        framealpha=0.0,
        frameon=False,
        loc='lower right',
        #bbox_to_anchor=(0, 0.98, 1, 0.2),
        ncol=len(legend_labels)
        )

    # adjust layout, save, and close
    g.fig.set_figwidth(7)
    g.fig.set_figheight(5)
    g.savefig(outfile)
    plt.close()

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-bn', '--blastn_file',
        help='Please specify the blastn file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-bx', '--blastx_file',
        help='Please specify the blastx file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-dx', '--dmndx_file',
        help='Please specify the Diamond blastx file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-fd', '--fasta_dir',
        help='Please specify the fasta directory!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-od', '--output_directory',
        help='Please specify a name for the output directory!',
        metavar='',
        type=str,
        required=True
        )

    args=vars(parser.parse_args())

    # define input params
    blastn_file = args['blastn_file']
    blastx_file = args['blastx_file']
    dmndx_file = args['dmndx_file']
    fasta_dir = args['fasta_dir']
    outdir = args['output_directory']

    # Do what you came here to do:
    print('\n\nRunning Script...')

    # parse the blast files
    blastn = parse_blast(blastn_file)
    blastx = parse_blast(blastx_file, blastx=True)
    dmndx = parse_blast(dmndx_file, blastx=True)

    # create output directory
    _ = Path(outdir).mkdir(parents=True, exist_ok=True)

    # create large data frame
    df = compile_data(fasta_dir, blastn, blastx, dmndx, outdir)

    # Create plots from the df
    ################

    # All simulated reads blastx vs blastn vs diamond x
    print('\n\nBuilding plot set 1:')
    print('\tAll simulated reads, blastx vs. blastn vs. diamondx')
    _ = all_reads_bx_bn_plots(df, outdir)

    # plot hmm failed positive read alignments
    print('\n\nBuilding plot set 2:')
    print('\tAlignment metrics for reads failing HMM alignment')
    _ = alignment_failed_plots(df, outdir, 'hmm')
    print('\tAlignment metrics for reads failing Diamond/BLASTx alignment')
    _ = alignment_failed_plots(df, outdir, 'DorBx')
    # plot filter results vs filter types by alignment position
    print('\n\nBuilding plot set 3:')
    print('\tAlignment metrics by filter type and filter result')
    _ = filter_result_vs_filter_type(df, outdir)
    
    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')

if __name__ == "__main__":
    main()

