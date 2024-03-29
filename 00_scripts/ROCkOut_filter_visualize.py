#!/usr/bin/env python

''' Visualize ROCkOut align and filter results.

Classify results as TP, FP, TN, FN and build a plot.

Reads passing and failing rockout alignments from mock metagenomes creating
from ROCkOut read simulation and labeling.

Creates plots to visualize the filtering labels and the alignment stats.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Nov 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns; sns.set(style="white", color_codes=True)


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


def parse_alignment(alignment):

    data = {}

    with open(alignment, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            read = X[0]
            pid = float(X[2]) # percent sequence identity
            alen = int(X[3]) # alignment length in amino acids
            mm = int(X[4]) # number of mismatches
            gp = int(X[5]) # number of gap openings
            ev = float(X[10]) # e-value
            bs = float(X[11]) # bitscore
            qlen = int(X[12]) / 3 # full length of read divided by 3 amino acids

            pmatch = alen / qlen # percent match length of read length
            data[read] = [pid, alen, mm, gp, ev, bs, pmatch]

    return data


def parse_blast(passing, failing, alignment):

    data = {
        'Classifier': [], 'Midpoint': [], 'pID': [], 'Bitscore': [],
        'Mismatches': [], 'Gaps': [], 'evalue': [], 'Alignment length': [],
        'Percent alignment': [], 'alignment pid': [], 'alignment bitscore': []
        }

    # read passing alignments file into data dict
    with open(passing, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            read = X[0]
            label = X[0].split(';')[-1]
            pid = float(X[2])
            bitscore = float(X[8])
            sstart = min([int(X[6]), int(X[7])])
            ssend = max([int(X[6]), int(X[7])])

            # get additional data from alignment file
            Y = alignment[read]
            ypid = Y[0]
            alen = Y[1]
            mm = Y[2]
            gp = Y[3]
            ev = Y[4]
            ybs = Y[5]
            pmatch = Y[6]

            #print(read, pid, ypid, bitscore, ybs)

            midp = (sstart + ssend) / 2
            classifier = 'TP' if label == 'Positive' else 'FP'

            data['Classifier'].append(classifier)
            data['Midpoint'].append(midp)
            data['pID'].append(pid)
            data['Bitscore'].append(bitscore)
            data['Mismatches'].append(mm)
            data['Gaps'].append(gp)
            data['evalue'].append(ev)
            data['Alignment length'].append(alen)
            data['Percent alignment'].append(pmatch)
            data['alignment pid'].append(ypid)
            data['alignment bitscore'].append(ybs)

    # read failing alignments file into data dict
    with open(failing, 'r') as file:
        for line in file:
            read = X[0]
            X = line.rstrip().split('\t')
            label = X[0].split(';')[-1]
            pid = float(X[2])
            bitscore = float(X[8])
            sstart = min([int(X[6]), int(X[7])])
            ssend = max([int(X[6]), int(X[7])])

            # get additional data from alignment file
            Y = alignment[read]
            ypid = Y[0]
            alen = Y[1]
            mm = Y[2]
            gp = Y[3]
            ev = Y[4]
            ybs = Y[5]
            pmatch = Y[6]

            #print(read, pid, ypid, bitscore, ybs)

            midp = (sstart + ssend) / 2
            classifier = 'FN' if label == 'Positive' else 'TN'

            data['Classifier'].append(classifier)
            data['Midpoint'].append(midp)
            data['pID'].append(pid)
            data['Bitscore'].append(bitscore)
            data['Mismatches'].append(mm)
            data['Gaps'].append(gp)
            data['evalue'].append(ev)
            data['Alignment length'].append(alen)
            data['Percent alignment'].append(pmatch)
            data['alignment pid'].append(ypid)
            data['alignment bitscore'].append(ybs)

    df = pd.DataFrame(data)

    return df


def build_plot(df, outpre):

    # Color dicts
    C = {'TP': '#3182bd', 'FP': '#de2d26', 'TN': '#7fbc41', 'FN': '#c51b7d'}
    labs = df['Classifier'].tolist()
    xs = df['Midpoint'].tolist()
    xlab = 'Gene position (amino acids)'
    cmap = [ C[i] for i in labs ]

    # Sequence Identity - a
    dfa = df[['Midpoint', 'pID', 'Classifier']]
    dfa.columns = ['xs', 'ys', 'hue']
    ptitle = f'Sequence Identity'
    ylab = 'Sequence identity (%)'
    ys = df['pID'].tolist()
    outfile = f'{outpre}_pid.pdf'
    _ = scatter_plot(dfa, ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # Bitscore - b
    dfb = df[['Midpoint', 'Bitscore', 'Classifier']]
    dfb.columns = ['xs', 'ys', 'hue']
    ys = df['Bitscore'].tolist()
    ptitle = f'Bitscore'
    ylab = 'Bitscore'
    outfile = f'{outpre}_bts.pdf'
    _ = scatter_plot(dfb, ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # Mismatches - c
    dfc = df[['Midpoint', 'Mismatches', 'Classifier']]
    dfc.columns = ['xs', 'ys', 'hue']
    ys = df['Mismatches'].tolist()
    ptitle = f'Mismatches'
    ylab = 'Mismatches'
    outfile = f'{outpre}_Mismatches.pdf'
    _ = scatter_plot(dfc, ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # Gaps - d
    dfd = df[['Midpoint', 'Gaps', 'Classifier']]
    dfd.columns = ['xs', 'ys', 'hue']
    ys = df['Gaps'].tolist()
    ptitle = f'Gaps'
    ylab = 'Gaps'
    outfile = f'{outpre}_Gaps.pdf'
    _ = scatter_plot(dfd, ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # evalue - e
    dfe = df[['Midpoint', 'evalue', 'Classifier']]
    dfe.columns = ['xs', 'ys', 'hue']
    ys = df['evalue'].tolist()
    ptitle = f'evalue'
    ylab = 'evalue'
    outfile = f'{outpre}_evalue.pdf'
    _ = scatter_plot(dfe, ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # Alignment length - f
    dff = df[['Midpoint', 'Alignment length', 'Classifier']]
    dff.columns = ['xs', 'ys', 'hue']
    ys = df['Alignment length'].tolist()
    ptitle = f'Alignment length'
    ylab = 'Alignment length'
    outfile = f'{outpre}_alen.pdf'
    _ = scatter_plot(dff, ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # Percent alignment - g
    dfg = df[['Midpoint', 'Percent alignment', 'Classifier']]
    dfg.columns = ['xs', 'ys', 'hue']
    ys = df['Percent alignment'].tolist()
    ptitle = f'Alignment length / read length'
    ylab = 'Alignment length (%)'
    outfile = f'{outpre}_prcnt-alen.pdf'
    _ = scatter_plot(dfg, ptitle, xlab, ylab, xs, ys, C, cmap, outfile)


    # alignment vs filter pid - h
    xs = df['alignment pid'].tolist()
    xlab = 'Sequence identity (%) from alignment'
    dfh = df[['alignment pid', 'pID', 'Classifier']]
    dfh.columns = ['xs', 'ys', 'hue']
    ptitle = f'Alignment vs. Filter: Sequence Identity'
    ylab = 'Sequence identity (%) from filter'
    ys = df['pID'].tolist()
    outfile = f'{outpre}_alignment-vs-pid.pdf'
    _ = scatter_plot(dfh, ptitle, xlab, ylab, xs, ys, C, cmap, outfile)

    # alignment vs filter bitscore - i
    xs = df['alignment bitscore'].tolist()
    xlab = 'Bitscore from alignment'
    dfi = df[['alignment bitscore', 'Bitscore', 'Classifier']]
    dfi.columns = ['xs', 'ys', 'hue']
    ptitle = f'Alignment vs. Filter: Bitscore'
    ylab = 'Bitscore from filter'
    ys = df['pID'].tolist()
    outfile = f'{outpre}_alignment-vs-bitscore.pdf'
    _ = scatter_plot(dfi, ptitle, xlab, ylab, xs, ys, C, cmap, outfile)


    return True


def scatter_plot(df, ptitle, xlab, ylab, xs, ys, C, cmap, outfile):

    # number of data points
    n = len(xs)
    # Set Colors and markers
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'
    color = '#252525'
    marker = '.'
    p = 10 # point size
    a = 0.25 # alpha

    # build plot
    g = sns.JointGrid(x=xs, y=ys, palette=cmap)
    # x margin kde plot
    sns.kdeplot(
            data=df,
            x='xs',
            hue='hue',
            palette=C,
            ax=g.ax_marg_x,
            legend=False,
            )
    # y margin kde plot
    sns.kdeplot(
            data=df,
            y='ys',
            hue='hue',
            palette=C,
            ax=g.ax_marg_y,
            legend=False,
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
    xmin, xmax = min(xs), max(xs)
    xoffset = round((xmax - xmin) / 20)
    g.ax_joint.set_xlim([xmin - xoffset, xmax + xoffset])
    ymin, ymax = min(ys), max(ys)
    yoffset = round((ymax - ymin) / 20)
    g.ax_joint.set_ylim([ymin - yoffset, ymax + yoffset+yoffset])
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
        loc='upper center',
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
        '-p', '--passing_file',
        help='Please specify the passing alignments file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-f', '--failing_file',
        help='Please specify the failing alignments file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-a', '--alignment_file',
        help='Please specify the original alignments file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the output file (use .pdf)!',
        metavar='',
        type=str,
        required=True
        )

    args=vars(parser.parse_args())

    # define input params
    passing = args['passing_file']
    failing = args['failing_file']
    alignment = args['alignment_file']
    outpre = args['output_file']

    # Do what you came here to do:
    print('\n\nRunning Script...')

    # parse the alignment files
    data = parse_alignment(alignment)
    df = parse_blast(passing, failing, data)

    _ = build_plot(df, outpre)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')

if __name__ == "__main__":
    main()

