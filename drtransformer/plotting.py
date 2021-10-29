#!/usr/bin/env python
#
# drtransformer.plotting
# 
# All plotting utilities should be available through
# the commandline script DrPlotter --help
#

import logging
drlog = logging.getLogger(__name__)

import numpy as np
import matplotlib.pyplot as plt

from . import __version__
from .utils import make_pair_table, make_loop_index

def plot_simulation(trajectories, basename, formats, 
                    lin_time, log_time, tlen = None, motifs = None, title = ''):
    """DrTransformer standard plotting function.
    """

    fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1.8, 1]})
    fig.set_size_inches(7, 3)
    plt.subplots_adjust(wspace=0)

    [ax1, ax2] = axs
    ax1.set_xscale('linear')
    ax1.spines['right'].set_visible(False)
    ax1.set_xlim((0, lin_time))
    ax1.set_ylim([-0.02, 1.02])
    # Tick business
    xticks = [round(x, 2) for x in np.linspace(0, lin_time, 6)]
    ax1.set_xticks(xticks)

    ax2.set_xscale('log')
    ax2.spines['left'].set_visible(False)
    ax2.set_xlim((lin_time, log_time))
    ax2.set_ylim([-0.02, 1.02])
    # Tick business
    ax2.yaxis.tick_right()
    ax2.set_yticklabels([])

    ax1t = ax1.twiny()
    ax1t.spines['right'].set_visible(False)
    ax1t.set_xlim(ax1.get_xlim())
    # Tick business
    att, atl = zip(*tlen)
    tl = [int(x) for x in np.linspace(atl[0], atl[-1], 10)]
    tt = [t for e, (t, l) in enumerate(tlen) if l in tl]
    ax1t.set_xticks(tt)
    ax1t.set_xticklabels(tl)
    ax1t.set_xticks(att, minor = True)

    # Mark the end of transcription.
    ax1.axvline(x = lin_time, linewidth = 2, color = 'black', linestyle = '-') 
    ax2.axvline(x = lin_time, linewidth = 2, color = 'black', linestyle = '-') 

    # Set the grid lines.
    ax1.grid(axis = 'y', which = 'major', alpha = 0.7,
             color='gray', linestyle='--', linewidth=0.5)
    ax2.grid(axis = 'y', which = 'major', alpha = 0.7, 
             color='gray', linestyle='--', linewidth=0.5)
    ax1t.grid(which = 'major', alpha = 0.7,
             color='gray', linestyle='--', linewidth=0.5)

    # Plot the data.
    if motifs:
        lw = {'pp': 2, 'pm': .5, 'mm': 1}
        dh = {'pp': '-', 'pm': ':', 'mm': '--'}
        zo = {'pp': 1,   'pm': 2,    'mm': 3}
        for i, (name, _) in enumerate(motifs):
            for suffix in ['pp', 'pm', 'mm']:
                line = trajectories[f'{name}-{suffix}']
                t, o = list(zip(*line))
                if max(o) >= 0.1:
                    p, = ax1.plot(t, o, dh[suffix], lw=lw[suffix], 
                                 color = f'C{i}', zorder = zo[suffix])
                    l, = ax2.plot(t, o, dh[suffix], lw=lw[suffix], 
                                    color = f'C{i}', zorder = zo[suffix])
                    l.set_label(f'{name}-{suffix}')

    else:
        for ni in sorted(trajectories):
            line = trajectories[ni]
            t, o = list(zip(*line))
            # Determine which lines are part of the legend:
            # like this, it is only those that are populated
            # at the end of transcription and if they reach
            # an occupancy of 10% or higher
            if t[-1] > lin_time:
                p, = ax1.plot(t, o, '-', lw=1.5)
                l, = ax2.plot(t, o, '-', lw=1.5)
                if max(o) >= 0.1:
                    l.set_label(f"ID {ni}")
            elif max(o) >= 0.1:
                p, = ax1.plot(t, o, '-', lw=1)
                l, = ax2.plot(t, o, '-', lw=1)

    # Legends and labels.
    ax1.set_ylabel('occupancy')
    ax1.set_xlabel('time (seconds)')
    ax1.xaxis.set_label_coords(0.7, -0.15)
    ax1t.set_xlabel(f'{title} transcript length')
    ax1t.xaxis.set_label_coords(0.7, 1.15)
    ax2.legend()

    # Save a file.
    for ending in formats:
        pfile = basename + '.' + ending
        plt.savefig(pfile, bbox_inches = 'tight')
    return

def motif_finder(ss, motif):
    pt = make_pair_table(ss, base = 1)
    li = make_loop_index(ss)

    # for each base-pair, check if it is
    # contained or possible.
    found = 0
    exact = True
    for (l, r) in motif:
        assert l < r
        # check if incompatible
        if r > pt[0]:
            exact = False
            if l > pt[0]:
                return 'mm' # incompatible
            elif li[l-1] != 0:
                return 'mm' # incompatible
            else:
                pass # at most compatible
        else:
            if li[l-1] != li[r-1]:
                return 'mm' # incompatible
            elif pt[l] == r and pt[r] == l:
                found += 1
            else:
                exact = False
    if exact:
        return 'pp' # present
    elif found > 2:
        return 'pm' # compatible
    else:
        return 'mm' # incompatible

def plot_xmgrace(trajectories, filename):
    """Plot trajectories into a file in grace format.
    """
    head = """
@with line
@line on
@line loctype world
@line g0
@line linewidth 2
@line linestyle 1
@line color 7
@line arrow 0
@line arrow type 0
@line arrow length 1.000000
@line arrow layout 1.000000, 1.000000
@line def
"""
    with open(filename, 'w') as gfh:
        gfh.write(head)
        for nid in sorted(trajectories):
            course = trajectories[nid]
            t, o = list(zip(*course))
            for i in range(len(t)):
                gfh.write("{:f} {:f}\n".format(t[i], o[i]))
            gfh.write("&\n")
    return

def parse_drf_motifs(stream, motifs):
    xydata = {}
    for (name, m) in motifs:
        xydata[f'{name}-pp'] = [[0, 0]] # present
        xydata[f'{name}-pm'] = [[0, 0]] # compatible
        xydata[f'{name}-mm'] = [[0, 1]] # incompatible
    time_len = []
    llen, lint, ltime, lltime = 0, 0, '0', '0'
    for e, line in enumerate(stream):
        if e == 0:
            continue
        [idx, stime, occu, ss, en] = line.split()
        time = float(stime)
        occu = float(occu)
        for (name, m) in motifs:
            if stime != ltime: # Initialzie all possibilities to 0
                xydata[f'{name}-pp'].append([time, 0])
                xydata[f'{name}-pm'].append([time, 0])
                xydata[f'{name}-mm'].append([time, 0])
            suffix = motif_finder(ss, m) # none, true, false?
            xydata[f'{name}-{suffix}'][-1][1] += occu
        if len(ss) > llen:
            llen = len(ss)
            lint = float(lltime)
            time_len.append((float(lltime), len(ss)))
        if ltime != stime:
            ltime = stime
            lltime = ltime
    logt = time if time > lint + 60 else lint + 60
    return xydata, lint, logt, time_len

def parse_drf(stream):
    xydata = {}
    time_len = []
    llen, lint, ltime, lltime = 0, 0, '0', '0'
    for e, line in enumerate(stream):
        if e == 0:
            continue
        [idx, stime, occu, ss, en] = line.split()
        time = float(stime)
        occu = float(occu)
        xydata[idx] = (xydata.get(idx, [])) + [(time, occu)]
        if len(ss) > llen:
            llen = len(ss)
            lint = float(lltime)
            time_len.append((float(lltime), len(ss)))
        if ltime != stime:
            ltime = stime
            lltime = ltime
    logt = time if time > lint + 60 else lint + 60
    return xydata, lint, logt, time_len

def main():
    """ DrPlotter -- visulization of cotranscriptional folding simulations.
    """
    import sys, argparse
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'DrPlotter: visualization of the DrForna file format.')

    parser.add_argument('--version', action = 'version', 
            version = '%(prog)s ' + __version__)
    parser.add_argument("-v", "--verbose", action = 'count', default = 0,
            help = """Track process using verbose output.""")
    parser.add_argument("-n", "--name", default = 'DrPlotter', metavar = '<str>',
            help = """Name your output files, name the header of your plots, etc.
            this option overwrites the fasta-header.""")
    parser.add_argument("-m", "--motifs", default = None, metavar = '<str>',
            help = """Specify motifs using base-pairs m1=5-15,6-14:m2=5-81""")
    parser.add_argument("-f", "--formats", nargs = '+', default = ['pdf'],
            choices = ('pdf', 'svg', 'png', 'gr', 'eps'),
            help = """Plot the simulation using matplotlib (pdf, svg, png)
            and/or write an input file for xmgrace (gr). The legend uses 
            the identities of structures as specifed in the logfile. """)
    parser.add_argument("--molecule", default = '', metavar = '<str>',
            help = """Include molecule name in the plot. """)
    args = parser.parse_args()

    if args.motifs:
        motifs = []
        mstrings = args.motifs.split(':')
        for e, bpl in enumerate(mstrings):
            mot = []
            name, bpl = bpl.split('=')
            bps = bpl.split(',')
            for bp in bps:
                [l, r] = [int(x) for x in bp.split('-')]
                mot.append((l, r))
            motifs.append((name, mot))
        xydata, lint, logt, tlen = parse_drf_motifs(sys.stdin, motifs)
    else:
        motifs = None
        xydata, lint, logt, tlen = parse_drf(sys.stdin)

    if 'gr' in args.formats:
        plot_xmgrace(xydata, args.name + '.gr')

    mplf = [f for f in args.formats if f != 'gr']
    if mplf:
        with plt.style.context('seaborn-ticks'):
            plot_simulation(xydata, args.name, mplf, lint, logt, tlen, motifs, 
                    title = args.molecule)

 
if __name__ == '__main__':
    main()
