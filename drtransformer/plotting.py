#!/usr/bin/env python
#
# drtransformer.plotting
# 
# All plotting utilities should be available through
# the commandline script DrPlotter --help
#

import logging
drlog = logging.getLogger(__name__)

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from . import __version__
from .utils import make_pair_table, make_loop_index

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

def plot_simulation(trajectories, basename, formats, 
                    lin_time, log_time, motifs = None, title = ''):
    """DrTransformer simulation plotting.
    """
    assert log_time >= lin_time * 10

    # Do the plotting
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['right'].set_visible(False)
    ax.set_ylim([-0.05, 1.05])
    ax.set_xscale('linear')

    # Make the second part of the plot logarithmic
    offset = 0.00001
    ax.set_xlim((0, lin_time + offset))
    divider = make_axes_locatable(ax)
    axLog = divider.append_axes("right", size=2.5, pad=0, sharey=ax)
    axLog.set_xscale('log')
    axLog.set_xlim((lin_time + offset, log_time))
    axLog.set_ylim([-0.05, 1.05])
    axLog.yaxis.set_visible(False)
    axLog.spines['left'].set_visible(False)

    if motifs:
        lw = {'pp': 2, 'pm': .5, 'mm': 1}
        dh = {'pp': '-', 'pm': ':', 'mm': '--'}
        zo = {'pp': 1,   'pm': 2,    'mm': 3}
        for i, (name, _) in enumerate(motifs):
            for suffix in ['pp', 'pm', 'mm']:
                line = trajectories[f'{name}-{suffix}']
                t, o = list(zip(*line))
                if max(o) >= 0.1:
                    p, = ax.plot(t, o, dh[suffix], lw=lw[suffix], 
                                 color = f'C{i}', zorder = zo[suffix])
                    L, = axLog.plot(t, o, dh[suffix], lw=lw[suffix], 
                                    color = f'C{i}', zorder = zo[suffix])
                    L.set_label(f'{name}-{suffix}')

    else:
        for ni in sorted(trajectories):
            line = trajectories[ni]
            t, o = list(zip(*line))
            # Determine which lines are part of the legend:
            # like this, it is only those that are populated
            # at the end of transcription and if they reach
            # an occupancy of 10% or higher
            if t[-1] > lin_time:
                p, = ax.plot(t, o, '-', lw=1.5)
                L, = axLog.plot(t, o, '-', lw=1.5)
                if max(o) >= 0.1:
                    L.set_label(f"ID {ni}")
            elif max(o) >= 0.1:
                p, = ax.plot(t, o, '-', lw=0.5)
                L, = axLog.plot(t, o, '-', lw=0.5)

    fig.set_size_inches(7, 3)
    fig.text(0.5, 0.95, title, ha='center', va='center')

    # Add ticks for 1 minute, 1 hour, 1 day, 1 year
    axLog.axvline(x = lin_time, linewidth = 3, color = 'black', linestyle = '-') 
    axLog.axvline(x = 60, linewidth = 1, color = 'black', linestyle = '--') # 1 minute
    axLog.axvline(x = 3600, linewidth = 1, color = 'black', linestyle = '--') # 1 hour
    axLog.axvline(x = 86400, linewidth = 1, color = 'black', linestyle = '--') # 1 day
    axLog.axvline(x = 31536000, linewidth = 1, color = 'black', linestyle = '--') # 1 year
    plt.legend()

    ax.set_ylabel('occupancy', fontsize = 11)
    ax.set_xlabel('time [s]', ha = 'center', va = 'center', fontsize = 11)
    ax.xaxis.set_label_coords(.9, -0.15)

    for ending in formats:
        pfile = basename + '.' + ending
        plt.savefig(pfile, bbox_inches = 'tight')
    return

def parse_drf_motifs(stream, motifs):
    xydata = {}
    for (name, m) in motifs:
        xydata[f'{name}-pp'] = [[0, 0]] # present
        xydata[f'{name}-pm'] = [[0, 0]] # compatible
        xydata[f'{name}-mm'] = [[0, 1]] # incompatible

    llen, lint, ltime = 0, 0, '0'
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
            lint = time
            llen = len(ss)
        ltime = stime
    logt = time if time > 10 * lint else 10 * lt
    return xydata, lint, logt

def parse_drf(stream):
    xydata = {}
    llen, lint = 0, 0
    for e, line in enumerate(stream):
        if e == 0:
            continue
        [idx, time, occu, ss, en] = line.split()
        time = float(time)
        occu = float(occu)
        xydata[idx] = (xydata.get(idx, [])) + [(time, occu)]
        if len(ss) > llen:
            lint = time
            llen = len(ss)
    logt = time if time > 10 * lint else 10 * lt
    return xydata, lint, logt

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
        xydata, lint, logt = parse_drf_motifs(sys.stdin, motifs)
    else:
        motifs = None
        xydata, lint, logt = parse_drf(sys.stdin)

    if 'gr' in args.formats:
        plot_xmgrace(xydata, args.name + '.gr')

    mplf = [f for f in args.formats if f != 'gr']
    if mplf:
        with plt.style.context('ggplot'):
            plot_simulation(xydata, args.name, mplf, lint, logt, motifs, title = args.name)

 
if __name__ == '__main__':
    main()
