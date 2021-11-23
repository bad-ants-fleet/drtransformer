#!/usr/bin/env python
#
# DrTransformer -- cotranscriptional folding.
#
import logging
import os, sys, argparse
import math
import numpy as np
from datetime import datetime # Performance report 

import RNA
from . import __version__
from .landscape import TrafoLandscape
from .rnafolding import top_down_coarse_graining, parse_model_details
from .utils import parse_vienna_stdin, get_tkn_simulation_files

def restricted_float(x):
    y = float(x)
    if y < 0.0 or y > 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0]")
    return y

def parse_drtrafo_args(parser):
    """ A collection of arguments that are used by DrTransformer """

    environ = parser.add_argument_group('DrTransformer dependencies')
    output = parser.add_argument_group('DrTransformer output')
    trans = parser.add_argument_group('Transcription parameters')
    algo  = parser.add_argument_group('DrTransformer algorithm')

    ###################
    # Default options #
    ###################
    parser.add_argument('--version', action = 'version', 
            version = '%(prog)s ' + __version__)

    parser.add_argument("-v", "--verbose", action = 'count', default = 0,
            help = """Track process by writing verbose output to STDOUT during
            calculations. Use --logfile if you want to see *just* verbose
            information via STDOUT.""")

    ########################
    # DrTransformer output #
    ########################
    output.add_argument("--name", default = '', metavar = '<str>',
            help = """Name your output files, name the header of your plots, etc.
            this option overwrites the fasta-header.""")

    output.add_argument("--stdout", default = None, action = 'store',
            choices = ('log', 'drf', 'OFF'),
            help = """Choose STDOUT formats to follow the cotranscriptional
            folding progress in real time: *log*: a human readable output
            format.  *drf*: DrForna visualization input format. *OFF*: actively
            suppress output. The default (None) switches between *OFF*, if --logfile is
            specified, or *log* otherwise.""")

    output.add_argument("--logfile", action = "store_true",
            help = """Write verbose information to a file:
            {--outdir}/{--name}.log""")

    output.add_argument("--outdir", default = '', action = 'store', metavar = '<str>',
            help = """Place regular output files, into this directory. Creates
            the directory if it does not exist. """)

    output.add_argument("--tmpdir", default = '', action = 'store', metavar = '<str>',
            help = """Specify path for storing landscape files for debugging. These
            files will not be removed when the program terminates. """)

    output.add_argument("--no-timecourse", action = "store_true",
            help = """Do not produce the time-course file (outdir/name.drf).""")

    output.add_argument("--plot-minh", type = float, default = None, metavar = '<flt>',
            help = """Coarsify the output based on a minimum barrier height. In contrast
            to t-fast, this does *not* affect the accuracy of the model.""")

    output.add_argument("--t-lin", type = int, default = 30, metavar = '<int>',
            help = """Evenly space output *--t-lin* times during transcription on a linear time scale.""")

    output.add_argument("--t-log", type = int, default = 300, metavar = '<int>',
            help = """Evenly space output *--t-log* times after transcription on a logarithmic time scale.""")

    output.add_argument("--performance-report", action = "store_true", 
            # prints datetimes and statprof data for performance analysis
            help = argparse.SUPPRESS)

    ############################
    # Transcription parameters #
    ############################
    trans.add_argument("--t-ext", type = float, default = 0.02, metavar = '<flt>',
            help = """Inverse of transcription rate, i.e. time per nucleotide extension
            [seconds per nucleotide].""")

    trans.add_argument("--t-end", type = float, default = 3600, metavar = '<flt>',
            help = "Post-transcriptional simulation time [seconds].")

    trans.add_argument("--pause-sites", nargs='+', metavar='<int>=<flt>',
            help="""Transcriptional pausing sites.  E.g. \"--pause-sites 82=2e3
            84=33\" alters the simulation time at nucleotides 82 and 84 to 2000
            and 33 seconds, respectively. """)

    trans.add_argument("--start", type = int, default = 1, metavar = '<int>',
            help = "Start transcription at this nucleotide.")

    trans.add_argument("--stop", type = int, default = None, metavar = '<int>',
            help = "Stop transcription at this nucleotide")

    ###########################
    # DrTransformer algorithm #
    ###########################
    algo.add_argument("--o-prune", type = restricted_float, default = 0.05, metavar = '<flt>',
            help = """Occupancy threshold to prune structures from the 
            network. The structures with lowest occupancy are removed until
            at most o-prune occupancy has been removed from the total population. """)

    algo.add_argument("--t-fast", type = float, default = 0.001, metavar = '<flt>',
            help = """Folding times faster than --t-fast are considered
            instantaneous.  Structural transitions that are faster than
            --t-fast are considered part of the same macrostate. Directly
            translates to an energy barrier separating conformations using:
            dG = -RT*ln((1/t-fast)/k0). None: t-fast = 1/k_0 """)

    algo.add_argument("--minh", type = float, default = None, metavar = '<flt>',
            # An alternative to specify --t-fast in terms of a barrier height.
            help = argparse.SUPPRESS)

    algo.add_argument("--force", action = "store_true", 
            # Enforce a setting against all warnings.
            help = argparse.SUPPRESS)

    algo.add_argument("--fpwm", type = int, default = 4, metavar = '<int>',
            help = """Findpath search width multiplier. This value times the
            base-pair distance results in the findpath search width. Higher
            values increase the chances to find energetically lower saddle
            energies for the transition. """)

    algo.add_argument("--mfree", type = int, default = 6, metavar = '<int>',
            help = """Minimum number of freed bases during helix fraying.
            Fraying helices can vary greatly in length, starting with at
            least two base-pairs. This parameter defines the minimum amount of
            bases freed by helix fraying. For example, 6 can correspond to a
            stack of two base-pairs and a loop region of 2 nucleotides. If less
            bases are freed and there exists a nested stacked helix, this helix
            is considered to fray as well.""")

    algo.add_argument("--delth", type = int, default = 10, metavar = '<int>',
            help = """Delete structures if inactive for more than *delth* rounds.""")

    algo.add_argument("--k0", type = float, default = 2e5, metavar = '<flt>',
            help = """Arrhenius rate constant (pre-exponential factor). Adjust
            this constant of the Arrhenius equation to relate free energy
            changes to experimentally determined folding time [atu/sec].""")
    return

def write_output(data, stdout = False, fh = None):
    # Helper function to print data to filehandle, STDOUT, or both.
    if stdout:
        sys.stdout.write(data)
    if fh:
        fh.write(data)
    return

def set_handle_verbosity(h, v):
    if v == 0:
        h.setLevel(logging.WARNING)
    elif v == 1:
        h.setLevel(logging.INFO)
    elif v == 2:
        h.setLevel(logging.DEBUG)
    elif v >= 3:
        h.setLevel(logging.NOTSET)

def main():
    """ DrTransformer - cotranscriptional folding. 
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'DrTransformer: RNA folding kinetics during transcription.')
    parse_drtrafo_args(parser)
    parse_model_details(parser)
    args = parser.parse_args()

    # ~~~~~~~~~~~~~
    # Logging Setup 
    # ~~~~~~~~~~~~~
    title = f'DrTransformer-v{__version__}: RNA folding kinetics during transcription.'
    logger = logging.getLogger('drtransformer')
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    formatter = logging.Formatter('# %(levelname)s - %(message)s')
    set_handle_verbosity(ch, args.verbose)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info(title)

    (name, fullseq) = parse_vienna_stdin(sys.stdin, chars='ACGUNacgun')

    # Adjust arguments, prepare simulation
    if args.plot_minh:
        args.plot_minh = int(round(args.plot_minh*100))

    if args.name == '':
        args.name = name
    else:
        name = args.name

    if os.path.split(args.name)[0]:
        raise SystemExit('ERROR: Argument "--name" must not contain file path.')

    if args.outdir:
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        filepath = args.outdir + '/' + args.name
    else:
        filepath = args.name

    dfh = open(filepath + '.drf', 'w') if not args.no_timecourse else None
    lfh = open(filepath + '.log', 'w') if args.logfile else None

    # Adjust file handle
    if args.stdout is None and lfh is None:
        args.stdout = 'log'

    if args.tmpdir:
        _tmpdir = args.tmpdir
        if not os.path.exists(_tmpdir):
            os.makedirs(_tmpdir)

    # Adjust simulation parameters
    _RT = 0.61632077549999997
    if args.temp != 37.0:
        kelvin = 273.15 + args.temp
        _RT = (_RT / 310.15) * kelvin

    if args.stop is None:
        args.stop = len(fullseq)
    elif args.stop > len(fullseq):
        raise SystemExit(f'Invalid argument {(args.stop > len(fullseq))=}')

    if not 1 <= args.start <= args.stop:
        raise SystemExit(f'Invalid argument {(1 <= args.start <= args.stop)=}')

    if args.t_end < args.t_ext:
        raise SystemExit('ERROR: Conflicting Settings: ' + \
                'Arguments must be such that "--t-end" >= "--t-ext"')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Adjust the simulation window for treekin:
    #
    #   k0
    #   2e5  atu/s    50 nuc/s               1e-inf /s
    #   4000 /ext     1  nuc/ext             1e-inf /ext
    #   |------$------|-------------------$----------->
    #          k-fast                     k-slow
    #          --minh                     --maxh
    #          |------|---------------|--->
    #          t0     t_ext 0.02 s    t_end = 86400 s
    #   <----->    simulation             <---------->
    #   instant                             rate = 0
    #
    # (1) The rate of a spontaneous folding event
    # (>= k-fast) has to be much faster than the
    # rate of chain elongation (kx).
    #
    # (2) The rate of an effectively 0 folding
    # event (< k-slow) has to be much slower than
    # the rate of chain elongation (kx), it
    # should also be much slower than one over the
    # post-transcriptional simulation time --t-end
    # NOTE: k-slow is not supported anymore. It was 
    # primarily motivated by numeric instabilities 
    # during simulations, caused by large differences 
    # in rate constants. Those seem to be resolved.
    #
    # Parameters:
    # k0 = maximum folding rate /s
    # t-ext = time for chain elongation
    # t-end = post-transcriptional simulation time
    # and a few more ...
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    if args.minh:
        logger.warning('Overwriting t-fast parameter.')
        args.t_fast = 1/(args.k0 * math.exp(-args.minh/_RT))
    else:
        args.minh = max(0, -_RT * math.log(1 / args.t_fast / args.k0))
    logger.info(f'--t-fast: {args.t_fast} s => {args.minh} kcal/mol barrier height ' + 
                f'and {1/args.t_fast} /s rate at k0 = {args.k0}')

    if not args.force and args.t_fast and args.t_fast * 10 > args.t_ext:
        raise SystemExit('ERROR: Conflicting Settings: ' + 
                'Arguments must be such that "--t-fast" * 10 > "--t-ext".\n' + 
                '       => An instant folding time must be at least 10x shorter than ' +
                'the time of nucleotide extension. You may use --force to ignore this setting.')

    ############################
    # ~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Start with DrTransformer #
    # ~~~~~~~~~~~~~~~~~~~~~~~~ #
    ############################

    if args.paramFile:
        RNA.read_parameter_file(args.paramFile)

    # Set model details.
    vrna_md = RNA.md()
    vrna_md.noLP = 1
    vrna_md.logML = 0
    vrna_md.temperature = args.temp
    vrna_md.dangles = args.dangles
    vrna_md.special_hp = not args.noTetra
    vrna_md.noGU = args.noGU
    vrna_md.noGUclosure = args.noClosingGU

    # Write logging output
    if args.stdout == 'log' or lfh:
        fdata  = "# File generated using DrTransformer v{}\n".format(__version__)
        fdata += "#\n"
        fdata += "# >{}\n# {} \n".format(name, fullseq)
        fdata += "#\n"
        fdata += "# Co-transcriptional folding parameters:\n"
        fdata += "# --t-ext: {} sec\n".format(args.t_ext)
        fdata += "# --t-end: {} sec\n".format(args.t_end)
        fdata += "# --start: {}\n".format(args.start)
        fdata += "# --stop: {}\n".format(args.stop)
        fdata += "#\n"
        fdata += "# Algorithm parameters:\n"
        fdata += "# --o-prune: {}\n".format(args.o_prune)
        fdata += "# --t-fast: {} sec\n".format(args.t_fast)
        fdata += "# --fpwm: {}\n".format(args.fpwm)
        fdata += "# --mfree: {} nuc\n".format(args.mfree)
        fdata += "# --k0: {}\n".format(args.k0)
        fdata += "#\n"
        fdata += "# ViennaRNA model details:\n"
        fdata += "# --temp: {} C\n".format(args.temp)
        fdata += "# --dangles: {}\n".format(args.dangles)
        fdata += "# --paramFile: {}\n".format(args.paramFile)
        fdata += "# --noTetra: {}\n".format(args.noTetra)
        fdata += "# --noGU: {}\n".format(args.noGU)
        fdata += "# --noClosingGU: {}\n".format(args.noClosingGU)
        fdata += "#\n"
        fdata += "#\n"
        fdata += "# Results:\n"
        fdata += "# Transcription Step | Energy-sorted structure count | Structure | Energy "
        fdata += "| [Occupancy-t0 Occupancy-t8] | Structure ID (-> Plotting ID)\n"
        write_output(fdata, stdout=(args.stdout == 'log'), fh = lfh)

    # Write DrForna output
    if args.stdout == 'drf' or dfh:
        # Dictionary to store time course data.
        all_courses = dict()
        fdata = "id time occupancy structure energy\n"
        write_output(fdata, stdout = (args.stdout == 'drf'), fh = dfh)

    # Initialize a directed conformation graph
    TL = TrafoLandscape(fullseq, vrna_md)
    TL.k0 = args.k0
    TL.fpwm = args.fpwm
    TL.mfree = args.mfree
    TL.minh = int(round(args.minh*100)) if args.minh else 0
    TL.transcript_length = args.start - 1

    psites = np.full(args.stop+1, args.t_ext, dtype = float)
    #NOTE: psites[0] is useless
    if args.pause_sites:
        for term in args.pause_sites:
            site, pause = term.split('=')
            psites[int(site)] = float(pause)

    #############
    # ~~~~~~~~~ #
    # Main loop #
    # ~~~~~~~~~ #
    #############

    if args.performance_report:
        import statprof
        statprof.start()

    time = 0
    for tlen in range(args.start, args.stop+1):
        logger.info(f'** Transcription step {tlen} **')
        logger.info(f'Before expansion:      {len(list(TL.active_local_mins)):3d} active lmins, {len(list(TL.active_nodes)):3d} active structures, {len(list(TL.inactive_nodes)):3d} inactive structures.')
        itime = datetime.now()

        # Get new nodes and connect them.
        nn, on, prep = TL.expand(args.performance_report)
        logger.info(f'After expansion:                         {len(list(TL.active_nodes)):3d} active structures, {len(list(TL.inactive_nodes)):3d} inactive structures.' +
                    f' (Found {len(nn)} new nodes and revisited {len(on)} pruned nodes.)')
        etime = datetime.now()

        cn, ce = TL.get_coarse_network()
        logger.info(f'After coarse graining: {len(list(TL.active_local_mins)):3d} active lmins, {len(list(TL.active_nodes)):3d} active structures, {len(list(TL.inactive_nodes)):3d} inactive structures.' +
                    f' (Simulation network size: nodes = {cn}, edges = {ce}.)')
        ctime = datetime.now()

        # Adjust the length of the lin-time simulation:
        t0, t1 = 0, psites[tlen]
        # Adjust the length of the log-time simulation:
        t8 = t1 + args.t_end if tlen == args.stop else t1 + sum(psites[tlen+1:])
        if np.isclose(t0, t1) or np.isclose(t1, t8): # only lin or log-part!
            if np.isclose(t1, t8):
                times = np.array(np.linspace(t0, t1, args.t_lin))
            else:
                times = np.array(np.logspace(np.log10(t1), np.log10(t8), num=args.t_log))
        else:
            lin_times = np.array(np.linspace(t0, t1, args.t_lin))
            log_times = np.array(np.logspace(np.log10(t1), np.log10(t8), num=args.t_log))
            log_times = np.delete(log_times, 0)
            times = np.concatenate([lin_times, log_times])
        if tlen != args.start:
            times = np.delete(times, 0)

        snodes, p0 = TL.get_occupancies()
        assert np.isclose(sum(p0), 1)
        if args.plot_minh:
            assert args.plot_minh > TL.minh, "Plot-minh must be greater than minh."
            # Provide coarse network to get even more coarse network
            ndata = {n: d for n, d in TL.nodes.items() if n in snodes} # only active.
            edata = TL.cg_edges
            assert (n in snodes for n in ndata)
            assert (n in ndata for n in snodes)
            _, _, plot_cgm = top_down_coarse_graining(ndata, edata, args.plot_minh)

            mapping = {hn: set() for hn in snodes}
            for lmin, hidden in plot_cgm.items():
                assert lmin in snodes
                for hn in hidden:
                    assert hn in snodes
                    mapping[hn].add(lmin)
        else:
            plot_cgm, mapping = None, None

        if args.tmpdir:
            _fname = _tmpdir + '/' + name + '-' + str(tlen)
            # Produce input for treekin simulation
            nlist, bbfile, brfile, bofile, callfile = get_tkn_simulation_files(TL, _fname) 

        ti, p8, pf = 0, np.zeros(len(p0)), np.zeros(len(p0))
        for (t, pt) in TL.simulate(snodes, p0, times, force=[t1]):
            if tlen < args.stop and t > t1:
                # NOTE: this is extra time to determine whether structures will
                # become important in the time frame of transcription.
                pf = np.maximum(pf, pt)
                continue
            if args.stdout == 'drf' or dfh:
                tt = time + t
                for e, occu in enumerate(pt):
                    node = snodes[e]
                    if plot_cgm and node not in plot_cgm:
                        continue
                    elif plot_cgm:
                        for n in plot_cgm[node]:
                            ci = snodes.index(n)
                            occu += pt[ci]/len(mapping[n])
                    ni = TL.nodes[node]['identity']
                    if occu < 0.001 and ni not in all_courses:
                        continue
                    ne = TL.nodes[node]['energy']/100
                    all_courses[ni] = [(tt, occu)]
                    fdata = f"{ni} {tt:03.9f} {occu:03.4f} {node[:tlen]} {ne:6.2f}\n"
                    write_output(fdata, stdout = (args.stdout == 'drf'), fh = dfh)
            ti, p8 = t, pt
        keep = [n for e, n in enumerate(snodes) if pf[e] > args.o_prune] if args.o_prune else []
        TL.set_occupancies(snodes, p8)
        time += ti

        # NOTE: just for debugging, check if p8 is calculated correctly.
        #tlp8 = TL.get_equilibrium_occupancies(snodes)
        #Z = sum(math.e**(-TL.nodes[n]['energy']/100/_RT) for n in snodes)
        #myp8 = np.array([math.e**(-TL.nodes[n]['energy']/100/_RT)/Z for n in snodes])
        #assert np.allclose(tlp8, myp8)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Output simulation results #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~ #

        # Print the current state *after* the simulation.
        if args.stdout == 'log' or lfh:
            for e, node in enumerate(snodes):
                ni = TL.nodes[node]['identity']
                ne = TL.nodes[node]['energy']/100
                po = p0[e]
                no = p8[e]
                if args.plot_minh:
                    lmins = [TL.nodes[n]['identity'] for n in sorted(mapping[node], 
                                                                     key = lambda x: TL.nodes[x]['identity'])]
                    if lmins:
                        ni +=  f" -> {', '.join(lmins)}"
                ax = ' ' if np.isclose(po, no, atol=1e-4) else '+' if no > po else '-'
                fdata = f"{tlen:4d} {e+1:4d} {node[:tlen]} {ne:6.2f} {ax}[{po:6.4f} -> {no:6.4f}] ID = {ni}\n"
                write_output(fdata, stdout = (args.stdout == 'log'), fh = lfh)
        stime = datetime.now()

        # ~~~~~ #
        # Prune #
        # ~~~~~ #
        if args.o_prune > 0:
            delth = args.delth
            pn, dn = TL.prune(args.o_prune, delth, keep) 
            logger.info(f'After pruning:         {len(list(TL.active_local_mins)):3d} active lmins,' + 
                        f' {len(list(TL.active_nodes)):3d} active structures, {len(list(TL.inactive_nodes)):3d} inactive structures.' +
                        f' (Pruned {len(pn)} nodes, kept {len(keep)} nodes due to lookahead, deleted {len(dn)} inactive nodes.)')

        ptime = datetime.now()
        exptime = (etime - itime).total_seconds() 
        cgntime = (ctime - etime).total_seconds()
        simtime = (stime - ctime).total_seconds()
        prntime = (ptime - stime).total_seconds()
        tottime = (ptime - itime).total_seconds()
        logger.info(f'{tlen=}, {tottime=}, {exptime=}, {cgntime=}, {simtime=}, {prntime=}.')
        if args.performance_report:
            (exptime, frayytime, guidetime, floodtime) = prep
            print(tlen, tottime, frayytime, guidetime, floodtime, exptime, cgntime, simtime, prntime)
            sys.stdout.flush()

    # Write the last results
    if args.stdout == 'log' or lfh:
        pe = TL.get_equilibrium_occupancies(snodes)
        fdata  = "# Distribution of structures at the end:\n"
        fdata += "#         {}\n".format(TL.transcript)
        lnodes, pX = TL.get_occupancies()
        if args.plot_minh:
            for e, node in enumerate([s for s in snodes if s in plot_cgm]):
                ne = TL.nodes[node]['energy']/100
                no = p8[e] + sum(p8[snodes.index(n)]/len(mapping[n]) for n in plot_cgm[node])
                eo = pe[e] + sum(pe[snodes.index(n)]/len(mapping[n]) for n in plot_cgm[node])
                ni = TL.nodes[node]['identity']
                nids = [TL.nodes[n]['identity'] for n in sorted(plot_cgm[node], key = lambda x: TL.nodes[x]['identity'])]
                if nids:
                    ni +=  f" + {' + '.join(nids)}"
                ax = ' ' if np.isclose(no, eo, atol=1e-4) else '+' if eo > no else '-'
                fdata += f"{tlen:4d} {e+1:4d} {node[:tlen]} {ne:6.2f} {no:6.4f} {ax}[t8: {eo:6.4f}] ID = {ni}\n"
            write_output(fdata, stdout = (args.stdout == 'log'), fh = lfh)
        else:
            for e, node in enumerate(snodes):
                ne = TL.nodes[node]['energy']/100
                no = p8[e]
                eo = pe[e]
                ni = TL.nodes[node]['identity']
                ax = ' ' if np.isclose(no, eo, atol=1e-4) else '+' if eo > no else '-'
                fdata += f"{tlen:4d} {e+1:4d} {node[:tlen]} {ne:6.2f} {no:6.4f} {ax}[t8: {eo:6.4f}] ID = {ni}\n"
            write_output(fdata, stdout=(args.stdout == 'log'), fh = lfh)

    # CLEANUP file handles
    if lfh: lfh.close()
    if dfh: dfh.close()

    if args.performance_report:
        statprof.stop()
        statprof.display()

    return

if __name__ == '__main__':
    main()
