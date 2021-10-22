# Call

import os

from struct import pack

def get_tkn_simulation_files(TL, basename):
    """ Print a rate matrix and the initial occupancy vector.

    This function prints files and parameters to simulate dynamics using the
    commandline tool treekin. A *.bar file contains a sorted list of present
    structures, their energy and their neighborhood and the corresponding
    energy barriers. A *.rts or *.rts.bin file contains the matrix of
    transition rates either in text or binary format. Additionaly, it returns
    a vector "p0", which contains the present occupancy of structures. The
    order or elements in p0 contains

    Note:
      A *.bar file contains the energy barriers to all other nodes in the
      graph. So it is not the same as a "classic" barfile produce by
      barriers.

    Args:
      basename (str): Basename of output files.

    Returns:

    """
    seq = TL.transcript
    nodes = TL.nodes
    snodes = sorted(TL.active_local_mins, key = lambda x: TL.nodes[x]['energy'])
    num = len(snodes) + 1

    bofile = basename + '_lands.bar'
    brfile = basename + '_rates.txt'
    bbfile = basename + '_rates.bin'
    p0file = basename + '_p0.txt'
    p0 = []

    with open(bofile, 'w') as bar, open(brfile, 'w') as rts, open(bbfile, 'wb') as brts, open(p0file, 'w') as pf:
        bar.write("  ID {}  Energy  {}\n".format(seq, 
            ' '.join(map("{:7d}".format, range(1, num)))))
        brts.write(pack("i", len(snodes)))
        for ni, node in enumerate(snodes, 1):
            ne = nodes[node]['energy']
            no = nodes[node]['occupancy']

            # Calculate barrier heights to all other basins.
            barstr = ''
            for other in snodes:
                oe = nodes[other]['energy']
                sE = TL.get_cg_saddle(node, other)
                if sE is not None:
                    barstr += ' {:7.2f}'.format((sE - ne)/100)
                else:
                    barstr += ' {:7.2f}'.format(float('nan'))

            # Print structures and neighbors to bfile:
            bar.write("{:4d} {} {:7.2f} {}\n".format(ni, node[:len(seq)], ne/100, barstr))

            # Add ni occupancy to p0
            if no > 0:
                p0.append("{}={}".format(ni, no))

            # Print rate matrix to rfile and brfile
            trates = []
            rates = []
            for other in snodes:
                if TL.has_cg_edge(node, other):
                    rates.append(TL.cg_edges[(node, other)]['weight'])
                    trates.append(TL.cg_edges[(other, node)]['weight'])
                else:
                    rates.append(0)
                    trates.append(0)
            line = "".join(map("{:10.4g}".format, rates))
            rts.write("{}\n".format(line))
            for r in trates:
                brts.write(pack("d", r))
        pf.write(f"cat {brfile} | treekin {' '.join(f'--p0 {p}' for p in p0)}\n")

    return snodes, bbfile, brfile, bofile, p0file

