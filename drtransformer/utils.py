#
# drtransformer.utils
# 
# All sorts of useful stuff.
#

import logging
drlog = logging.getLogger(__name__)

import re
from struct import pack

class RiboUtilsError(Exception):
    pass

def make_pair_table(ss, base=0, chars=['.']):
    """ Return a secondary struture in form of pair table.

    Args:
      ss (str): secondary structure in dot-bracket format
      base (int, optional): choose between a pair-table with base 0 or 1
      chars (list, optional): a list of characters to be are ignored, default:
        ['.']

    **Example:**
       base=0: ((..)). => [5,4,-1,-1,1,0,-1]
        i.e. start counting from 0, unpaired = -1
       base=1: ((..)). => [7,6,5,0,0,2,1,0]
        i.e. start counting from 1, unpaired = 0, pt[0]=len(ss)

    Returns:
      [list]: A pair-table
    """
    stack = []
    if base == 0:
        pt = [-1] * len(ss)
    elif base == 1:
        pt = [0] * (len(ss) + base)
        pt[0] = len(ss)
    else:
        raise RiboUtilsError(f"unexpected value in make_pair_table: (base = {base})")

    for i, char in enumerate(ss, base):
        if (char == '('):
            stack.append(i)
        elif (char == ')'):
            try:
                j = stack.pop()
            except IndexError as e:
                raise RiboUtilsError("Too many closing brackets in secondary structure")
            pt[i] = j
            pt[j] = i
        elif (char not in set(chars)):
            raise RiboUtilsError(f"unexpected character in sequence: '{char}'")
    if stack != []:
        raise RiboUtilsError("Too many opening brackets in secondary structure")
    return pt

def make_loop_index(ss):
    """Returns a list of loop indices to see where new base-pairs can be formed.
    """
    loop, stack = [], []
    nl, L = 0, 0
    for i, char in enumerate(ss):
        if char == '(':
            nl += 1
            L = nl
            stack.append(i)
        loop.append(L)
        if char == ')':
            _ = stack.pop()
            try:
                L = loop[stack[-1]]
            except IndexError:
                L = 0
    return loop

def parse_vienna(filename, **kwargs):
    with open(filename, 'r') as fh:
        return parse_vienna_stdin(fh, **kwargs)

def parse_vienna_stdin(stdin, chars = 'ACUGN&', skip = '-'):
    """Parse name and sequence information from file with fasta format.

    Only one input-sequence is allowed at a time.

    Args:
      stdin (list): Input to parse, ususally :obj:`sys.stdin`
      chars (string, optional): Allowed characters in a sequence.

    Returns:
      str, str: name and sequence.
    """
    name = 'NoName'
    seq = ''
    for line in stdin:
        if re.match('>', line):
            if name != 'NoName':
                raise NotImplementedError(
                    'Only single-sequence fasta format supported!')
            else:
                name = line.strip().split()[0][1:]
        else:
            seq += line.strip()
    seq = seq.translate({ord(c): None for c in skip})
    m = re.search('[^' + chars + ']', seq)
    if m:
        raise RiboUtilsError("Does not look like RNA: ('{}' in '{}')".format(
            m.string[m.span()[0]], seq))
    return name, seq

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

