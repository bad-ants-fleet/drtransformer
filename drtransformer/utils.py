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

def get_simulation_files(TL, basename, treekin = False):
    """Print a simulation files for rate matrix, the initial occupancy vector and a call for DrSimulate.

    NOTE: the rate matrix is given in transposed form if you select treekin output!

    This function prints files and parameters to simulate dynamics using the
    commandline tool DrSimulate (or treekin). 
     - A *_lands.txt file contains a sorted list of present structures,
        their energy and their neighborhood with corresponding energy barriers. 
     - A *_rates.txt file contains the matrix of transition rates. 
     - A *_pinit.txt file contains the vector "p0" with initial conditions for the
        rate matrix.
     - A *_drsim.txt contains the commandline-call quivalent.

    Args:
        basename (str): Basename of output files.
        treekin (bool, optional): Return output for treekin.

    Returns:
        A list of file names.
    """
    seq = TL.transcript
    nodes = TL.nodes
    snodes = sorted(TL.active_local_mins, key = lambda x: TL.nodes[x]['energy'])
    num = len(snodes) + 1

    lfile = basename + '_lands.txt'
    rfile = basename + '_rates.txt'
    pfile = basename + '_pinit.txt'
    sfile = basename + '_drsim.txt'

    p0 = []
    with open(lfile, 'w') as bar, open(rfile, 'w') as rts, \
         open(pfile, 'w') as ini, open(sfile, 'w') as sim:
        bar.write("  ID {}  Energy  {}\n".format(seq, 
            ' '.join(map("{:7d}".format, range(1, num)))))
        for ni, node in enumerate(snodes, 1):
            ne = nodes[node]['energy']
            no = nodes[node]['occupancy']
            p0.append(no)

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

            # Print rate matrix to rfile and brfile
            rates = []
            trates = [] # treekin format
            for other in snodes:
                if TL.has_cg_edge(node, other):
                    rates.append(TL.cg_edges[(other, node)]['weight'])
                    trates.append(TL.cg_edges[(node, other)]['weight'])
                else:
                    rates.append(0)
                    trates.append(0)
            if treekin:
                line = "".join(map("{:10.4g}".format, trates))
            else:
                line = "".join(map("{:10.4g}".format, rates))
            rts.write("{}\n".format(line))

        p0s = [f'{i:10.4g}' for i in p0]
        ini.write(f"{' '.join(p0s)}\n")
        if treekin:
            p0s = ' --p0 '.join(f'{i}={p}' for i, p in enumerate(p0, 1) if p > 0)
            sim.write(f"cat {rfile} | treekin --p0 {p0s}\n")
        else:
            p0s = ' '.join(f'{i}={p}' for i, p in enumerate(p0, 1) if p > 0)
            sim.write(f"DrSimulate -r {rfile} --p0 {p0s}\n")

    return snodes, lfile, rfile, pfile, sfile

