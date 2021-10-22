#
# drtransformer.utils
# 
# All sorts of useful stuff.
#

import logging
rlog = logging.getLogger(__name__)
import re

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

