#
# drtransformer.linalg
# 
# Calculate matrix exponentials using numpy.
# Most of this code is inpired by the program treekin:
#  https://www.tbi.univie.ac.at/RNA/Treekin/ 
#

import logging
drlog = logging.getLogger(__name__)

import numpy as np
import scipy
import scipy.linalg as sl
from itertools import combinations

def get_p8_detbal(A):
    """ Calculate the equilibrium distribution vector "p8".

    Given a the rate matrix A, calculate the equilibrium distribution vector
    "p8" using detailed balance: P_i k_{ij} = P_j k_{ji}.

    Args:
        A: A numpy matrix with entries A[i][j] containing the rate constant for
            reaction i->j. 

    WARNING: This function requires a matrix with ZERO entries. Do not use it
    with numerical result matrices where some A[i][j] may be "almost" zero.

    Returns:
        [np.array] The equilbrium distribution vector p8.
    """
    dim = len(A)
    p8 = np.zeros(dim, dtype=np.float64)
    p8[0] = 1.

    i, count = 0, 1  # current index, number of solved states
    done = np.zeros(dim, dtype=int)
    while (count != dim):
        # Rewrite entries of p8
        for j in range(dim):
            if i == j:
                continue
            if A[i][j] != 0 or A[j][i] != 0:
                assert A[i][j] > 0 and A[j][i] > 0, "unbalanced matrix"
                if p8[j] > 0:
                    p8[j] += p8[i] * (A[i][j] / A[j][i])
                    p8[j] /= 2
                else:
                    p8[j] = p8[i] * (A[i][j] / A[j][i])
                    count += 1
                assert 0 <= p8[j] <= p8[0], "unreasonable equilibrium occupancy in p8"
        # Given i, we have looped through all j
        done[i] = 1
        # Determine next i for solving
        for c, p in enumerate(p8):
            if not done[c] and p > 0:
                i = c
                break
        assert not (i == dim and count < dim), "non-ergodic chain!"
    return p8 / sum(p8) # make the vector sum = 1.0

def mx_symmetrize(A, p8):
    """ Symmetrize a rate matrix by using its equilibrium distribution.

    Generating a symmetric matrix U throught the following formula:
        U = _sqrtP8 * A * sqrtP8_ 
    where _sqrtP8 and sqrtP8 are diagonal matrices constructed from the 
    equilibrium distribution vector.

    Args:
        A (np.array): A numpy matrix with entries A[i][j] containing the rate
            constant for reaction j<-i. (Note that this is the transpose of the
            "original" rate matrix.)
        p8 (np.array): The equilibrium distribution vector.

    Returns:
        U, _sqrP8, sqrP8_
    """
    dim = len(A)
    _sqrtP8 = np.diag(1/np.sqrt(p8))
    sqrtP8_ = np.diag(np.sqrt(p8))
    U = np.matmul(np.matmul(_sqrtP8, A), sqrtP8_)

    # force correction of numerical errors
    err = 0
    for (i, j) in combinations(range(dim), 2):
        err += abs((U[i][j] + U[j][i]) / 2 - U[i][j])
        U[i][j] = (U[i][j] + U[j][i]) / 2
        U[j][i] = U[i][j]
    drlog.debug(f'Corrected numerical error: {err} ({err/(dim*dim-dim)} per number).')
    return U, _sqrtP8, sqrtP8_

def mx_decompose_sym(U):
    """ Decompose a symmetric matrix into eigenvectors and eigenvalues.

    Note, we use scipy.linalg.eigh here instead of numpy.linalg.eig because the
    former is much more robust for symmetric matrices. It also should be
    faster, but I have not verified that.

    Args:
        U (np.array): A numpy matrix with entries U[i][j] containing the
            symmetric flux of reactions i->j and j<-i. 

    Returns:
        evecs, evals, evecs.T: A sorted matrix of eigenvectors, and the
            corresponing vector of eigenvalues.  
    """
    # Check if matrix has issues
    evals, evecs = sl.eigh(U)

    # Sort them by largest eigenvalue
    idx = evals.argsort()[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]
    return evecs, evals, evecs.T # eves.T == np.linalg.inv(evecs) due to symmetry

def mx_print(A):
    """ A helper function to pretty print matrices.
    """
    return '\n'.join(' '.join([f'{x:10.4f}' for x in row]) for row in A)

def main():
    """ A python implementation of (some parts of) the treekin program.

    https://www.tbi.univie.ac.at/RNA/Treekin/
    """
    import sys
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='%(prog)s [options] filename')

    parser.add_argument("-v", "--verbose", action='count', default=0,
            help = """Track process using verbose output.""")

    parser.add_argument("-r", "--rate-matrix", default=None, nargs='?', metavar='<str>', 
            help = "Path to the input file containing the rate matrix (A[i][j] = i -> j).")

    parser.add_argument("--p0", nargs='+', metavar='<int/str>=<flt>', required=True,
            help="""Vector of initial species concentrations. 
            E.g. \"--p0 1=0.5 3=0.7\" stands for 1st species at a concentration of 0.5 
            and 3rd species at a concentration of 0.7. You may chose to address species
            directly by name, e.g.: --p0 C=0.5.""")
 
    parser.add_argument("--t-lin", type = int, default = 10, metavar = '<int>',
            help = """Evenly space output *t-lin* times [t0, t1] at start of simulation.""")

    parser.add_argument("--t-log", type = int, default = 10, metavar = '<int>',
            help = """Evenly space output *t-log* times (t1, t8] at end of simulation.""")

    parser.add_argument("--t0", type = float, default = 0, metavar = '<flt>',
            help = """Give t0.""")

    parser.add_argument("--t1", type = float, default = 0.1, metavar = '<flt>',
            help = """Give t1.""")

    parser.add_argument("--t8", type = float, default = 1e5, metavar = '<flt>',
            help = """Give t8.""")

    args = parser.parse_args()

    # ~~~~~~~~~~~~~
    # Logging Setup 
    # ~~~~~~~~~~~~~
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()

    if args.verbose == 0:
        ch.setLevel(logging.WARNING)
    elif args.verbose == 1:
        ch.setLevel(logging.INFO)
    elif args.verbose == 2:
        ch.setLevel(logging.DEBUG)
    elif args.verbose >= 3:
        ch.setLevel(logging.NOTSET)

    logger.addHandler(ch)

    # Adjust simulation output times.
    t0, t8 = args.t0, args.t8
    if args.t_log:
        assert args.t1, "--t-log can only be given in combination with --t1 <float> !"
    if args.t_lin:
        t1 = args.t1 if args.t1 else args.t8
        lin_times = np.array(np.linspace(t0, t1, args.t_lin, dtype='float128'))
    else:
        t1 = args.t1 if args.t1 else args.t0
        lin_times = np.array([t1], dtype='float128')
    if args.t_log:
        log_times = np.array(np.logspace(np.log10(t1), np.log10(t8), num=args.t_log, dtype='float128'))
        log_times = np.delete(log_times, 0)
    else:
        log_times = []

    # Parse matrix input.
    A = np.loadtxt(args.rate_matrix, dtype=np.float64)
    np.fill_diagonal(A, -np.einsum('ji->j', A))
    dim = len(A)

    # Parse p0 input.
    p0 = np.zeros(dim, dtype = np.float64)
    for term in args.p0:
        p, o = term.split('=')
        p0[int(p) - 1] = float(o)
    assert np.isclose(sum(p0), 1)

    # TODO: assert ergodicity

    #
    # Start the linalg stuff.
    #
    logger.info(f"Initial occupancy vector ({sum(p0)=}):\n{mx_print([p0])}")
    logger.info(f"Input matrix A including diagonal elements:\n{mx_print(A)}")

    # Get equilibrium distribution.
    p8 = get_p8_detbal(A)
    logger.info(f"Equilibrium distribution vector p8 ({sum(p8)=}):\n{mx_print([p8])}")

    logger.info(f"\nSolving: U = _P8 * A^T * P8_\n")
    U, _P, P_= mx_symmetrize(A.T, p8)
    logger.info("Symmetrized matrix U:\n" + mx_print(U))
    logger.info("1 / sqrt(P8) matrix:\n" + mx_print(_P))
    logger.info("sqrt(P8):\n" + mx_print(P_))
    pU8 = get_p8_detbal(U)
    logger.info(f"Equilibrium distribution vector pU8 ({sum(pU8)=}):\n{mx_print([pU8])}")

    logger.info(f"\nDecomposing: U = _S * L * S_\n")
    _S, L, S_ = mx_decompose_sym(U)
    logger.info("Eigenvalues L:\n" + mx_print([L]))
    logger.info("Eigenvectors _S:\n" + mx_print(_S))
    logger.info("Eigenvectors S_:\n" + mx_print(S_)) # It's just the transpose of _S!
    
    CL = np.matmul(P_, _S)
    logger.info("Left correction matrix CL:\n" + mx_print(CL))
    CR = np.matmul(S_, _P)
    logger.info("Right correction matrix CR:\n" + mx_print(CR))

    logger.info("CL * CR:\n" + mx_print(np.matmul(CL, CR)))
    assert np.allclose(np.matmul(CL, CR), np.identity(dim), atol = 1e-4, rtol = 1e-4)

    tV = CR.dot(p0)
    logger.info(f"Correction vector tV ({sum(tV)=}):\n{mx_print([tV])}")
    assert np.allclose(tV, S_.dot(_P.dot(p0)))
    assert isinstance(tV[0], np.float64)

    times = np.concatenate([lin_times, log_times])
    # Simulation of decomposed matrix.
    for t in times:
        eLt = np.diag(np.exp(L*t))
        pt = CL.dot(eLt.dot(tV))
        pt = CL.dot(np.multiply(np.exp(L*t), tV))
        print(f"{t:8.6e} {' '.join([f'{x:8.6e}' for x in abs(pt)])}")
        if np.allclose(pt, p8, rtol = 1e-4, atol = 1e-4):
            logger.info(f'# Reached equilibrium at time {t=}.')
            print(f"{times[-1]:8.6e} {' '.join([f'{x:8.6e}' for x in p8])}")
            break
        if not np.isclose(sum(pt), 1., atol = 1e-4):
            logger.error(f'# Unstable simulation at time {t=}: {sum(pt)=}')
            break

if __name__ == '__main__':
    main()
