#
# drtransformer.linalg
# 
# Calculate matrix exponentials using numpy.
#
import logging

import argparse
import numpy as np
import scipy
import scipy.linalg as sl

def get_p8_detbal(A):
    """ Calculate the equilibrium distribution vector p8.

    Given a the rate matrix A, calculate the equilibrium distribution using
    detailed balance. 
    WARNING: This function requires a matrix with ZERO entries. Do not use it
    with numerical result matrices where some A[i][j] may be "almost" zero.

    A matrix entry A[i][j] must be the rate i->j. 
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
                assert A[i][j] > 0 and A[j][i] > 0
                if p8[j] > 0:
                    p8[j] += p8[i] * (A[i][j] / A[j][i])
                    p8[j] /= 2
                else:
                    p8[j] = p8[i] * (A[i][j] / A[j][i])
                    count += 1
                assert 0 <= p8[j] <= p8[0] # or np.isclose(p8[j], 0) or np.isclose(p8[j], 1)
        # Given i, we have looped through all j
        done[i] = 1
        # Determine next i for solving
        for c, p in enumerate(p8):
            if not done[c] and p > 0:
                i = c
                break
        assert not (i == dim and count < dim), "non-ergodic chain!"
    return p8 / sum(p8) # make the vector sum = 1.0

def symmetrize(A, p8):
    """
    Takes the transposed matrix: rates i <- j

    U = _sqrPI * A * sprPI_
    """
    dim = len(A)
    _sqrP8 = np.diag(1/np.sqrt(p8))
    sqrP8_ = np.diag(np.sqrt(p8))
    U = np.matmul(np.matmul(_sqrP8, A), sqrP8_)

    # correct for numerical errors
    err = 0
    for (i, j) in zip(*np.triu_indices(dim)):
        if i == j:
            continue
        err += abs((U[i][j] + U[j][i]) / 2 - U[i][j])
        U[i][j] = (U[i][j] + U[j][i]) / 2
        U[j][i] = U[i][j]

    #print(f'Corrected numerical error: {err} ({err/(dim*dim-dim)} per number)')
    return _sqrP8, U, sqrP8_

def decompose_sym(U):
    # Check if matrix has issues
    #evals, evecs = np.linalg.eig(U)
    evals, evecs = sl.eigh(U)

    # Sort them by largest eigenvalue
    idx = evals.argsort()[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]

    return evecs, evals, evecs.T # eves.T == np.linalg.inv(evecs) due to symmetry

def mxprint(A):
    return '\n'.join(' '.join([f'{x:10.4f}' for x in row]) for row in A)


def main():
    """ A python implementation of (some parts of) the treekin program.

    https://www.tbi.univie.ac.at/RNA/Treekin/
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='%(prog)s [options] filename')

    parser.add_argument("-v", "--verbose", action='count', default=0,
            help="""Track process by writing verbose output to STDOUT during calculations.""")

    parser.add_argument('input_filename', default=None, nargs='?', metavar='<str>', 
            help="Path to the input file.")

    parser.add_argument("--p0", nargs='+', metavar='<int/str>=<flt>', required=True,
            help="""Vector of initial species concentrations. 
            E.g. \"--p0 1=0.5 3=0.7\" stands for 1st species at a concentration of 0.5 
            and 3rd species at a concentration of 0.7. You may chose to address species
            directly by name, e.g.: --p0 C=0.5.""")
 
    parser.add_argument("--t-lin", type = int, default = None, metavar = '<int>',
            help = """Evenly space output *t-lin* times [t0, t1] at start of simulation.""")

    parser.add_argument("--t-log", type = int, default = None, metavar = '<int>',
            help = """Evenly space output *t-log* times (t1, t8] at end of simulation.""")

    parser.add_argument("--t0", type = float, default = 0, metavar = '<flt>',
            help = """Give t0.""")

    parser.add_argument("--t1", type = float, default = None, metavar = '<flt>',
            help = """Give t1.""")

    parser.add_argument("--t8", type = float, default = 10, metavar = '<flt>',
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
    A = np.loadtxt(args.input_filename, dtype=np.float64)
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
    logger.info(f"Initial occupancy vector ({sum(p0)=}):\n{mxprint([p0])}")
    logger.info(f"Input matrix A including diagonal elements:\n{mxprint(A)}")

    # Get equilibrium distribution.
    p8 = get_p8_detbal(A)
    logger.info(f"Equilibrium distribution vector p8 ({sum(p8)=}):\n{mxprint([p8])}")

    logger.info(f"\nSolving: U = _P8 * A^T * P8_\n")
    _P, U, P_= symmetrize(A.T, p8)
    logger.info("Symmetrized matrix U:\n" + mxprint(U))
    logger.info("1 / sqrt(P8) matrix:\n" + mxprint(_P))
    logger.info("sqrt(P8):\n" + mxprint(P_))
    pU8 = get_p8_detbal(U)
    logger.info(f"Equilibrium distribution vector pU8 ({sum(pU8)=}):\n{mxprint([pU8])}")

    logger.info(f"\nDecomposing: U = _S * L * S_\n")
    _S, L, S_ = decompose_sym(U)
    logger.info("Eigenvalues L:\n" + mxprint([L]))
    logger.info("Eigenvectors _S:\n" + mxprint(_S))
    logger.info("Eigenvectors S_:\n" + mxprint(S_)) # It's just the transpose of _S!
    
    CL = np.matmul(P_, _S)
    logger.info("Left correction matrix CL:\n" + mxprint(CL))
    CR = np.matmul(S_, _P)
    logger.info("Right correction matrix CR:\n" + mxprint(CR))

    logger.info("CL * CR:\n" + mxprint(np.matmul(CL, CR)))
    assert np.allclose(np.matmul(CL, CR), np.identity(dim), atol = 1e-4, rtol = 1e-4)

    tV = CR.dot(p0)
    logger.info(f"Correction vector tV ({sum(tV)=}):\n{mxprint([tV])}")
    assert np.allclose(tV, S_.dot(_P.dot(p0)))
    assert isinstance(tV[0], np.float64)

    # Simulation of decomposed matrix.
    for t in np.concatenate([lin_times, log_times]):
        eLt = np.diag(np.exp(L*t))
        pt = CL.dot(eLt.dot(tV))
        pt = CL.dot(np.multiply(np.exp(L*t), tV))
        print(f"{t:22.20e} {' '.join([f'{x:8.6e}' for x in pt])}")
        if np.allclose(pt, p8, rtol = 1e-4, atol = 1e-4):
            print(f'# Reached equilibrium at time {t=}.')
            #break
        if not np.isclose(sum(pt), 1., atol = 1e-4):
            print(f'# Unstable simulation at time {t=}: {sum(pt)=}')
            break

if __name__ == '__main__':
    main()

