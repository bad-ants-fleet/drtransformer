#!/usr/bin/env python

import numpy as np
import scipy
import unittest
from itertools import combinations

from drtransformer.linalg import (get_p8_detbal, 
                                  mx_print,
                                  mx_simulate,
                                  mx_symmetrize,
                                  mx_decompose_sym)

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class MatrixOperations(unittest.TestCase):
    def test_simulation_SIR(self):
        # The model
        # S -> I k = .1
        # I -> R k = .3 
        # A hack for detailed balance.
        # I -> S k = .00001
        # R -> I k = .00001
        mx = """
            -0.1000   0.0001   0.0000
             0.1000  -0.3001   0.0001
             0.0000   0.3000  -0.0001
             """
        R = np.array([row.split() for row in mx.split('\n') if len(row.split())], dtype = np.float64) 
        p0 = np.zeros(len(R), dtype = np.float64)
        p0[1] = 1
        times = [0, 1, 10]
        out = list([t, *pt] for t, pt in mx_simulate(R, p0, times))
        exp = [[0, 3.1139235611686275e-19, 0.9999999999999998, 2.688530694733204e-16], 
               [1, 8.200619191284835e-05, 0.7407604171205738, 0.2591575766875133], 
               [10, 0.00015911946877757312, 0.050058708350295376, 0.949782172180927]]
        assert np.allclose(out, exp) # if this breaks, use atol and rtol

    def test_simulation_stockmarket(self):
        mx = """
             0.900 0.075 0.025
             0.150 0.800 0.050
             0.250 0.250 0.500 
             """
        R = np.array([row.split() for row in mx.split('\n') if len(row.split())], dtype = np.float64).T # --transpose!
        p0 = np.zeros(len(R), dtype = np.float64)
        p0[1] = 1
        times = [0, 1, 10]
        out = list([t, *pt] for t, pt in mx_simulate(R, p0, times))
        exp = [[0,  4.8971659961325386e-17,  1.0,                 3.849911819966555e-17], 
               [1,  0.1345978156175058,      0.8284029189628698,  0.03699926541962439], 
               [10, 0.5750044673082901,      0.35872364214983066, 0.06627189054187904]]
        assert np.allclose(out, exp) # if this breaks, use atol and rtol

    def test_get_p8_detbal_00(self):
        A = np.array([[-0.2,  0.2],
                      [ 0.2, -0.2]], dtype=np.float64)
        p8 = np.array([0.5, 0.5])
        assert np.allclose(p8, get_p8_detbal(A))

    def test_get_p8_detbal_01(self):
        mx = """
                    0  0.004171
             3.31e-06         0
             """
        A = np.array([row.split() for row in mx.split('\n') if len(row.split())], dtype = np.float64)
        p8 = np.array([9.99207055e-01, 7.92945421e-04])
        assert np.allclose(p8, get_p8_detbal(A))

    def test_mx_symmetrize_00(self):
        mx = """
                 0     124.9     683.5
             43.33         0         0
             13.92         0         0
        """
        A = np.array([row.split() for row in mx.split('\n') if len(row.split())], dtype = np.float64)
        np.fill_diagonal(A, -np.einsum('ij->j', A))
        dim = len(A)
        # Make sure p8 is correct.
        #p8 = np.array([0.74539718, 0.23942224, 0.01518058])
        p8 = np.array([0.73137732, 0.25372762, 0.01489506])
        #print(get_p8_detbal(A))
        assert np.allclose(p8, get_p8_detbal(A))
        # Test symmetrization
        myU = np.array([[ -57.25      ,   73.56573251,   97.54137584],
                        [  73.56573251, -124.9       ,    0.        ],
                        [  97.54137584,    0.        , -683.5       ]])

        #print(A)
        _, U, _ = mx_symmetrize(A, p8)
        # The weak check.
        assert np.allclose(myU, U)
        # The forced symmetrization check.
        assert all(U[i][j] == U[j][i] for (i, j) in combinations(range(dim), 2))

    def test_mx_decompose_sym_00(self):
        # NOTE: I think this test may break, since different eigenvectors may
        # be returned. At some point these test should contain examples
        # calculated by pen and paper.
        dim = 3
        U = np.array([[ 0, 76, 97],
                      [76,  0,  0],
                      [97,  0,  0]])
        assert all(U[i][j] == U[j][i] for (i, j) in combinations(range(dim), 2))

        eva = np.array([ 123.22743201, 0.        , -123.22743201])
        eve = np.array([[-0.70710678,  0.        ,  0.70710678],
                        [-0.43610513, -0.78716239, -0.43610513],
                        [-0.55660786,  0.61674579, -0.55660786]]), 

        evecs, evals, _ = mx_decompose_sym(U)
        assert np.allclose(eve, evecs)
        assert np.allclose(eva, evals)

if __name__ == '__main__':
    unittest.main()
