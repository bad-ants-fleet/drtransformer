#!/usr/bin/env python

import numpy as np
import scipy
import unittest
from itertools import combinations

from drtransformer.linalg import (get_p8_detbal, 
                                  mx_print,
                                  mx_symmetrize,
                                  mx_decompose_sym)

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class MatrixOperations(unittest.TestCase):
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
        print(get_p8_detbal(A))
        assert np.allclose(p8, get_p8_detbal(A))
        # Test symmetrization
        myU = np.array([[ -57.25      ,   73.56573251,   97.54137584],
                        [  73.56573251, -124.9       ,    0.        ],
                        [  97.54137584,    0.        , -683.5       ]])

        print(A)
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
