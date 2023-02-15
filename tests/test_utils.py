#!/usr/bin/env python

import unittest

from drtransformer.utils import parse_vienna_stdin


class test_io(unittest.TestCase):
    # fasta-like input parser that takes single name/sequence input
    def test_vienna_stdin(self):
        # single sequence input
        regular_input = ['>testseq with additional info in header',
                         'CUUGUUCGAUGUUGUCUUUCACGAAUCCGGUCCGCAACCCAUUAUAAGCAGAUGUGUAG']

        fasta_input = ['>testseq with additional info in header',
                       'CUUGUUCGAUGU ',
                       'UGUCUUUCACGAAUCCGGUCC',
                       ' GCAACCCAUUAUAAGCAGAUGUGUAG',
                       '']
        result = (
            'testseq',
            'CUUGUUCGAUGUUGUCUUUCACGAAUCCGGUCCGCAACCCAUUAUAAGCAGAUGUGUAG')
        self.assertEqual(parse_vienna_stdin(regular_input), result)
        self.assertEqual(parse_vienna_stdin(fasta_input), result)

    def test_vienna_stdin_cofold(self):
        # cofold input with multiple cutpoints
        multifold_input = ['>testseq with additional info in header',
                           'CUUGUUCGAUGU& ',
                           'UGUCUUUCACGAAUCCGGUCC&GCAACCCAUUAUAAGCAGAUGUGUAG']

        result = (
            'testseq',
            'CUUGUUCGAUGU&UGUCUUUCACGAAUCCGGUCC&GCAACCCAUUAUAAGCAGAUGUGUAG')
        self.assertEqual(parse_vienna_stdin(multifold_input, chars = 'ACGU&'), result)

if __name__ == '__main__':
    unittest.main()
