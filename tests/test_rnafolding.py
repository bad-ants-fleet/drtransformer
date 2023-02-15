#!/usr/bin/env python

import RNA
import unittest

from drtransformer.rnafolding import (# FrayingNeighborsTests
                                      find_fraying_neighbors,
                                      open_fraying_helices,
                                      fold_exterior_loop,
                                      # ConstrainedFoldingTests
                                      get_basepairs, 
                                      mfe_intersect,
                                      # GuideGraphTests 
                                      get_bpd_cache, 
                                      get_bpd_i_cache,
                                      guiding_edge_search,
                                      shortcut_edge_search,
                                      guiding_node_search,
                                      get_guide_graph,
                                      # FindpathTests
                                      findpath_split,
                                      common_basepairs,
                                      common_exterior_bases, 
                                      split_struct, 
                                      merge_struct,
                                      # FloodingTests
                                      path_flooding, 
                                      edge_flooding,
                                      neighborhood_flooding,
                                      # TopDownCoarseGrainTests
                                      top_down_coarse_graining)

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class FrayingNeighborsTests(unittest.TestCase):
    def setUp(self):
        self.md = RNA.md()

    def test_find_fraying_neighbors_00(self):
        vmd = self.md
        seq =  "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        sss = ["......((((......)).))..(((((((....))))....)))..."]
        nbs = [".....(((((......)).)))((((((((....))))....)))).."]
        res = {sss[0]:set(nbs)}
        #print(find_fraying_neighbors(seq, vmd, sss, mfree = 6))
        assert res == find_fraying_neighbors(seq, vmd, sss, mfree = 6)

    def test_find_fraying_neighbors_01(self):
        vmd = self.md
        seq =  "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        sss = ["...((.(...........)))..(((((((....))))....)))..."]
        nbs = [".....(((((......)).)))((((((((....))))....))))..",
               "..(((.(...........))))((((((((....))))....)))).."]
        res = {sss[0]:set(nbs)}
        #print(find_fraying_neighbors(seq, vmd, sss, mfree = 6))
        assert res == find_fraying_neighbors(seq, vmd, sss, mfree = 6)

    def test_find_fraying_neighbors_02(self):
        vmd = self.md
        seq =  "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        sss = ["...((.(...........)))..(((((((....))))....)))...",
               "......((((......)).))..(((((((....))))....)))..."]
        nb0 = [".....(((((......)).)))((((((((....))))....))))..",
               "..(((.(...........))))((((((((....))))....)))).."]
        nb1 = [".....(((((......)).)))((((((((....))))....)))).."]
        res = {sss[0]:set(nb0),
               sss[1]:set(nb1)}
        #print(find_fraying_neighbors(seq, vmd, sss, mfree = 6))
        assert res == find_fraying_neighbors(seq, vmd, sss, mfree = 6)

    def test_open_fraying_helices_00(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        ss = ".....(((((......)).)))((((((((....))))....)))).."

        out = open_fraying_helices(se, ss, mfree=6)

        res = [
            "........((......))....((((((((....))))....))))..",
            ".....(((((......)).)))....((((....))))..........",
            "........((......))........((((....)))).........."]

        self.assertEqual(sorted(out), sorted(res))

        out = open_fraying_helices(se, ss, mfree=8)

        res = [
            "......................((((((((....))))....))))..",
            ".....(((((......)).)))....((((....))))..........",
            "..........................((((....)))).........."]

        self.assertEqual(sorted(out), sorted(res))

    def test_open_fraying_helices_multi(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCG"
        ss = "..((.(((((......)).)))((((((((....))))....)))).))"

        out = list(open_fraying_helices(se, ss, mfree=6))
        res = [".....(((((......)).)))((((((((....))))....))))..."]
        self.assertEqual(sorted(out), sorted(res))

        out = list(open_fraying_helices(se, ss, mfree=7))
        res = [
            "........((......))....((((((((....))))....))))...",
            ".....(((((......)).)))....((((....))))...........",
            "........((......))........((((....))))..........."]
        self.assertEqual(sorted(out), sorted(res))

    def test_fold_exterior_loop_00(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        ss = ".....(((((......)).)))((((((((....))))....)))).."
        md = self.md
        nbr, ecache = fold_exterior_loop(se, ss, md)

        self.assertEqual(nbr, '.....(((((......)).)))((((((((....))))....))))..')
        self.assertEqual(ecache, {'CUCGUCGNNNCGGGNNNCCGC': '.....((...))((...))..'})

@unittest.skipIf(SKIP, "skipping tests")
class ConstrainedFoldingTests(unittest.TestCase):
    def test_mfe_intersect(self):
        seq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCG"
        ss1 = "..((....(((........(((.....)))....)))..(((........)))..))." 
        ss2 = ".((((....))))..(((.(.((...)).)..)))......(((.......)))...."
        bps = get_basepairs([ss1, ss2])
        mss, mfe = mfe_intersect(seq, RNA.md(), bps)
        assert mss != ss1
        assert mss != ss2

    def test_mfe_intersect(self):
        seq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCG"
        ss1 = "..((....(((........(((.....)))....)))..(((........)))..))." 
        ss2 = ".((((....))))..(((.((((...))))..)))....(((........)))....."
        bps = get_basepairs([ss1, ss2])
        mss, mfe = mfe_intersect(seq, RNA.md(), bps)
        assert mss == ss2
 
    def test_mfe_intersect(self):
        seq = 'UCCGACAUUAAGACACACCAGGGUCUCGUAUCCCUAGGGUAAGGUACGCGCGGACCGGCCAAUCGGGUAUUGCUGCAAACUAUGGCAAUAGUGCAUAGGUUCAGACGAAGUACGGGUGGAUAUUUGUAGCCAGUAUGCUGGGUCUCCGGG'
        ss1 = '((((..........((.((((((........)))).))))..........))))((((...(((.(((((.((((((((((((.(((....)))))))(((((..((.....))..))))).))))))))....))))).)))..)))).'
        ss2 = '((((..........((.((((((........)))).))))..........))))((((....((((((((((((((((((((((((....)).)))))(((((..((.....))..))))).)))))))).))))).))))....)))).'
        ss3 = '((((..........((.((((((........)))).))))..........))))((((....((((((((((((((((((((((.((....)))))))(((((..((.....))..))))).)))))))).))))).))))....)))).'
        bps = get_basepairs([ss1, ss2, ss3])
        mss, mfe = mfe_intersect(seq, RNA.md(), bps)
        assert mss == ss3

@unittest.skipIf(SKIP, "skipping tests")
class GuideGraphTests(unittest.TestCase):
    def test_bpd_i_cache(self):
           p = '.((((.((((.((...(((...).))...))))))..))))...'
           i = '.((((.((((........(...)........))))..))))...'
           q = '.((((.((((.((...))(...).((...))))))..))))...'
           assert get_bpd_cache(p, i) == get_bpd_i_cache(p, q)
           assert get_bpd_cache(i, q) == get_bpd_i_cache(q, p)
           assert get_bpd_i_cache(p, i) == get_bpd_i_cache(p, q)
           assert get_bpd_i_cache(q, i) == get_bpd_i_cache(q, p)
           assert get_bpd_i_cache(i, p) == 0
           assert get_bpd_i_cache(i, q) == 0

    def test_guiding_edge_search_01(self):
        sss = ["..........((((....))))...",
               "...((((...))))...........",
               ".(((......)))............"]
        edges = guiding_edge_search(sss)
        assert len(edges) == 4
        edges = shortcut_edge_search(sss, edges)
        assert len(edges) == 6

    def test_guiding_edge_search_02(self):
        # NOTE(SB): I did not actually check if this is correct.
        """
              AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUGACAACAUU
            1 ..........((((.((((.((.((.......)).))))))..))))...  -6.70    0  13.00
            2 ..........((((.((((.((...((.....)).))))))..))))...  -6.10    1   2.10
            3 ..........((((.....((((.((.........)).)))).))))...  -5.90    1   6.30
            4 ((((.....(((........)))....))))....(((...)))......  -5.70    1   8.50
            5 ...((((...)))).....((((.((.........)).))))........  -5.60    3   4.80
            6 .(((......)))......((((.((.........)).))))........  -5.50    5   4.30
            7 ..........((((..((((.....((.....)).....))))))))...  -5.50    3   5.30
            8 ((((.....(((........)))....))))...................  -5.00    4   3.40
            9 ((((.....((.((.....))))....))))....(((...)))......  -5.00    4   2.80
           10 ((((.....((.((.....)).))...))))....(((...)))......  -4.90    4   3.50
        """
        sss = [
           "..........((((.((((.((.((.......)).))))))..))))...",
           "..........((((.((((.((...((.....)).))))))..))))...",
           "..........((((.....((((.((.........)).)))).))))...",
           "((((.....(((........)))....))))....(((...)))......",
           "...((((...)))).....((((.((.........)).))))........",
           ".(((......)))......((((.((.........)).))))........",
           "..........((((..((((.....((.....)).....))))))))...",
           "((((.....(((........)))....))))...................",
           "((((.....((.((.....))))....))))....(((...)))......",
           "((((.....((.((.....)).))...))))....(((...)))......"]
 
        sti = {x: e for e, x in enumerate(sss, 1)}
        edges = guiding_edge_search(sti.keys())
        #print()
        #for (x, y) in sorted(edges, key=lambda x: (sti[x[0]], sti[x[1]])):
        #    if sti[x] < sti[y]:
        #        print(sti[x], sti[y], x, y)
        assert len(edges) == 22

    def test_guiding_edge_search_03(self):
        seq = "AGACGACAAGGUUGAAUCGCA"
        sss = """(.((......)))........
                 .((((((...))))..))...
                 .(((......)))........
                 .((...((....))..))...
                 .(.(((((....))..)))).
                 .(.((((...))))..)....
                 .(.(((..........)))).
                 .(.((............))).
                 ...(((((....))..)))..
                 ...((((...)))).......
                 ...((((......)..)))..
                 ...(((..........)))..
                 ...((............))..
                 ....................."""
 
        sti = {s: e for e, s in enumerate(sss.split(), 1)}
        #print(len(sti))
        #for ss in sti:
        #    print(sti[ss], ss)
        edges = guiding_edge_search(set(sti.keys()))
        #print()
        #for (x, y) in sorted(edges, key=lambda x: (sti[x[0]], sti[x[1]])):
        #    if sti[x] < sti[y]:
        #        print(sti[x], sti[y], x, y)
        assert len(edges) == 36

    def test_guide_graph_construction_00(self):
        """
              AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUGACAACAUU
            #1 ..........((((.((((.((.((.......)).))))))..))))...  -6.70    0  13.00
            2 ..........((((.((((.((...((.....)).))))))..))))...  -6.10    1   2.10
            3 ..........((((.....((((.((.........)).)))).))))...  -5.90    1   6.30
            #4 ((((.....(((........)))....))))....(((...)))......  -5.70    1   8.50
            5 ...((((...)))).....((((.((.........)).))))........  -5.60    3   4.80
            6 .(((......)))......((((.((.........)).))))........  -5.50    5   4.30
            #7 ..........((((..((((.....((.....)).....))))))))...  -5.50    3   5.30
            8 ((((.....(((........)))....))))...................  -5.00    4   3.40
            #9 ((((.....((.((.....))))....))))....(((...)))......  -5.00    4   2.80
           10 ((((.....((.((.....)).))...))))....(((...)))......  -4.90    4   3.50
        """
        seq =   "AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUGACAACAUU"
        sss = [
          #( 1, "..........((((.((((.((.((.......)).))))))..))))...",  -6.70),
           ( 2, "..........((((.((((.((...((.....)).))))))..))))...",  -6.10),
           ( 3, "..........((((.....((((.((.........)).)))).))))...",  -5.90),
          #( 4, "((((.....(((........)))....))))....(((...)))......",  -5.70),
           ( 5, "...((((...)))).....((((.((.........)).))))........",  -5.60),
           ( 6, ".(((......)))......((((.((.........)).))))........",  -5.50),
          #( 7, "..........((((..((((.....((.....)).....))))))))...",  -5.50),
           ( 8, "((((.....(((........)))....))))...................",  -5.00),
          #( 9, "((((.....((.((.....))))....))))....(((...)))......",  -5.00),
           (10, "((((.....((.((.....)).))...))))....(((...)))......",  -4.90)]
 
        md = RNA.md()
        ndata = {ssx: {'energy': int(round(enx*100)), 'identity': idx} for (idx, ssx, enx) in sss}

        nodes, edges = get_guide_graph(seq, md, ndata.keys())
        assert all(n not in ndata for n in nodes)
        assert len(nodes) == 1 
        assert len(edges) == 18

    def test_guide_graph_construction_01(self):
        """
              CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCAACGUCCCUUGAUACCCCAGCCAACAAGCGG
              ((((.(((.........))))))).......................................................
              ((((.(((.........)))))))..((((....)))).........................................
              .....(((((...........)))))((((....)))).........................................
        """
        seq =  "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAG"
        sss = ["((((.(((.........)))))))..............",
               "((((.(((.........)))))))..((((....))))",
               ".....(((((...........)))))((((....))))"]
        md = RNA.md()
        nodes, edges = get_guide_graph(seq, RNA.md(), sss)
        assert len(nodes) == 0 # no new nodes
        assert len(edges) == 4 # two reversible edges
        assert (sss[0], sss[2]) not in edges

    def test_guide_graph_construction_02(self):
        seq =  "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAG"
        sss = ["((((.(((.........)))))))..............",
               "((((.(((.........)))))))..((((....))))",
               ".....(((((...........)))))((((....))))",
               "((.(.(((((...........)))))).))........"]
        md = RNA.md()
        nodes, edges = get_guide_graph(seq, RNA.md(), sss)
        assert len(nodes) == 0 # no new nodes
        assert (sss[0], sss[2]) in edges
        assert len(edges) == 10 # two reversible edges

    def test_guide_graph_construction_03(self):
        seq =    "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAG"
        sss =   ['((((.(((.........)))))))..............', # a 
                 '((((.(((.........)))))))..((((....))))', # b
                 '..((((((((......)).)))).))............', # c
                 '..((((((((......)).)))).))((((....))))', # d
                 '..((((((((......)).)).))))((((....))))', # e
                 '((.(.(((((...........)))))).))........', # f
                 '.....(((((...........)))))((((....))))'] # g

        neds = {('..((((((((......)).)).))))((((....))))', '((((.(((.........)))))))..((((....))))'), # e -> b 
                ('((((.(((.........)))))))..((((....))))', '..((((((((......)).)).))))((((....))))'), # b -> e 
                ('((((.(((.........)))))))..............', '((((.(((.........)))))))..((((....))))'), # a -> b 
                ('((((.(((.........)))))))..((((....))))', '((((.(((.........)))))))..............'), # b -> a 
                ('((((.(((.........)))))))..............', '..((((((((......)).)))).))............'), # a -> c 
                ('..((((((((......)).)))).))............', '((((.(((.........)))))))..............'), # c -> a 
                ('..((((((((......)).)))).))............', '..((((((((......)).)))).))((((....))))'), # c -> d 
                ('..((((((((......)).)))).))((((....))))', '..((((((((......)).)))).))............'), # d -> c 
                ('((((.(((.........)))))))..((((....))))', '..((((((((......)).)))).))((((....))))'), # b -> d 
                ('..((((((((......)).)))).))((((....))))', '((((.(((.........)))))))..((((....))))'), # d -> b 
                ('((((.(((.........)))))))..............', '((.(.(((((...........)))))).))........'), # a -> f 
                ('((.(.(((((...........)))))).))........', '((((.(((.........)))))))..............'), # f -> a 
                ('..((((((((......)).)))).))............', '((.(.(((((...........)))))).))........'), # c -> f 
                ('((.(.(((((...........)))))).))........', '..((((((((......)).)))).))............'), # f -> c 
                ('..((((((((......)).)))).))((((....))))', '..((((((((......)).)).))))((((....))))'), # d -> e
                ('..((((((((......)).)).))))((((....))))', '..((((((((......)).)))).))((((....))))'), # e -> d 
                ('((((.(((.........)))))))..............', '.....(((((...........)))))((((....))))'), # a -> g 
                ('.....(((((...........)))))((((....))))', '((((.(((.........)))))))..............'), # g -> a 
                ('((((.(((.........)))))))..((((....))))', '.....(((((...........)))))((((....))))'), # b -> g 
                ('.....(((((...........)))))((((....))))', '((((.(((.........)))))))..((((....))))'), # g -> b
                ('..((((((((......)).)))).))((((....))))', '.....(((((...........)))))((((....))))'), # d -> g 
                ('.....(((((...........)))))((((....))))', '..((((((((......)).)))).))((((....))))'), # g -> d 
                ('..((((((((......)).)).))))((((....))))', '.....(((((...........)))))((((....))))'), # e -> g 
                ('.....(((((...........)))))((((....))))', '..((((((((......)).)).))))((((....))))'), # g -> e 
                ('((.(.(((((...........)))))).))........', '.....(((((...........)))))((((....))))'), # f -> g 
                ('.....(((((...........)))))((((....))))', '((.(.(((((...........)))))).))........')} # g -> f 

        md = RNA.md()
        nodes, edges = get_guide_graph(seq, RNA.md(), sss)
        assert len(nodes) == 0
        assert neds == edges 
        assert len(edges) == 26

@unittest.skipIf(SKIP, "skipping tests")
class FindpathTests(unittest.TestCase):
    # NOTE: The path returned by findpath is sometimes different on ubuntu and macos!
    def test_split_merge_01(self):
        seq = 'GCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUU'
        ss1 = '(......((((((((......)))).......))))..)'
        ss2 = '((((.......((((......)))).......)))...)'
        pt1 = RNA.ptable(ss1)
        pt2 = RNA.ptable(ss2)

        ceb = []
        assert list(common_exterior_bases(pt1, pt2)) == ceb

        cbp = [(1, 39), (12, 25), (13, 24), (14, 23), (15, 22)]
        assert list(common_basepairs(pt1, pt2)) == cbp
        for x in common_basepairs(pt1, pt2):
            assert seq == merge_struct(*split_struct(seq, *x), *x)
            assert ss1 == merge_struct(*split_struct(ss1, *x), *x)
            assert ss2 == merge_struct(*split_struct(ss2, *x), *x)

    def test_split_merge_02(self):
        seq = 'AAAGCCGCCUUAAACGGGUAUUGGUACCNNNGGCAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUAC'
        ss1 = '...((((.(((....)))...)))).((...))...((((.(((..((.(((((......)))))..))..)))...))))..'
        ss2 = '....((((((.....))))...))..((...))....((((((.......((((......)))).......))))...))...'
        pt1 = RNA.ptable(ss1)
        pt2 = RNA.ptable(ss2)

        ceb = [1, 2, 3, 26, 34, 35, 36, 82, 83]
        assert list(common_exterior_bases(pt1, pt2)) == ceb
        for x in common_exterior_bases(pt1, pt2):
            assert seq == merge_struct(*split_struct(seq, x, None), x, None)
            assert ss1 == merge_struct(*split_struct(ss1, x, None), x, None)
            assert ss2 == merge_struct(*split_struct(ss2, x, None), x, None)

        cbp = [(5, 24), (6, 23), (27, 33), (28, 32), (38, 80), (39, 79), (51, 64), (52, 63), (53, 62), (54, 61)]
        assert list(common_basepairs(pt1, pt2)) == cbp
        for x in common_basepairs(pt1, pt2):
            assert seq == merge_struct(*split_struct(seq, *x), *x)
            assert ss1 == merge_struct(*split_struct(ss1, *x), *x)
            assert ss2 == merge_struct(*split_struct(ss2, *x), *x)

    def test_findpath_split_00(self):
        seq = 'GCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUU'
        ss1 = '(......((((((((......)))).......))))..)'
        ss2 = '((((.......((((......)))).......)))...)'
        # Set model details.
        path, barrier = findpath_split(seq, ss1, ss2, RNA.md())
        assert barrier == 810
        #assert path == [('(......((((((((......)))).......))))..)', 400),
        #                ('(.......(((((((......)))).......)))...)', 600),
        #                ('(.......((.((((......))))........))...)', 860),
        #                ('(........(.((((......))))........)....)', 1210),
        #                ('(..........((((......)))).............)', 680),
        #                ('((.........((((......)))).........)...)', 770),
        #                ('(((........((((......))))........))...)', 430),
        #                ('((((.......((((......)))).......)))...)', 280)]
        assert barrier == max((en for ss, en in path)) - path[0][1]

    def test_findpath_split_01(self):
        seq = 'AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUACAAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC'
        ss1 = '...((((.(((..((.(((((......)))))..))..)))...))))...............((((.(((..((.(((((......)))))..))..)))...))))............'
        ss2 = '...(((((((.......((((......)))).......))))...)))...............(((((((.......((((......)))).......))))...)))............'
        vrna_md = RNA.md()
        #print(f'\n{seq}\n{ss1}\n{ss2}')
        path, barrier = findpath_split(seq, ss1, ss2, vrna_md, th = 99)
        assert barrier == 410
        assert barrier == max((en for ss, en in path)) - path[0][1]
        #for (ss, en) in path:
        #    print(f'{ss} {en:>5d}')

        path, barrier = findpath_split(seq, ss1, ss2, vrna_md, th = 5)
        assert barrier == 410
        #for (ss, en) in path:
        #    print(f'{ss} {en:>5d}')

        path, barrier = findpath_split(seq, ss1, ss2, vrna_md, th = 1)
        assert barrier == 410
        #for (ss, en) in path:
        #    print(f'{ss} {en:>5d}')

    def test_findpath_split_02(self):
        seq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCG"
        ss1 = "..((....(((........(((.....)))....)))..(((........)))..))." 
        ss2 = ".((((....))))..(((.(.((...)).)..)))......(((.......)))...."
        path, barrier = findpath_split(seq, ss1, ss2, RNA.md(), th = 1)
        assert barrier == 200 # may change with a new findpath version?

@unittest.skipIf(True, "skipping tests")
class FloodingTests(unittest.TestCase):
    # TODO: needs tests.
    def test_path_flooding(self):
        pass

    def test_edge_flooding(self):
        pass

    def test_neighborhood_flooding(self):
        pass

@unittest.skipIf(SKIP, "skipping tests")
class TopDownCoarseGrainTests(unittest.TestCase):
    def test_top_down_coarse_graining_direct_path(self):
        seq = "UCUACUAUUCCGGCUUGACAUAAAUAUCGAGUGCUCGACCGCUAUUAUGGUACUUUCCAGCGUUUUGAUUGGUGGAUAAUAUCCCCCAAAAACGCGAGUC"
        path = [('............(((((..........)))))((((..((........)).........((((((...((((.((((...)))).)))))))))))))).', -1850), # lmin
                ('............(((((..........))))).(((..((........)).........((((((...((((.((((...)))).)))))))))))))..', -1710),
                ('............(((((..........)))))..((..((........)).........((((((...((((.((((...)))).))))))))))))...', -1440),
                ('............(((((..........)))))...(..((........)).........((((((...((((.((((...)))).)))))))))))....', -1300), # saddle
                ('............(((((..........)))))......((........)).........((((((...((((.((((...)))).)))))))))).....', -1660), # lmin
                ('............(((((..........))))).......(........)..........((((((...((((.((((...)))).)))))))))).....', -1370),
                ('............(((((..........)))))...........................((((((...((((.((((...)))).)))))))))).....', -1620),
                ('...........((((((..........))))).......)...................((((((...((((.((((...)))).)))))))))).....', -1230),
                ('..........(((((((..........))))).......))..................((((((...((((.((((...)))).)))))))))).....', -1390),
                ('..........((.((((..........))))........))..................((((((...((((.((((...)))).)))))))))).....', -1110), # saddle
                ('..........(((((((..........)))).......)))..................((((((...((((.((((...)))).)))))))))).....', -1520), # lmin
                ('....(.....(((((((..........)))).......)))........).........((((((...((((.((((...)))).)))))))))).....', -1100), # saddle
                ('....((....(((((((..........)))).......))).......)).........((((((...((((.((((...)))).)))))))))).....', -1260),
                ('....(((...(((((((..........)))).......)))......))).........((((((...((((.((((...)))).)))))))))).....', -1380),
                ('....((((..(((((((..........)))).......))).....)))).........((((((...((((.((((...)))).)))))))))).....', -1580),
                ('...(((((..(((((((..........)))).......))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1720), # lmin
                ('...(((((..((((((............))).......))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1440),
                ('...(((((..(((((..............)).......))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1280),
                ('...(((((..((((................).......))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1140), # saddle
                ('...(((((..(((.........................))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1560),
                ('...(((((..(((...(..................)..))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1380),
                ('...(((((..(((..((..................)).))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1480),
                ('...(((((..(((.(((..................)))))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1750),
                ('...(((((..(((.((((................))))))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1870),
                ('...(((((..(((.(((((.............).))))))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1860),
                ('...(((((..(((.((((((...........)).))))))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1950),
                ('...(((((..(((.(((((((.........))).))))))).....)))))........((((((...((((.((((...)))).)))))))))).....', -2150),
                ('..((((((..(((.(((((((.........))).))))))).....)))))).......((((((...((((.((((...)))).)))))))))).....', -2270), # lmin
                ('..((((((..(((.(((((((.........))).))))))).....)))))).....(.((((((...((((.((((...)))).)))))))))))....', -2140),
                ('..((((((..(((.(((((((.........))).))))))).....)))))).....(..(((((...((((.((((...)))).))))))))).)....', -1710)]

        lmins = ['............(((((..........)))))((((..((........)).........((((((...((((.((((...)))).)))))))))))))).', 
                 '............(((((..........)))))......((........)).........((((((...((((.((((...)))).)))))))))).....',
                 '..........(((((((..........)))).......)))..................((((((...((((.((((...)))).)))))))))).....',
                 '...(((((..(((((((..........)))).......))).....)))))........((((((...((((.((((...)))).)))))))))).....',
                 '..((((((..(((.(((((((.........))).))))))).....)))))).......((((((...((((.((((...)))).)))))))))).....']

        ndata = {ss: {'identity': e, 'energy': en} for e, (ss, en) in enumerate(path)}
        gedges = guiding_edge_search(set(ndata.keys()))
        edata = dict()
        for (x, y) in gedges:
            assert get_bpd_cache(x, y) == 1
            edata[(x,y)] = {'saddle_energy': max(ndata[x]['energy'], ndata[y]['energy'])}

        cgn, cge, cgm = top_down_coarse_graining(ndata, edata, minh = 300)
        for n in sorted(cgn, key = lambda x: cgn[x]['identity']):
            assert n in lmins

        # sum of basins can be larger, whenever structures are saddle points
        assert sum(len(cgm[n]) for n in cgn) >= len(ndata)-len(cgn)

        for i, node in enumerate(sorted(cgn, key = lambda x: cgn[x]['identity'])):
            ne = cgn[node]['energy']
            # Calculate barrier heights to all other basins.
            barstr = ''
            for j, other in enumerate(sorted(cgn, key = lambda x: cgn[x]['identity'])):
                oe = cgn[other]['energy']
                sE = cge[(node, other)]['saddle_energy'] if (node, other) in cge else None
                if sE is not None:
                    barstr += ' {:7.2f}'.format((sE - ne)/100)
                    assert abs(i-j)==1
                else:
                    barstr += ' {:7.2f}'.format(float('nan'))
                    assert abs(i-j)!=1

        #from drtransformer.rnafolding import as_barfile
        #print('\n' + as_barfile(seq, cgn, cge, cgm))

    def test_example_workflow(self):
        # A random set of sequences returned by randseq -l 100 | RNAsubopt --stochBT_en=50 | sort -u
        seq =       "CAAAGGCGACUCUCCUUAGACUCUAUAAAUAGUAAAUAGCUCCUAGGGACAAGGCUUACGUCCGCGUUAUUCACAUAAGUCGUCGCUUCAGUUGUGCGAC"
        sss = [( 1, "............((((((((..((((.........)))).)).)))))).............................(((((.((.......)))))))", -10.40),
               ( 2, ".....((((((.(((((((...((.............))...)))))))...(((....)))...((.....))...)))))).................", -10.80),
               ( 3, "....((((((..((((((((..((((.........)))).)).))))))...((.(.(((....))).).))......))))))(((......).))...", -11.20),
               ( 4, "....(((((((.(((((((...((((.........))))...)))))))..(.((........)).)..........)))))))..(.(.((...)))).", -11.60),
               ( 5, ".....(((((..(((((((...(((((.....)..))))...)))))))...((((((.((..(.....)..)).)))))))))))..............", -13.80),
               ( 6, "...(((((((..((((((....(((...........)))....))))))...((((((.((...........)).)))))))))))))((....))....", -14.20),
               ( 7, "...(.(((.....((((((...((((.........))))...))))))....(((....)))))).)...............((((.........)))).", -14.60),
               ( 8, "....((((((...((((((...((((.........))))...))))))....((((((.((...........)).))))))))))))..(....).....", -15.00),
               ( 9, "..((.(((.(..(((((((...((((.........))))...)))))))...)((....)).))).)).............(((((.........)))))", -15.10),
               (10, "...((....))..((((((...((((.........))))...))))))....((((((.................))))))(((((.........)))))", -15.80),
               (11, "....((((((..(((((((...((((.........))))...)))))))....((........))((.....)).......)))))).............", -16.00),
               (12, "....((((((..(((((((...((((.........))))...)))))))...((((((.(.............).))))))))))))...((......))", -16.10),
               (13, "....(((((((.(((((((...((.............))...))))))).((.((........)).)).........)))))))((.........))...", -16.20),
               (14, "..........(.(((((((...((((.........))))...))))))).).((((((.(.............).))))))(((((.........)))))", -16.20),
               (15, "...((.(.....(((((((...((((.........))))...)))))))...).))......................((((.(((.......)))))))", -16.60),
               (16, "...(((((((...((((((...((((..(...)..))))...))))))....((((((.((...........)).)))))))))))))............", -17.00),
               (17, "....((((((...((((((...((((.........))))...))))))....((((((.((...........)).))))))))))))((.((...)))).", -17.00),
               (18, "((((((((((..(((((((...((((........).)))...)))))))...((((((.((...........)).)))))))))))))...)))......", -17.20),
               (19, "...(((((((..(((((((...((((.........))))...)))))))...((((((.((....(.....))).)))))))))))))((....))....", -17.30),
               (20, "...(((((((...(((((((..((((.........)))).).))))))....((((((.((...........)).)))))))))))))............", -17.40),
               (21, "((((((((((..(((((((...((((.........))))...)))))))...((((((.((..(......).)).)))))))))))))...)))......", -17.70),
               (22, "....((((((..((((((((..((((.........)))).).)))))))...((((((.((...........)).)))))))))))).............", -17.90),
               (23, ".....(((((...((((((...((((.........))))...))))))....((((((.((...........)).)))))))))))..............", -17.90),
               (24, "(..(((((((..(((((((...(((...........)))...)))))))...((((((.((...........)).)))))))))))))..).........", -18.10),
               (25, "....(((((((.(((((((...((((.........))))...)))))))..(.((........)).)..........)))))))((.........))...", -18.20),
               (26, "(..(((((((..(((((((...((((.........))))...)))))))...((((((((....))).........))))))))))))..).........", -18.30),
               (27, "...(((((((..((((((((..((((.........)))).).)))))))...((((((.((...........)).)))))))))))))............", -18.50),
               (28, "((((((((((..(((((((...((((.........))))...)))))))...((((((.(.............).)))))))))))))...)))......", -18.70),
               (29, "...((((.....(((((((...((((.........))))...)))))))....))))........................(((((.........)))))", -18.70),
               (30, "....((((((..(((((((...((((.........))))...)))))))...((((((.((...........)).))))))))))))...(((....)))", -18.70),
               (31, "...(((((((..(((((((...((((.........))))...)))))))...((((((.((...........)).)))))))))))))..((...))...", -18.90),
               (32, "....((((((..(((((((...((((.........))))...)))))))...((((((.((...........)).))))))))))))((........)).", -18.90),
               (33, ".....((((((.(((((((...((((.........))))...))))))).((.((........)).)).........))))))(((.........)))..", -18.90),
               (34, "....((((((...((((((...((((.........))))...))))))....((((((.((...........)).)))))))))))).............", -19.20),
               (35, "..((.(((....(((((((...((((.........))))...)))))))...(((....)))))).)).............(((((.........)))))", -19.30),
               (36, "...(((((((...((((((...((((.........))))...))))))....((((((.((...........)).)))))))))))))............", -19.80),
               (37, "....((((((..(((((((...((((.........))))...)))))))...((((((.((...........)).)))))))))))).............", -20.30),
               (38, "...(((((((..(((((((...((((.........))))...)))))))...((((((.((...........)).)))))))))))))............", -20.90),
               (39, "....(((((((.(((((((...((((.........))))...)))))))....((........))............)))))))................", -17.60)]
        md = RNA.md()
        fpwm = 4
        myminh = 300

        ndata = {ssx: {'energy': int(round(enx*100)), 'identity': idx} for (idx, ssx, enx) in sss}

        gnodes, gedges = get_guide_graph(seq, md, ndata.keys())
        assert len(gnodes) == 25
        assert len(gedges) == 240
        for nid, (ss, en) in enumerate(gnodes, 40):
            ndata[ss] = {'energy': en, 'identity': 40}

        # NOTE: Results differ between ubuntu and macos due to findpath results!
        ndata, edata = neighborhood_flooding((seq, md, fpwm), ndata, gedges, minh = myminh)
        cgn, cge, cgm = top_down_coarse_graining(ndata, edata, minh = myminh)
        assert len(set([x for val in cgm.values() for x in val])) == len(ndata) - len(cgn)
        #from drtransformer.rnafolding import as_barfile
        #print('\n' + as_barfile(seq, cgn, cge, cgm))

if __name__ == '__main__':
    unittest.main()
