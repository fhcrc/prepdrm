import unittest

import pysam

from .. import samutil

OPNAMES = "MIDNSH"
OP_LOOKUP = {c: i for i, c in enumerate(OPNAMES)}
S = OP_LOOKUP['S']
M = OP_LOOKUP['M']
I = OP_LOOKUP['I']
D = OP_LOOKUP['D']

class AllPairsIterTestCase(unittest.TestCase):

    def test_simple(self):
        a = pysam.AlignedRead()
        a.seq = 'TAAAA'
        a.pos = 4
        a.cigarstring = 'M5'
        reference = 'AAAATAAGA'
        a.qual = ';' * len(a.seq)
        result = list(samutil.all_pairs_iter(a, reference))
        self.assertEqual(5, len(result))
        self.assertEqual(list('TAAGA'), [i.rbase for i in result])
        self.assertEqual(list(a.seq), [i.qbase for i in result])
        self.assertEqual([0] * 5, [i.cigar_op for i in result])
        self.assertEqual(range(a.pos, a.pos + 5), [i.rpos for i in result])
        self.assertEqual(range(5), [i.qpos for i in result])

    def test_more_ops(self):
        a = pysam.AlignedRead()
        #   ggTAA---Ca
        # AAAAC-AGTAA
        a.seq = 'GGTAACA'
        a.pos = 4
        a.cigarstring = 'S2M1I1M1D3M1S1'
        reference = 'AAAACAGTAA'
        a.qual = ';' * len(a.seq)
        ap = samutil.AlignedPair
        expected = [ap(0, None, 'G', None, ';', S),
                    ap(1, None, 'G', None, ';', S),
                    ap(2, 4, 'T', 'C', ';', M),
                    ap(3, None, 'A', None, ';', I),
                    ap(4, 5, 'A', 'A', ';', M),
                    ap(None, 6, None, 'G', None, D),
                    ap(None, 7, None, 'T', None, D),
                    ap(None, 8, None, 'A', None, D),
                    ap(5, 9, 'C', 'A', ';', M),
                    ap(6, None, 'A', None, ';', S),
                    ]
        result = list(samutil.all_pairs_iter(a, reference))
        self.assertEqual(expected, result)


def suite():
    suite = unittest.TestSuite()
    for c in [AllPairsIterTestCase]:
        suite.addTest(unittest.makeSuite(c))
    return suite
