import unittest

import pysam

from ..scripts import classify_mutations
from ..scripts.classify_mutations import classify_read

class ClassifyReadTestCase(unittest.TestCase):
    def setUp(self):
        self.mutation = classify_mutations.MutationLocation(reference='MB2059_pol', name='K65R', codon_start1=953, wt='K', mutations=[])
        self.mutation.mutations.append(classify_mutations.Mutation(self.mutation, 'R'))

    def test_non_coding(self):
        a = pysam.AlignedRead()
        a.qname = 'H3QYW4L01AXO44'
        a.seq = 'AACAATACTACAATATTGACCAAGTAAAAAGAAAGACAGTAACGTAGTAGGAAGTAAATTAGTAGTATTCAGAGAGACTTAATAAGAGAACTACAAGTACTTCGTAGGGTAAGTTACAACGTAGGAACTACCACACTCCTGACGGGGCTAAAAAGAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATACTTTTCAGTTCCCTTATATGAAGAATTTAGAAAATACACCGCATTCACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAAGGATGGAAAGGATCACCGGCAGTATTCCAAAGTAGCATGACAAAAATCTTAGAACCTTTTAGAAAACAAAATCCAGAAATGGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAATAAAAATAGAGGAATTAAGGGGACACCTATTGAAGTGGGGATTTACCACACCAGACAAAAAGCATCAGAAAGAACCACCATTT'
        a.cigarstring = 'S1M26I1M12I3M33I3M1D2M9I3M3I3M25I2M33D1M173D1M198'
        a.pos = 926

        result = classify_mutations.classify_read(a, self.mutation)
        self.assertEqual('non_coding', result.classification)

        # Try without context
        result = classify_read(a, self.mutation, context=0)
        self.assertEqual('wt', result.classification)

        # Allow mutations
        result = classify_read(a, self.mutation, 3, allow_indels_in_context=True)
        self.assertEqual('wt', result.classification)

    def test_wt(self):
        a = pysam.AlignedRead()
        a.seq = 'TAAAAGAAA'
        a.pos = 949
        a.cigarstring = 'M9'
        result = classify_read(a, self.mutation)
        self.assertEqual('wt', result.classification)

    def test_mut(self):
        a = pysam.AlignedRead()
        a.seq = 'AGAAGAAAG'
        a.pos = 949
        a.cigarstring = 'M9'
        result = classify_read(a, self.mutation)
        self.assertEqual('mut', result.classification)

    def test_leading_gap(self):
        a = pysam.AlignedRead()
        a.seq = 'GGAAGAAAG'
        a.pos = 948
        a.cigarstring = 'M1D1M8'
        result = classify_read(a, self.mutation)
        self.assertEqual('non_coding', result.classification)
        self.assertEqual('-GAAGAAAG', result.sequence)
        self.assertEqual('-gaAGAaag', result.codon_sequence)
        self.assertTrue(result.any_indels)

    def test_leading_gaps(self):
        a = pysam.AlignedRead()
        a.seq = 'AGAAAG'
        a.pos = 952
        a.cigarstring = 'M6'
        result = classify_read(a, self.mutation)
        self.assertEqual('non_coding', result.classification)
        self.assertEqual('---AGAAAG', result.sequence)
        self.assertEqual('---AGAaag', result.codon_sequence)
        self.assertTrue(result.any_indels)

    def test_not_covered(self):
        a = pysam.AlignedRead()
        a.seq = 'AGAAGAAAG'
        a.pos = 11000
        a.cigarstring = 'M9'
        result = classify_read(a, self.mutation, False)
        self.assertIsNone(result)


class ClassifyInsertionTestCase(unittest.TestCase):

    def setUp(self):
        mutation = classify_mutations.MutationLocation(reference='MB2059_pol',
                                                       name='K65R',
                                                       codon_start1=953,
                                                       wt='K', mutations=[])
        mutation.mutations.append(classify_mutations.Insertion(mutation, '+2'))
        self.mutation = mutation

    def test_no_insertion(self):
        a = pysam.AlignedRead()
        a.seq = 'TAAAAGAAA'
        a.pos = 949
        a.cigarstring = 'M9'

        self.assertEqual("wt", classify_read(a, self.mutation).classification)

    def test_correct_insertion(self):
        a = pysam.AlignedRead()
        a.seq = 'TAAAAGTTTCCCAAAT'
        a.pos = 949
        a.cigarstring = 'M6I6M3S1'
        self.assertEqual("mut", classify_read(a, self.mutation).classification)

    def test_noncoding_insertion(self):
        a = pysam.AlignedRead()
        a.seq = 'TAAAAGTTTCCAAAT'
        a.pos = 949
        a.cigarstring = 'M6I5M3S1'
        self.assertEqual("non_coding", classify_read(a, self.mutation).classification)

    def test_insuff_codon_insertion(self):
        a = pysam.AlignedRead()
        a.seq = 'TAAAAGTTTAAAT'
        a.pos = 949
        a.cigarstring = 'M6I3M3S1'
        r = classify_read(a, self.mutation)
        self.assertEqual("non_coding", r.classification)


def suite():
    suite = unittest.TestSuite()
    for c in [ClassifyReadTestCase, ClassifyInsertionTestCase]:
        suite.addTest(unittest.makeSuite(c))
    return suite
