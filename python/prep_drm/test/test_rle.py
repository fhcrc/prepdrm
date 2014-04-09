import unittest
from .. import rle

def rle_iter(rles):
    for i in rles:
        yield i.c, i.length

class RLETestCase(unittest.TestCase):

    def _roundtrip(self, s, rle_version):
        r = rle.encode(s)
        self.assertEqual(rle_version, list(rle_iter(r)))
        self.assertEqual(s, rle.decode(r))

    def test_empty(self):
        self._roundtrip('', [])

    def test_no_lengths(self):
        self._roundtrip('aca', [('a', 1), ('c', 1), ('a', 1)])

    def test_single_run(self):
        self._roundtrip('aaaaaaaa', [('a', 8)])

    def test_basic(self):
        self._roundtrip('gaaccccaaaaa', [('g', 1), ('a', 2), ('c', 4), ('a', 5)])

def suite():
    s = unittest.TestSuite()
    classes = [RLETestCase]
    for c in classes:
        s.addTests(unittest.makeSuite(c))
    return s
