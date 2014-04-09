import unittest

from ..pwalign import pw_global

class PairwiseGlobalTestCase(unittest.TestCase):

    def test_basic(self):
        res = pw_global('TACAGTTTA', 'ACCGT')[:3]
        expected = 'TACAGTTTA', '-ACCGT---', 2
        self.assertEqual(expected, res)


def suite():
    s = unittest.TestSuite()
    s.addTests(unittest.makeSuite(PairwiseGlobalTestCase))
    return s
