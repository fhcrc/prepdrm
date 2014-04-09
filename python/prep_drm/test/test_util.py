import unittest

from .. import util


class IsolateRegionTestCase(unittest.TestCase):
    def setUp(self):
        self.ap = [(1L, 199L),
                   (2L, 200L),
                   (3L, None),
                   (4L, 201L),
                   (5L, 202L),
                   (6L, 203L),
                   (7L, 204L),
                   (None, 205L),
                   (None, 206L),
                   (None, 207L),
                   (None, 208L),
                   (8L, 209L),
                   (9L, 210L)]

    def test_before_region(self):
        self.assertEqual([], util.isolate_region(self.ap, 50, 75))

    def test_after_region(self):
        self.assertEqual([], util.isolate_region(self.ap, 215, 224))

    def test_whole_region(self):
        self.assertEqual(self.ap, util.isolate_region(self.ap, 199, 211))
        self.assertEqual(self.ap, util.isolate_region(self.ap, 0, 1000))

    def test_partial_rdel(self):
        self.assertEqual(self.ap[1:6], util.isolate_region(self.ap, 200, 204))

    def test_partial_qdel(self):
        self.assertEqual(self.ap[6:11], util.isolate_region(self.ap, 204, 209))
        self.assertEqual(self.ap[9:11], util.isolate_region(self.ap, 207, 209))


def suite():
    suite = unittest.TestSuite()
    for c in [IsolateRegionTestCase]:
        suite.addTest(unittest.makeSuite(c))
    return suite
