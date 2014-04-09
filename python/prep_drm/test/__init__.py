import unittest

modules = [
        'prep_drm.test.test_pwalign',
        'prep_drm.test.test_rle',
        'prep_drm.test.test_samutil',
        'prep_drm.test.test_scripts_classify_mutations',
        'prep_drm.test.test_scripts_hypermutation',
        'prep_drm.test.test_ssw',
        'prep_drm.test.test_util',
]

def suite():
    suite = unittest.TestSuite()
    for m in modules:
        mod = __import__(m, fromlist=[m.rsplit('.', 1)[-1]])
        suite.addTest(mod.suite())

    return suite

suite = suite()
