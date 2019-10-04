__all__ = ['test_all']

from unittest import TestSuite, defaultTestLoader

import tests.test_primefield as t1
import tests.test_polynomial as t2
import tests.test_finitefield as t3
import tests.test_polyparse as t4
import tests.test_transtype as t5
import tests.test_pfpoly as t6
import tests.test_ffpoly as t7
import tests.test_factorize as t8
import tests.test_ff_factorize as t9
import tests.test_gaussint as t10

NB_TESTS = 10


def test_all():
    suite = TestSuite()
    loader = defaultTestLoader
    for i in range(1, NB_TESTS+1):
        alias = f't{i}'
        module = globals()[alias]
        suite.addTest(loader.loadTestsFromModule(module))

    return suite

