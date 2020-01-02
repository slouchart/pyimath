__all__ = ['load_tests']

from importlib import import_module
from unittest import TestSuite
from unittest import main as run_tests



NB_TESTS = 10

TEST_CASES = [
    'primefield',
    'polynomial',
    'finitefield',
    'polyparse',
    'transtype',
    'pfpoly',
    'ffpoly',
    'factorize',
    'ff_factorize',
    'gaussint',
    'functions',
]


def load_tests(loader, standard_tests, pattern):
    suite = TestSuite()
    for module_stem in TEST_CASES:
        module = import_module(f'tests.test_{module_stem}')
        suite.addTest(loader.loadTestsFromModule(module))

    return suite


if __name__ == "__main__":
    run_tests(verbosity=2)
