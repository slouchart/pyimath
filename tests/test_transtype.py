from unittest import TestCase
from unittest import main as run_tests

from pyimath.primefield import PrimeField
from pyimath.finitefield import FiniteField


class TestTranstyping(TestCase):

    def testPrimeFieldAndInt(self):
        """Checking correct transtyping int to PFElement for + and - in F7
        """
        f7 = PrimeField(7)
        self.assertTrue(f7(3) + 1 == -3)
        self.assertTrue(1 - f7(3) == -2)
        self.assertTrue(1 * f7(3) == 3)

        self.assertTrue(2 * f7(3) == f7(2) * 3 == -1)
        self.assertTrue(4 * f7(3) == f7(3) * 4 == -2)

    def testPrimeAndFinite(self):
        """Checking correct transtyping int, PFElement to FFElement for + and - in F9
        """
        f3 = PrimeField(3)
        f9 = FiniteField(3, 2, f3.polynomial(-1, 1, 1))

        self.assertTrue(f9(1, 0) + 1 == f3(-1))
        self.assertTrue(1 + f3(1) == f9(-1, 0))
        self.assertTrue(f9(-1, 0) == -1)

        self.assertTrue(3 * f9(-1, 0) == f9(-1, 0) * 3)
        self.assertTrue(f3(1) * f9(1, 1) == f3(1) + f9(0, 1))

    def testPrimeAndPolynomial(self):
        """Checking correct transtyping int to Polynomial for + and == in F2
        """
        f2 = PrimeField(2)
        p = f2.polynomial(1, 1, 1)
        self.assertTrue(p + 1 == f2.polynomial(0, 1, 1))

        self.assertTrue(p.monic(0) + 1 == 0)


if __name__ == '__main__':
    run_tests()
