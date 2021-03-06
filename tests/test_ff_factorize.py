from unittest import TestCase
from unittest import main as run_tests


from pyimath.finitefield import finite_field
from pyimath.factorize import factorize


class TestFiniteFieldFactorization(TestCase):

    def testFiniteFieldF4(self):
        """Check full Cantor-Zassenhaus factorization over F4
        """
        f4 = finite_field(4)

        p1 = f4.polynomial(1, f4(1, 1), 1)
        p2 = f4.polynomial(1, f4(1, 1), 0, 1)
        p3 = f4.polynomial(f4(0, 1), 1) ** 2
        p4 = f4.polynomial(f4(0, 1), 1, f4(0, 1), 0, f4(1, 1)) / f4(1, 1)

        ps = [p1, p2, p3, p4]

        for p in ps:
            with self.subTest(f'Testing factorizing {str(p)}'):
                factors, c = factorize(p).cantor_zassenhaus()
                self.assertTrue(p == factorize(p).factors_product(factors) * c)

    def testFiniteFieldF9(self):
        """Check full Cantor-Zassenhaus factorization over F9
        """
        f9 = finite_field(9)

        p1 = f9.polynomial(1, f9(1, 1), 1)
        p2 = f9.polynomial(1, f9(1, 1), 0, 1)
        p3 = f9.polynomial(f9(0, 1), 1) ** 2
        p4 = f9.polynomial(f9(0, 1), 1, f9(0, 1), 0, f9(1, 1)) / f9(1, 1)
        p5 = f9.polynomial(f9(-1, 1), f9(1, -1), 1) ** 3

        ps = [p1, p2, p3, p4, p5]

        for p in ps:
            with self.subTest(f'Testing factorizing {str(p)}'):
                factors, c = factorize(p).cantor_zassenhaus()
                self.assertTrue(p == factorize(p).factors_product(factors) * c)

    def testFiniteFieldF8(self):
        """Check full Cantor-Zassenhaus factorization over F8
        """
        f8 = finite_field(8)

        p1 = f8.polynomial(1, f8(1, 1), 1)
        p2 = f8.polynomial(1, f8(1, 1), 0, 1)
        p3 = f8.polynomial(f8(0, 1), 1) ** 2
        p4 = f8.polynomial(f8(0, 1), 1, f8(0, 1), 0, f8(1, 1)) / f8(1, 1)
        p5 = f8.polynomial(f8(1, 1), f8(1, 1, 1), 1) ** 3

        ps = [p1, p2, p3, p4, p5]

        for p in ps:
            with self.subTest(f'Testing factorizing {str(p)}'):
                factors, c = factorize(p).cantor_zassenhaus()
                self.assertTrue(p == factorize(p).factors_product(factors) * c)


if __name__ == '__main__':
    run_tests()
