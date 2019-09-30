from unittest import TestCase
from unittest import skip
from unittest import main as run_tests

from imath.finitefield import FiniteField, FFElement
from imath.polynomial import Polynomial
from imath.primefield import PrimeField, PFElement


class TestF4(TestCase):

    def setUp(self):
        self.f4 = FiniteField(2, 2, PrimeField(2).polynomial(1, 1, 1))

    def testStr(self):
        f4 = self.f4
        self.assertEqual(str(f4), 'Finite field of order 2^2')

    @skip
    def testRepr(self):
        f4 = self.f4
        s0 = repr(f4)
        self.assertEqual(f4, eval(s0))

    def testSum(self):
        f4 = self.f4
        self.assertTrue(f4(0, 1) + f4(1, 0) == f4(1, 1))

    def testDiff(self):
        f4 = self.f4
        self.assertTrue(f4(1, 1) - f4(1, 0) == f4(0, 1))

    def testProd(self):
        f4 = self.f4
        p = f4(1, 1) * f4(1, 1)
        self.assertTrue(p == f4(0, 1))


class TestF9(TestCase):
    def testPow(self):
        f9 = FiniteField(3, 2, PrimeField(3).polynomial(1, 0, 1))
        self.assertTrue(f9(1, 1) ** 2 == f9(0, -1))


class TestF8(TestCase):
    def setUp(self):
        self.f8 = FiniteField(2, 3, PrimeField(2).polynomial(1, 1, 0, 1))

    def testPower(self):
        f8 = self.f8
        self.assertTrue(f8(1, 1, 1) ** 2 == f8(1, 1, 0))

    def testDiv(self):
        f8 = self.f8
        n = f8(1, 1, 1)
        d = f8(1, 1, 0)

        self.assertTrue(n / d == f8(0, 0, 1))

    def testStr(self):
        f8 = self.f8
        e = f8(1, 1, 1)
        self.assertEqual(str(e), '1 + j + j^2')


class TestF25(TestCase):
    def setUp(self):
        self.f25 = FiniteField(5, 2, Polynomial([1, 1, 1], base_field=PrimeField(5)))

    def testMul(self):
        f25 = self.f25
        a = f25(1, 2)
        b = f25(1, 1)
        self.assertEqual(a * b, f25(-1, 1))

    def testDiv(self):
        f25 = self.f25
        n = f25(-1, 1)
        d = f25(1, 1)
        self.assertEqual(n / d, f25(1, 2))

    def testDivComplete(self):
        f25 = self.f25
        for e in f25:
            if e == f25.zero:
                continue
            inv_e = f25.one / e
            self.assertEqual(e * inv_e, f25.one)

    def testDivByZero(self):
        f25 = self.f25
        with self.assertRaises(ZeroDivisionError):
            f25.div(f25.one, f25.zero)

        with self.assertRaises(ZeroDivisionError):
            f25.multiplicative_inverse(f25.zero)


class TestF27(TestCase):
    def setUp(self):
        self.f27 = FiniteField(3, 3, Polynomial([-1, -1, 0, 1], base_field=PrimeField(3)))

    def testStr(self):
        f27 = self.f27
        self.assertEqual(str(f27(1, -1, -1)), '1 - j - j^2')


class TestTypingSanity(TestCase):
    def test1(self):
        f2 = PrimeField(2)
        f4 = FiniteField(2, 2, f2.polynomial(1, 1, 1))

        b = f4(0, 1) ** 0
        self.assertIsInstance(b, FFElement)
        self.assertTrue(b == 1)

        b = f2(1) ** 0
        self.assertIsInstance(b, PFElement)
        self.assertTrue(b == 1)


class TestGenerator(TestCase):
    pass


if __name__ == '__main__':
    run_tests()
