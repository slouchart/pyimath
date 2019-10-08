from unittest import TestCase
from unittest import skip
from unittest import main as run_tests

from pyimath.finitefield import FiniteField, FFElement, finite_field
from pyimath.polynomial import Polynomial
from pyimath.primefield import PrimeField, PFElement


class TestF4(TestCase):

    def setUp(self):
        self.f4 = FiniteField(2, 2, PrimeField(2).polynomial(1, 1, 1))

    def testStr(self):
        """Check printable representation of a finite field
        """
        f4 = self.f4
        self.assertEqual(str(f4), 'Finite field of order 2^2')

    @skip
    def testRepr(self):
        """Check evaluable representation of a finite field (skipped)
        """
        f4 = self.f4
        s0 = repr(f4)
        self.assertEqual(f4, eval(s0))

    def testSum(self):
        """Check addition in a finite field
        """
        f4 = self.f4
        self.assertTrue(f4(0, 1) + f4(1, 0) == f4(1, 1))

    def testDiff(self):
        """Check additive inverse and subtraction in a finite field
        """
        f4 = self.f4
        self.assertTrue(f4(1, 1) - f4(1, 0) == f4(0, 1))

    def testProd(self):
        """Check multiplication in a finite field
        """
        f4 = self.f4
        p = f4(1, 1) * f4(1, 1)
        self.assertTrue(p == f4(0, 1))


class TestF9(TestCase):
    def testPow(self):
        """Check exponentiation in F9
        """
        f9 = FiniteField(3, 2, PrimeField(3).polynomial(1, 0, 1))
        self.assertTrue(f9(1, 1) ** 2 == f9(0, -1))


class TestF8(TestCase):
    def setUp(self):
        self.f8 = FiniteField(2, 3, PrimeField(2).polynomial(1, 1, 0, 1))

    def testPower(self):
        """Check exponentiation in F8
        """
        f8 = self.f8
        self.assertTrue(f8(1, 1, 1) ** 2 == f8(1, 1, 0))

    def testDiv(self):
        """Check division in F9
        """
        f8 = self.f8
        n = f8(1, 1, 1)
        d = f8(1, 1, 0)

        self.assertTrue(n / d == f8(0, 0, 1))

    def testStr(self):
        """Check printable representation of an element of in F9
        """
        f8 = self.f8
        e = f8(1, 1, 1)
        self.assertEqual(str(e), '1 + j + j^2')


class TestF25(TestCase):
    def setUp(self):
        self.f25 = FiniteField(5, 2, Polynomial([1, 1, 1], base_field=PrimeField(5)))

    def testMul(self):
        """Check multiplication in F25
        """
        f25 = self.f25
        a = f25(1, 2)
        b = f25(1, 1)
        self.assertEqual(a * b, f25(-1, 1))

    def testDiv(self):
        """Check true division in F25
        """
        f25 = self.f25
        n = f25(-1, 1)
        d = f25(1, 1)
        self.assertEqual(n / d, f25(1, 2))

    def testDivComplete(self):
        """Check equivalence of true division and multiplication with the multiplicative inverse in F25
        """
        f25 = self.f25
        for e in f25:
            if e == f25.zero:
                continue
            inv_e = f25.one / e
            self.assertEqual(e * inv_e, f25.one)

    def testDivByZero(self):
        """Sanity check: Division by zero"""
        f25 = self.f25
        with self.assertRaises(ZeroDivisionError):
            f25.div(f25.one, f25.zero)

        with self.assertRaises(ZeroDivisionError):
            f25.multiplicative_inverse(f25.zero)


class TestF27(TestCase):
    def setUp(self):
        self.f27 = FiniteField(3, 3, Polynomial([-1, -1, 0, 1], base_field=PrimeField(3)))

    def testStr(self):
        """Check printable representation of element of F27
        """
        f27 = self.f27
        self.assertEqual(str(f27(1, -1, -1)), '1 - j - j^2')


class TestTypingSanity(TestCase):
    def test1(self):
        """Check return type of power in F4 and F2
        """
        f2 = PrimeField(2)
        f4 = FiniteField(2, 2, f2.polynomial(1, 1, 1))

        b = f4(0, 1) ** 0
        self.assertIsInstance(b, FFElement)
        self.assertTrue(b == 1)

        b = f2(1) ** 0
        self.assertIsInstance(b, PFElement)
        self.assertTrue(b == 1)


class TestGenerator(TestCase):

    def setUp(self) -> None:

        self.f2 = finite_field(2)
        self.f4 = FiniteField(2, 2, self.f2.polynomial(1, 1, 1))
        self.f4g = FiniteField(2, 2, self.f2.polynomial(1, 1, 1), generator=(0, 1))

    def test0(self):
        """Check generator order
        """
        self.assertEqual(self.f4g.generator ** (self.f4g.order-1), self.f4g.one)

    def test1(self):
        """Check equivalence between multiplication with and w/o use of generator
        """
        for a in self.f4:
            for b in self.f4:
                self.assertEqual(a * b, self.f4g.element(a) * self.f4g.element(b))

    def test2(self):
        """Check equivalence between true division with and w/o use of generator
        """
        for a in self.f4:
            for b in self.f4:
                if not b.null:
                    self.assertEqual(a / b, self.f4g.element(a) / self.f4g.element(b))

    def test3(self):
        """Check equivalence between multiplication with and w/o use of generator in F16
        """
        f16g = finite_field(16)
        f2 = finite_field(2)
        f16 = FiniteField(2, 4, f2.polynomial(1, 1, 0, 0, 1))
        for a in f16:
            for b in f16:
                self.assertEqual(a * b, f16g.element(a) * f16g.element(b))
                if not b.null:
                    self.assertEqual(a / b, f16g.element(a) / f16g.element(b))

    def test4(self):
        """Check equivalence between multiplication with and w/o use of generator in F25
        """
        f25g = finite_field(25)
        f5 = finite_field(5)
        f25 = FiniteField(5, 2, f5.polynomial(2, 0, 1))
        for a in f25:
            for b in f25:
                self.assertEqual(a * b, f25g.element(a) * f25g.element(b))
                if not b.null:
                    self.assertEqual(a / b, f25g.element(a) / f25g.element(b))

    def test5(self):
        """Check equivalence between multiplication with and w/o use of generator in F27
        """
        f27g = finite_field(27)
        f3 = finite_field(3)
        f27 = FiniteField(3, 3, f3.polynomial(-1, -1, 0, 1))
        for a in f27:
            for b in f27:
                self.assertEqual(a * b, f27g.element(a) * f27g.element(b))
                if not b.null:
                    self.assertEqual(a / b, f27g.element(a) / f27g.element(b))


if __name__ == '__main__':
    run_tests()
