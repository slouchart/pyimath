from unittest import TestCase
from unittest import main as run_tests

from imath.primefield import PrimeField, Polynomial


class TestPFPolynomial(TestCase):
    def setUp(self):
        self.f2 = PrimeField(2)
        self.f3 = PrimeField(3)
        self.f5 = PrimeField(5)

    def testStrAndRepr(self):
        f5 = self.f5
        p1 = f5.polynomial(-2, 1, 2, indeterminate='z')
        p2 = Polynomial([-2, 1, 2], f5, indeterminate='z')
        self.assertTrue(p1 == p2 == eval(repr(p1)))

        self.assertEqual(str(p1), '-2 + z + 2z^2')

    def testPolynomial(self):
        f5 = self.f5

        p = f5.polynomial(1, 1, 1)
        self.assertTrue(p.evaluate(f5(2)) == f5(2))

        self.assertTrue(str(f5.linear_polynomial(f5(-2))) == '2 + X')

    def testDivision(self):
        f5 = self.f5
        p = f5.polynomial(2, 1, 1)
        q, r = p.long_division(f5.polynomial(1, 1))
        self.assertTrue(str(q) == 'X' and str(r) == '2')

    def testAddConstant(self):
        f2 = self.f2
        p = f2.polynomial(f2(1), f2(1))
        p += f2(1)
        self.assertEqual(p, p.monic(1))

    def testAddConstantWithCast(self):
        f2 = self.f2
        p = f2.polynomial(1, 1)
        self.assertIsInstance(p[0], type(f2(0)))

        p += 1
        self.assertEqual(p, p.monic(1))
        self.assertIsInstance(p[1], type(f2(0)))

    def testMulConstant(self):
        f2 = self.f2
        p = f2.polynomial(f2(1), f2(1))
        p *= f2(1)
        self.assertEqual(p, f2.polynomial(f2(1), f2(1)))

        p *= f2.zero
        self.assertEqual(p, p.null)

    def testMulConstantWithCast(self):
        f2 = self.f2
        p = f2.polynomial(1, 1)
        self.assertIsInstance(p[0], type(f2(0)))

        p *= 1
        self.assertEqual(p, f2.polynomial(f2(1), f2(1)))
        self.assertIsInstance(p[1], type(f2(0)))


class TestIrreducibility(TestCase):

    def assertIsIrreducible(self, p):
        self.assertTrue(p.is_irreducible)

    def assertIsNotIrreducible(self, p):
        self.assertFalse(p.is_irreducible)

    def test1(self):
        f2 = PrimeField(2)
        p = f2.polynomial(1, 1, 1)
        self.assertIsIrreducible(p)

        p = f2.polynomial(1, 1)
        self.assertIsIrreducible(p)

        p = f2.polynomial(1, 0, 1)
        self.assertEqual(p, f2.polynomial(1, 1)**2)
        self.assertIsNotIrreducible(p)

    def test2(self):
        f3 = PrimeField(3)
        p = f3.polynomial(-1, 1, 1)
        self.assertIsIrreducible(p)

        p = f3.polynomial(1, 1, 1)
        self.assertIsNotIrreducible(p)

    def test3(self):
        f3 = PrimeField(3)
        p = f3.polynomial(-1, -1, -1, 1)
        self.assertIsIrreducible(p)

        p = f3.polynomial(1, 1, 1, 1)
        self.assertIsNotIrreducible(p)

    def test4(self):
        f3 = PrimeField(3)
        p1 = f3.polynomial(1, 0, 1)
        p2 = f3.polynomial(-1, 1, 1)
        p3 = f3.polynomial(-1, -1, 1)

        p = p1 * p2 * p3
        self.assertIsNotIrreducible(p)

        p = p ** 4
        self.assertIsNotIrreducible(p)

    def test5(self):
        f2 = PrimeField(2)
        p = f2.polynomial(1)
        p += p.monic(3)
        p += p.monic(32)
        self.assertEqual(str(p), '1 + X^3 + X^32')
        self.assertIsNotIrreducible(p)

    def test6(self):
        f3 = PrimeField(3)
        p0 = f3.polynomial(1, 0, 1, 0, 1, 0, 0, 0, 1)
        p1 = f3.polynomial(-1, 0, 0, 0, 1, 0, 0, -1, 1)
        p2 = f3.polynomial(1, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 1)
        p = p0 * p1 * p2

        self.assertIsNotIrreducible(p)


if __name__ == '__main__':
    run_tests()
