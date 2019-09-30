from unittest import TestCase
from unittest import main as run_tests


from imath.polyparse import symbolic_polynomial
from imath.finitefield import PrimeField, FiniteField
from imath.integer import IntegerRing
from imath.polynomial import Polynomial


class TestParserBasic(TestCase):

    def test1(self):
        expr = '1 + X + X^2 + X^3'
        p = symbolic_polynomial(expr, PrimeField(2))
        self.assertEqual(str(p), expr)

    def test2(self):
        expr = 'X^3 + 2X - 1'
        p = symbolic_polynomial(expr, PrimeField(5))
        self.assertEqual(p, p.base_field.polynomial(-1, 2, 0, 1))


class TestParserAdvanced(TestCase):
    def setUp(self):
        f2 = PrimeField(2)
        p_gen = f2.polynomial(1, 1, 1)
        self.f4 = FiniteField(2, 2, p_gen)

    def test1(self):
        p = Polynomial([0, 1, 1], IntegerRing())
        self.assertTrue(str(p) == 'X + X^2')

        f5 = PrimeField(5)
        p = f5.linear_polynomial(f5(-2))
        self.assertTrue(str(p) == '2 + X')

        p = f5.polynomial(0, -1, 1)
        self.assertTrue(str(p) == '-X + X^2')

        p = symbolic_polynomial('-2 - X^16', PrimeField(5))
        self.assertTrue(str(p) == '-2 - X^16')

    def test2(self):
        f4 = self.f4

        p = symbolic_polynomial('1', f4)
        self.assertTrue(str(p) == '1')

        p = symbolic_polynomial('X+1', f4)
        self.assertTrue(str(p) == '1 + X')

        p = symbolic_polynomial('(j)', f4)
        self.assertTrue(str(p) == '(j)')

        p = symbolic_polynomial('(1+j) + X^2', f4)
        self.assertTrue(str(p) == '(1+j) + X^2')

        p = symbolic_polynomial('(1+j)X^2', f4)
        self.assertTrue(str(p) == '(1+j)X^2')

        p = symbolic_polynomial('1 + (j)X', f4)
        self.assertTrue(str(p) == '1 + (j)X')

        p = symbolic_polynomial('1 + (1+j)X + (j)X^2', f4)
        self.assertTrue(str(p) == '1 + (1+j)X + (j)X^2')


class TestParserComplete(TestCase):

    def test1(self):
        f2 = PrimeField(2)
        f4 = FiniteField(2, 2, f2.polynomial(1, 1, 1))

        expr = '+(j+1)X + (1+j)X^2 + (j)'
        p = symbolic_polynomial(expr, f4)

        self.assertTrue(p == f4.polynomial(f4(0, 1), f4(1, 1), f4(1, 1)))

    def test2(self):
        expr = '-1'
        p = symbolic_polynomial(expr, IntegerRing())
        self.assertTrue(p == IntegerRing().polynomial(-1))

    def test3(self):
        expr = '-1 - X'
        p = symbolic_polynomial(expr, PrimeField(3))
        self.assertTrue(p == PrimeField(3).polynomial(-1, -1))

    def test4(self):
        expr = '-(1+w+w^2) + X'
        f2 = PrimeField(2)
        g = f2.polynomial(1, 1, 0, 1)
        f8 = FiniteField(2, 3, g, root_symbol='w')

        p = symbolic_polynomial(expr, f8)
        self.assertTrue(p == f8.polynomial(f8(1, 1, 1), 1))

    def test5(self):
        expr = '1 + X - X^2 + X^6'
        p = symbolic_polynomial(expr, IntegerRing())
        self.assertTrue(p == IntegerRing().polynomial(1, 1, -1, 0, 0, 0, 1))

    def test6(self):
        expr = '-'
        with self.assertRaises(SyntaxError):
            symbolic_polynomial(expr, IntegerRing())

    def test7(self):
        expr = '++X'
        with self.assertRaises(SyntaxError):
            symbolic_polynomial(expr, IntegerRing())

    def test8(self):
        expr = '(j+)X'

        f2 = PrimeField(2)
        p = f2.polynomial(1, 1, 1)
        f4 = FiniteField(2, 2, p)

        with self.assertRaises(SyntaxError):
            symbolic_polynomial(expr, f4)

    def test9(self):
        expr = 'X^'
        with self.assertRaises(SyntaxError):
            symbolic_polynomial(expr, IntegerRing())

    def test10(self):
        expr = 'X^ + 1'
        with self.assertRaises(SyntaxError):
            symbolic_polynomial(expr, IntegerRing())

    def test11(self):
        expr = '-X + X^2 X^3'
        with self.assertRaises(SyntaxError):
            symbolic_polynomial(expr, IntegerRing())

    def test12(self):
        expr = ''
        p = symbolic_polynomial(expr, IntegerRing())
        self.assertTrue(p == IntegerRing().polynomial(IntegerRing().zero))


if __name__ == '__main__':
    run_tests()
