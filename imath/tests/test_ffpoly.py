from unittest import TestCase
from unittest import main as run_tests

from imath import FiniteField, FFElement, PrimeField, gcd
from imath.polyparse import symbolic_polynomial


class TestStr(TestCase):

    def testF4Str(self):
        f2 = PrimeField(2)
        p_gen = f2.polynomial(1, 1, 1)
        f4 = FiniteField(2, 2, p_gen)

        p = f4.polynomial(f4(0, 1), f4(1, 0), f4(1, 1))
        self.assertTrue(str(p) == '(j) + X + (1+j)X^2')

        p *= f4(0, 1)
        self.assertTrue(str(p) == '(1+j) + (j)X + X^2')

    def testF25Repr(self):
        f5 = PrimeField(5)
        f25 = FiniteField(5, 2, f5.polynomial(2, 0, 1))

        p = symbolic_polynomial('(1+j)X^2 + (j)X^5 + X + 2X^3 + (2+j)', f25)
        self.assertEqual(p, eval(repr(f25.polynomial(f25(2, 1), 1, f25(1, 1), 2, 0, f25(0, 1)))))


class TestAdd(TestCase):

    def testF8Add(self):
        f2 = PrimeField(2)
        p = f2.polynomial(1, 1, 0, 1)
        f8 = FiniteField(2, 3, p, root_symbol='w')

        a = symbolic_polynomial('1 + X + (1+w^2)X^3', f8)
        b = symbolic_polynomial('(1+w)X^2 + X^3', f8)

        self.assertTrue(a+b == f8.polynomial(1, 1, f8(1, 1, 0), f8(0, 0, 1)))


class TestMul(TestCase):
    def testF9Product(self):
        f3 = PrimeField(3)
        p = f3.polynomial(1, 0, 1)
        f9 = FiniteField(3, 2, p)
        a = symbolic_polynomial('-1 + X - (1+j)X^2', f9)
        b = symbolic_polynomial('1 + X^2 + (j)X', f9)

        self.assertTrue(a * b == f9.polynomial(-1, f9(1, -1), 1, f9(-1, -1), f9(-1, -1)))


class TestDiv(TestCase):
    def testF16Div(self):
        f2 = PrimeField(2)
        p = f2.polynomial(1, 1, 0, 0, 1)
        f16 = FiniteField(2, 4, p)

        a = symbolic_polynomial('X^6 + X^3 + 1', f16)
        b = symbolic_polynomial('X^3 + 1', f16)

        q = a / b
        r = a % b

        self.assertEqual(q, f16.polynomial(0, 0, 0, 1))
        self.assertEqual(r, a.unit)

    def testF16gcd(self):
        f2 = PrimeField(2)
        p = f2.polynomial(1, 1, 0, 0, 1)
        f16 = FiniteField(2, 4, p)

        a = symbolic_polynomial('X^6 + X^3 + 1', f16)
        b = symbolic_polynomial('X^3 + 1', f16)

        self.assertEqual(gcd(a, b), f16.polynomial(1))


if __name__ == '__main__':
    run_tests()
