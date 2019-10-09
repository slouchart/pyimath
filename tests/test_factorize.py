from unittest import TestCase
from unittest import main as run_tests

from pyimath.factorize import factorize
from pyimath.primefield import PrimeField
from pyimath.polynomial import symbolic_polynomial


class TestSquareFreeFactorization(TestCase):

    def test1(self):
        """Check square free factorization of an irreducible polynomial over F2
        """
        f2 = PrimeField(2)
        p1 = f2.polynomial(1, 0, 1, 1)

        sqf, factors = factorize(p1).square_free()

        self.assertTrue(len(factors) == 0)
        self.assertTrue(sqf.is_irreducible)
        self.assertTrue(sqf == p1)

    def test2(self):
        """Check square free factorization over F2
        """
        f2 = PrimeField(2)
        p1 = f2.polynomial(1, 0, 1, 0, 1)

        sqf, factors = factorize(p1).square_free()

        self.assertTrue(sqf.is_unit)
        self.assertTrue(len(factors) == 1)
        self.assertTrue(factors[0].multiplicity == 2)
        self.assertTrue(factors[0].value == f2.polynomial(1, 1, 1))

    def test3(self):
        """Check square free factorization of a quadratic polynomial over F2
        """
        f2 = PrimeField(2)
        p1 = f2.polynomial(1, 1, 1, 1, 1, 1)

        sqf, factors = factorize(p1).square_free()

        self.assertTrue(sqf == f2.polynomial(1, 1))
        self.assertTrue(len(factors) == 1)
        self.assertTrue(factors[0].value.is_irreducible)
        self.assertTrue(factors[0].multiplicity == 2)
        self.assertTrue(factors[0].value ** 2 * sqf == p1)

    def test4(self):
        """Check square free factorization of a cubic polynomial over F2
        """
        f2 = PrimeField(2)
        p0 = f2.polynomial(1, 1, 1)
        p0 **= 3
        p0 *= f2.polynomial(1, 1)
        p1 = p0.copy

        sqf, factors = factorize(p1).square_free()
        self.assertTrue(sqf == f2.polynomial(1, 1))
        self.assertTrue(sqf.is_irreducible)
        self.assertTrue(len(factors) == 1)
        self.assertTrue(factors[0].value == f2.polynomial(1, 1, 1) and factors[0].multiplicity == 3)

    def test5(self):
        """Check square free factorization over F3
        """
        f3 = PrimeField(3)
        p0 = f3.polynomial(1, -1, 0, 0, -1, 0, 0, 1, 1)
        p1 = f3.polynomial(-1, 1, -1, -1, 1)
        self.assertTrue(p0.is_irreducible and p1.is_irreducible)

        p = p0 * p1

        sqf, factors = factorize(p).square_free()
        self.assertTrue(len(factors) == 0)
        self.assertTrue(sqf == p)

    def test6(self):
        """Check square free factorization over F5
        """
        f5 = PrimeField(5)
        p = f5.polynomial(-1, 2, 1, 1)
        sqf, factors = factorize(p).square_free()
        self.assertTrue(sqf == p)

    def test7(self):
        """Check square free factorization of a quadratic polynomial over F2
        """
        f5 = PrimeField(5)
        p = f5.polynomial(1, 1)**2 * f5.polynomial(2, 1)
        sqf, factors = factorize(p).square_free()
        self.assertTrue(sqf == f5.polynomial(2, 1))
        self.assertTrue(len(factors) == 1)
        fact = factors[0]
        self.assertTrue(fact.multiplicity == 2)
        self.assertTrue(fact.value == f5.polynomial(1, 1))

    def test8(self):
        """Check square free factorization of a cubic polynomial over F3
        """
        f3 = PrimeField(3)
        p = f3.polynomial(1, 1)**3 * f3.polynomial(-1, 1)
        sqf, factors = factorize(p).square_free()
        self.assertTrue(sqf == f3.polynomial(-1, 1))
        self.assertTrue(len(factors) == 1)
        fact = factors[0]
        self.assertTrue(fact.multiplicity == 3)
        self.assertTrue(fact.value == f3.polynomial(1, 1))

    def test9(self):
        """Check square free factorization of a polynomial of degree 9 over F3
        """
        f3 = PrimeField(3)
        p = f3.polynomial(-1, 1)**3 * f3.polynomial(1, 1)**2 * f3.polynomial(-1, 1, 1)
        sqf, factors = factorize(p).square_free()
        self.assertTrue(sqf == f3.polynomial(-1, 1, 1))
        self.assertTrue(len(factors) == 2)
        for fct in factors:
            self.assertTrue(fct.multiplicity == 2 and fct.value == f3.polynomial(1, 1)
                            or fct.multiplicity == 3 and fct.value == f3.polynomial(-1, 1))


class TestDistinctDegree(TestCase):
    def test1(self):
        """Check distinct degree factorization over F2
        """
        f2 = PrimeField(2)
        p0 = f2.polynomial(1, 1)
        p1 = f2.polynomial(1, 1, 1)
        p = p0 * p1

        factors = factorize(p).distinct_degree()
        self.assertTrue(all([f.is_irreducible for f in factors]))
        self.assertTrue(len(factors) == 2)
        self.assertTrue(p == factors[0].value * factors[1].value)

    def test2(self):
        """Check distinct degree factorization of a non square free polynomial over F2
        """
        f2 = PrimeField(2)
        p = f2.polynomial(1, 0, 1)  # non square free
        factors = factorize(p).distinct_degree()
        self.assertTrue(all([f.is_irreducible for f in factors]))
        self.assertTrue(len(factors) == 2)
        self.assertTrue(p == factors[0].value * factors[1].value)
        self.assertTrue(factors[0].value == factors[1].value)

    def test3(self):
        """Check distinct degree factorization of a square free equal degree polynomial over F3
        """
        f3 = PrimeField(3)
        p1 = f3.polynomial(1, 1)
        p2 = f3.polynomial(-1, 1)
        p = p1 * p2  # square free equal degree
        factors = factorize(p).distinct_degree()
        self.assertTrue(len(factors) == 1)
        self.assertTrue(factors[0].value == p)

    def test4(self):
        """Check distinct degree factorization over F3
        """
        f3 = PrimeField(3)
        p1 = f3.polynomial(1, 1)
        p2 = f3.polynomial(-1, 1, 1)

        p = p1 * p2
        factors = factorize(p).distinct_degree()
        self.assertTrue(len(factors) == 2)
        self.assertTrue(all([f.is_irreducible for f in factors]))
        self.assertTrue(factors[0].value * factors[1].value == p)
        self.assertTrue(factors[0].value.degree != factors[1].value.degree)

    def test5(self):
        """Check distinct degree factorization of a equal degree polynomial over F5
        """
        f5 = PrimeField(5)
        p1 = f5.polynomial(2, 0, 1)
        p2 = f5.polynomial(1, 1, 1)

        p = p1 * p2
        factors = factorize(p).distinct_degree()
        self.assertTrue(len(factors) == 1)
        self.assertTrue(factors[0].value == p)

    def test6(self):
        """Check distinct degree factorization over F5
        """
        f5 = PrimeField(5)
        p1 = f5.polynomial(2, 0, 1)
        p2 = f5.polynomial(1, 1, 0, 1)

        p = p1 * p2
        factors = factorize(p).distinct_degree()
        self.assertTrue(len(factors) == 2)
        self.assertTrue(factors[0].max_degree == factors[0].value.degree)
        self.assertTrue(factors[1].max_degree == factors[1].value.degree)

    def test7(self):
        """Check distinct degree factorization over F5 with:
        - two irreducible factors of the same degree
        - one irreducible factor of distinct degree
        """
        f5 = PrimeField(5)
        p1 = f5.polynomial(2, 0, 1)
        p2 = f5.polynomial(1, -1, 1)
        p3 = f5.polynomial(1, 1, 0, 1)
        p = p1 * p2 * p3

        factors = factorize(p).distinct_degree()
        self.assertTrue(len(factors) == 2)

        fdeg4 = [f for f in factors if f.max_degree == 2][0]
        fdeg3 = [f for f in factors if f.max_degree == 3][0]

        self.assertTrue(fdeg4.value.degree == 2*fdeg4.max_degree == 4)
        self.assertTrue(fdeg3.value.degree == fdeg3.max_degree == 3)

        self.assertTrue(fdeg4.value == p1 * p2)
        self.assertTrue(fdeg3.value == p3)


class TestEqualDegree(TestCase):
    def test1(self):
        """Check equal degree factorization over F5
        (non-deterministic, might be skipped)
        """
        f5 = PrimeField(5)
        p = f5.polynomial(1, 1) * f5.polynomial(2, 1)

        try:
            factors = factorize(p).equal_degree(nb_factors=2, max_degree=1)
            self.assertTrue(len(factors) == 2)
            self.assertTrue(all([f.value.degree == 1 for f in factors]))
            self.assertTrue(factors[0].value * factors[1].value == p)
        except RuntimeError as e:
            self.skipTest(f'Non-deterministic factorization failed: {e}')

    def test2(self):
        """Check equal degree factorization over F2
        (non-deterministic, might be skipped)
        """
        f2 = PrimeField(2)
        p = f2.polynomial(0, 1, 1)

        try:
            factors = factorize(p).equal_degree(nb_factors=2, max_degree=1)
            self.assertTrue(len(factors) == 2)
            self.assertTrue(all([f.value.degree == 1 for f in factors]))
            self.assertTrue(factors[0].value * factors[1].value == p)
        except RuntimeError as e:
            self.skipTest(f'Non-deterministic factorization failed: {e}')

    def test3(self):
        """Check equal degree factorization over F3
        (non-deterministic, might be skipped)
        """
        f3 = PrimeField(3)
        p1 = f3.polynomial(-1, 1, 1)
        p2 = f3.polynomial(1, 0, 1)

        p = p1 * p2
        try:
            factors = factorize(p).equal_degree(nb_factors=2, max_degree=2)
            self.assertTrue(len(factors) == 2)
            self.assertTrue(all([f.value.degree == 2 for f in factors]))
            self.assertTrue(factors[0].value * factors[1].value == p)
        except RuntimeError as e:
            self.skipTest(f'Non-deterministic factorization failed: {e}')

    def test4(self):
        """Check equal degree factorization over F3: parameter mismatch
        """
        f3 = PrimeField(3)
        p1 = f3.polynomial(-1, 1, 1)
        p2 = f3.polynomial(1, 0, 1)

        p = p1 * p2
        with self.assertRaises(AssertionError) as a_exc:
            _ = factorize(p).equal_degree(nb_factors=3, max_degree=2)

        self.assertEqual(str(a_exc.exception), f'Degree of {str(p)} mismatches input parameters 3, 2')

    def test5(self):
        """Check equal degree factorization over F2: parameter mismatch
        """
        f2 = PrimeField(2)
        p = f2.polynomial(1, 1, 1)

        # bad parameters
        with self.assertRaises(RuntimeError) as a_exc:
            _ = factorize(p).equal_degree(2, 1)

        self.assertEqual(str(a_exc.exception), 'unable to find 2 degree 1 factors for 1 + X + X^2')

    def test6(self):
        """Check equal degree factorization over F2: non square free polynomial
        """
        f2 = PrimeField(2)
        p = f2.polynomial(1, 1)
        p **= 3

        # non square-free polynomial
        with self.assertRaises(RuntimeError) as a_exc:
            _ = factorize(p).equal_degree(3, 1)

        self.assertEqual(str(a_exc.exception), 'unable to find 3 degree 1 factors for 1 + X + X^2 + X^3')

    def test7(self):
        """Check equal degree factorization over F3: non equal degree
        Raises most of the time but due to the non-deterministic nature of the Cantor-Zassenhaus algorithm,
        it might find a factor anyway
        """
        f3 = PrimeField(3)
        p1 = f3.polynomial(1, 1)
        p2 = f3.polynomial(1, -1, 0, 1)

        p = p1 * p2

        try:
            factors = factorize(p).equal_degree(2, 2)
            self.assertTrue(len(factors) == 2)
            self.assertTrue(factors[0].value.degree + factors[1].value.degree == p.degree)
        except RuntimeError as exc:
            self.assertEqual(str(exc), 'unable to find 2 degree 2 factors for 1 - X^2 + X^3 + X^4')


class TestFullFactorization(TestCase):
    def test1(self):
        """Check full factorization over F3
        """
        f3 = PrimeField(3)

        c0 = f3.polynomial(1, 1)  # irred deg 1
        c1 = f3.polynomial(0, 0, 1)  # trivial square
        c2 = f3.polynomial(-1, 1) ** 3
        c3 = f3.polynomial(-1, 1, 1)

        p = c0 * c1 * c2 * c3

        factors, c = factorize(p).cantor_zassenhaus()
        self.assertTrue(c == 1)
        self.assertTrue(len(factors) == 4)
        self.assertEqual(factorize(p).factors_product(factors), p)

    def test2(self):
        """Check full factorization over F5
        """
        f5 = PrimeField(5)
        p = f5.polynomial(1, 2, 0, 2, -1, 1)

        try:
            factors, c = factorize(p).cantor_zassenhaus()
            self.assertTrue(c == f5(1))
            self.assertTrue(factorize(p).factors_product(factors) == p)
        except RuntimeError as e:
            self.skipTest(f'Non-deterministic factorization failed: {e}')

    def test3(self):
        """Check full factorization over F2, square of degree 8
        """
        f2 = PrimeField(2)
        p = symbolic_polynomial('1 + X^2 + X^4 + X^8', f2)

        factors, c = factorize(p).cantor_zassenhaus()
        self.assertTrue(c == 1)
        self.assertTrue(len(factors) == 2)
        self.assertTrue(factorize(p).factors_product(factors) == p)

    def test4(self):
        """Check full factorization over F3, non monic
        """
        f3 = PrimeField(3)
        p = symbolic_polynomial('-X^3 - X^2 + X + 1', f3)

        factors, c = factorize(p).cantor_zassenhaus()
        self.assertTrue(factorize(p).factors_product(factors) * c == p)

    def test5(self):
        """Check full factorization over F7, non monic
        """
        f7 = PrimeField(7)
        p = symbolic_polynomial('-3X^3 - 2X^2 + X + 3', f7)

        factors, c = factorize(p).cantor_zassenhaus()
        self.assertTrue(factorize(p).factors_product(factors) * c == p)


if __name__ == '__main__':
    run_tests()
