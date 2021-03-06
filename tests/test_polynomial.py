from unittest import TestCase
from unittest import main as run_tests

from pyimath.integer import IntegerRing


def polynomial(*coeffs, indeterminate='X'):
    return IntegerRing().polynomial(*coeffs, indeterminate=indeterminate)


class TestPolynomialProperties(TestCase):
    def setUp(self):
        self.p = polynomial(1, 1)

    def test_monic(self):
        """Check polynomial property: is_monic
        """
        self.assertTrue(self.p.is_monic)

    def test_degree_one(self):
        """Check polynomial property: degree
        """
        self.assertTrue(self.p.degree == 1)

    def test_valuation_zero(self):
        """Check polynomial property: valuation
        """
        self.assertTrue(self.p.valuation == 0)

    def test_identity(self):
        """Check polynomial property: operator ==
        """
        self.assertEqual(self.p, self.p.copy)


class TestPrototyping(TestCase):
    def testProto(self):
        """Check prototyping from base field
        """
        proto = polynomial(1)
        p = proto(1, 2, 1)
        self.assertEqual(p.coefficients, [1, 2, 1])


class TestNullPolynomial(TestCase):
    def setUp(self):
        self.p = polynomial(0)

    def testDegreeZero(self):
        """Check degree 0"""
        self.assertEqual(self.p.degree, 0)

    def testUndefinedValuation(self):
        """Check undefined valuation for degree 0
        """
        def _caller():
            return self.p.valuation
        with self.assertRaises(ValueError):
            _caller()

    def testEqualNull(self):
        """Check equals null for degree 0
        """
        self.assertEqual(self.p, self.p.null)


class TestMonicPolynomial(TestCase):
    def setUp(self):
        self.p = polynomial(0, 1)

    def testDegreeOne(self):
        """Check monic polynomial of degree 1
        """
        self.assertEqual(self.p.degree, 1)

    def testValuationNonZero(self):
        """Check non zero valuation of monic polynomial
        """
        self.assertEqual(self.p.valuation, 1)

    def testMonicPoly(self):
        """Check making a monic polynomial w/ monic method
        """
        self.assertTrue(hasattr(self.p, 'monic'))
        self.assertEqual(self.p.monic(3).degree, 3)
        self.assertEqual(self.p.monic(1).valuation, 1)


class TestValuationAndDegree(TestCase):
    def setUp(self):
        self.pdeg2 = polynomial(0, 0, 1)
        self.pdeg3 = polynomial(0, 0, 2, 1)

    def testDegree(self):
        """Extended check: degree
        """
        self.assertEqual(self.pdeg2.degree, 2)
        self.assertEqual(self.pdeg3.degree, 3)

    def testValuation(self):
        """Extended check: valuation
        """
        self.assertEqual(self.pdeg2.valuation, 2)
        self.assertEqual(self.pdeg3.valuation, 2)


class TestConstantTerm(TestCase):
    def setUp(self):
        self.p = polynomial(1)

    def testConstant(self):
        """Check `constant` property
        """
        self.assertEqual(self.p.constant, 1)
        self.assertEqual(self.p.null.constant, 0)
        self.assertEqual(self.p.monic(3).constant, 0)


class TestLeadAndTrailCoefficients(TestCase):
    def setUp(self):
        self.p = polynomial(1, 2, 3)

    def testLeadingCoefficient(self):
        """Check `leading` property
        """
        self.assertTrue(hasattr(self.p, 'leading'))
        self.assertTrue(self.p.leading == 3)

    def testTrailingCoefficient(self):
        """Check `trailing` property
        """
        self.assertTrue(hasattr(self.p, 'trailing'))
        self.assertTrue(self.p.trailing == 1)

    def testEquality(self):
        """Check operator == for degree > 1
        """
        p = polynomial(0, 1, 1)
        self.assertEqual(p, polynomial(0, 1, 1))


class TestSumAndDiff(TestCase):
    def testSumSameDegree(self):
        """Check sum of polynomials of equal degree
        """
        self.assertEqual(polynomial(1, 2) + polynomial(0, 2), polynomial(1, 4))

    def testSumMonic(self):
        """Check sum of monic polynomials
        """
        p = polynomial(1)
        self.assertEqual(p + p.monic(1) + p.monic(2), polynomial(1, 1, 1))

    def testSum(self):
        """Check sum of polynomials of distinct degree
        """
        self.assertEqual(polynomial(1, 2) + polynomial(2, 3, 4), polynomial(3, 5, 4))

    def testDiffNull(self):
        """Subtraction of polynomial: sanity check
        """
        p = polynomial(1, 1)
        self.assertEqual(p - p, p.null)
        self.assertEqual(p - p.null, p)

    def testDiffSameDegree(self):
        """Check subtraction of polynomials of equal degree
        """
        self.assertEqual(polynomial(1, 1, 1) - polynomial(0, 0, 1), polynomial(1, 1))

    def testDiff(self):
        """Check subtraction of polynomials of distinct degree
        """
        p = polynomial(1, 2, 1)
        p -= polynomial(2, 1, 0)
        self.assertEqual(p, polynomial(-1, 1, 1))

    def testAddConstant(self):
        """Check adding a constant term
        """
        p = polynomial(1, 2, 1)
        self.assertTrue(1 + p == p + 1 == polynomial(2, 2, 1))


class TestProductAndPow(TestCase):

    def setUp(self):
        self.p = polynomial(1, 1)

    def testUnitProduct(self):
        """Check unitary product
        """
        unit = polynomial(1)
        self.assertEqual(self.p*unit, self.p)

    def testNullProduct(self):
        """Check null product
        """
        self.assertEqual(self.p * self.p.null, self.p.null)

    def testSquare(self):
        """Check squaring a monic polynomial of degree 1
        """
        self.assertEqual(self.p * self.p, polynomial(1, 2, 1))

    def testCubePow(self):
        """Check cubing a monic polynomial of degree 1
        """
        self.assertEqual(self.p**3, polynomial(1, 3, 3, 1))

    def testProductPowerCompatibility(self):
        """Check product <-> power compatibility
        """
        p2 = self.p**2
        p3 = self.p**3
        self.assertEqual(p2*p3, polynomial(1, 5, 10, 10, 5, 1))

    def testProduct(self):
        """Extended check: product
        """
        n = polynomial(1, 3, 3, 1)
        d = polynomial(1, 1)
        t = polynomial(0, 0, 1)
        self.assertTrue(d * t == t * d == polynomial(0, 0, 1, 1))
        self.assertEqual(n - d * t, polynomial(1, 3, 2))

    def testMulConstant(self):
        """Check multiplying by a constant
        """
        p = polynomial(1, 2, 3)
        self.assertTrue(p * 3 == 3 * p == polynomial(3, 6, 9))


class TestLongDivision(TestCase):

    def testDividePowers(self):
        """Check long division cube/factor and square/factor
        """
        self.assertEqual(polynomial(1, 3, 3, 1) / polynomial(1, 1), polynomial(1, 2, 1))
        self.assertEqual(polynomial(1, 2, 1) / polynomial(1, 1), polynomial(1, 1))

        self.assertTrue(polynomial(1, 3, 3, 1) % polynomial(1, 1) == polynomial(0))
        self.assertTrue(polynomial(1, 2, 1) % polynomial(1, 1) == polynomial(0))

    def testDivision(self):
        """Check div and mod of polynomials
        """
        p = polynomial(2, 0, 1, 1)

        d = polynomial(2, 0, 1)
        self.assertEqual(p // d, polynomial(1, 1))
        self.assertEqual(p % d, polynomial(0, -2))

        q, r = p // d, p % d
        self.assertEqual(q, polynomial(1, 1))
        self.assertEqual(r, polynomial(0, -2))


class TestEvaluatePolynomial(TestCase):
    def setUp(self):
        self.p = polynomial(1, 2, 1)

    def testEvaluate(self):
        """Check: evaluation of a polynomial
        """
        self.assertTrue(hasattr(self.p, 'evaluate'))
        self.assertTrue(self.p.evaluate(0) == 1)
        self.assertTrue(self.p.evaluate(1) == 4)


class TestFormalDerivative(TestCase):
    def testFormalDerivative(self):
        """Check: formal derivative of a polynomial
        """
        p = polynomial(1, 2, 3, 4)
        assert p.formal_derivative() == polynomial(2, 6, 12)


class TestStrRepr(TestCase):
    def setUp(self):
        self.p0 = polynomial(0)
        self.p1 = polynomial(-1, -1)
        self.p2 = polynomial(1, 0, -2, indeterminate='z')

    def testStr(self):
        """Extended check: printable representation
        """
        self.assertTrue(str(self.p0) == '0')
        self.assertTrue(str(self.p1) == '-1 - X')
        self.assertTrue(str(self.p2) == '1 - 2z^2')

    def testRepr(self):
        """Extended check: evaluable representation
        """
        self.assertEqual(self.p0, eval(repr(self.p0)))
        self.assertEqual(self.p1, eval(repr(self.p1)))
        self.assertEqual(self.p2, eval(repr(self.p2)))


class TestGCD(TestCase):
    def setUp(self):
        self.p1 = polynomial(1, 2, 1)
        self.p2 = polynomial(1, 3, 3, 1)
        self.p3 = polynomial(2, 0, 2)

    def testGCD_1(self):
        """Extended check: GCD of polynomials over Z
        """
        self.assertEqual(self.p1.gcd(polynomial(1)), polynomial(1))
        self.assertEqual(self.p1.gcd(self.p2), polynomial(1, 2, 1))

        with self.assertRaises(ValueError):
            self.p2.gcd(self.p3)


if __name__ == '__main__':
    run_tests()
