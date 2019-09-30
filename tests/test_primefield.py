from unittest import TestCase
from unittest import main as run_tests


from imath.primefield import PrimeField


class TestAdditiveGroup(TestCase):
    def setUp(self):
        self.gr2 = PrimeField.additive_group_representation(2)
        self.gr3 = PrimeField.additive_group_representation(3)
        self.gr5 = PrimeField.additive_group_representation(5)
        self.gr7 = PrimeField.additive_group_representation(7)

    def testGR2(self):
        self.assertEqual(self.gr2, [0, 1])
        self.assertTrue(PrimeField._pf_add(1, 1, self.gr2) == 0)
        self.assertTrue(PrimeField._pf_add(0, 1, self.gr2) == PrimeField._pf_add(1, 0, self.gr2) == 1)
        self.assertTrue(PrimeField._pf_add(0, 0, self.gr2) == 0)

    def testGR3(self):
        self.assertEqual(self.gr3, [-1, 0, 1])
        self.assertTrue(PrimeField._pf_add(1, -1, self.gr3) == 0)
        self.assertTrue(PrimeField._pf_add(1, 1, self.gr3) == -1)

    def testGR5(self):
        self.assertEqual(self.gr5, [-2, -1, 0, 1, 2])
        self.assertTrue(PrimeField._pf_add(2, 2, self.gr5) == -1)
        self.assertTrue(PrimeField._pf_add(-2, -2, self.gr5) == 1)
        self.assertTrue(PrimeField._pf_add(1, -2, self.gr5) == -1)
        for e in self.gr5:
            for f in self.gr5:
                self.assertTrue(PrimeField._pf_add(e, f, self.gr5) == PrimeField._pf_add(f, e, self.gr5))

    def testGR7(self):

        self.assertEqual(self.gr7, [-3, -2, -1, 0, 1, 2, 3])
        self.assertTrue(PrimeField._pf_add(2, 3, self.gr7) == -2)
        self.assertTrue(PrimeField._pf_add(3, 1, self.gr7) == -3)
        self.assertTrue(PrimeField._pf_add(-1, -3, self.gr7) == 3)


class TestMultiplicativeGroup(TestCase):
    def setUp(self):
        self.gr2 = PrimeField.additive_group_representation(2)
        self.gr3 = PrimeField.additive_group_representation(3)
        self.gr5 = PrimeField.additive_group_representation(5)
        self.gr7 = PrimeField.additive_group_representation(7)

        self.mulgr3 = PrimeField.multiplicative_group_representation(self.gr3)
        self.mulgr5 = PrimeField.multiplicative_group_representation(self.gr5)
        self.mulgr7 = PrimeField.multiplicative_group_representation(self.gr7)

    def testMULGR3(self):
        self.assertTrue(PrimeField._pf_mul(-1, -1, self.mulgr3) == 1)

    def testMULGR5(self):
        self.assertTrue(PrimeField._pf_mul(2, 2, self.mulgr5) == -1)
        self.assertTrue(PrimeField._pf_mul(2, -2, self.mulgr5) == 1)
        self.assertTrue(PrimeField._pf_ext_mul(2, -2, self.gr5) == 1)

    def testMULGR7(self):
        self.assertTrue(PrimeField._pf_mul(2, 2, self.mulgr7) == -3)
        self.assertTrue(PrimeField._pf_mul(2, 3, self.mulgr7) == -1)
        self.assertTrue(PrimeField._pf_mul(-3, -3, self.mulgr7) == 2)

        self.assertTrue(PrimeField._pf_div(-3, 2, self.mulgr7) == 2)
        self.assertTrue(PrimeField._pf_div(-2, 3, self.mulgr7) == -3)


class TestFieldOps(TestCase):
    def setUp(self):
        self.f5 = PrimeField(5)

    def TestAllOperators(self):
        f5 = self.f5
        u = f5(2)
        v = f5(-1)

        self.assertTrue(u + v == 1)
        self.assertTrue(u * u == -1)
        self.assertTrue(-u * -u == -1)
        self.assertTrue(u != v)
        self.assertTrue(+u == u)

        self.assertTrue(u == f5(2))
        self.assertTrue(v == f5(-1))

        self.assertTrue(u * v == f5(-2))
        self.assertTrue(4 * u == u + u + u + u)

        self.assertTrue(u ** 4 == 1)

        self.assertTrue(f5(2) / f5(-2) == f5(-1))
        self.assertTrue(f5(0) / f5(1) == f5(0))

        self.assertTrue(abs(f5(-2)) == f5(2))


if __name__ == '__main__':
    run_tests()
