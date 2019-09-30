from unittest import TestCase
from unittest import main as run_tests


from imath.gaussint import GaussianIntegerRing


class TestGaussInt(TestCase):
    def setUp(self):
        self.fg = GaussianIntegerRing()

    def testAddSub(self):
        fg = self.fg
        elt = fg.element(1, 1)
        z0 = fg(0, 0)
        self.assertTrue(fg.null == z0)
        self.assertTrue(elt + fg.zero == elt)

        elt2 = 1 + elt + (0, 2)
        self.assertTrue(elt2 == elt + (1, 2) == (1, 2) + elt)

        self.assertTrue(elt2-elt2 == fg.zero)

        self.assertTrue(fg(2, 3) - fg(1, 0) == -(1 - fg(2, 3)) == fg(1, 3))

    def testStr(self):
        fg = self.fg
        elt2 = fg(2, 3)
        self.assertEqual(str(elt2), '2 + 3i')
        self.assertEqual(str(fg(0, -1)), '-i')

    def testParse(self):
        elt = self.fg.parse('-1 - 3i')
        self.assertEqual(elt, self.fg(-1, -3))

        with self.assertRaises(SyntaxError):
            self.fg.parse('-X + 5')

        self.assertEqual(self.fg.parse(''), self.fg.zero)

        with self.assertRaises(ValueError):
            self.fg.parse('-6 + 87i^2')

    def testMul(self):
        fg = self.fg

        elt2 = fg(2, 3)
        self.assertTrue(elt2 * fg.neutral == elt2 * fg.one == elt2)

        self.assertTrue(elt2 == fg.one * (2, 3) == (2, 3) * fg.neutral)
        self.assertTrue(elt2 == 1 * fg(2, 3) == fg(2, 3) * 1)

    def testRepr(self):
        elt2 = self.fg(2, 3)

        try:
            elt = eval(repr(elt2))
        except BaseException as e:
            self.fail(f'eval({repr(elt2)}) unexpectedly raised an exception {e}')

        self.assertEqual(elt, elt2)

    def testConjugate(self):
        fg = self.fg
        self.assertTrue(fg(1, 0).conjugate == fg.one == ~fg(1, 0))
        self.assertTrue(fg(0, 1).conjugate == fg(0, -1))

    def testDiv(self):
        fg = self.fg
        n = fg(4, 2)
        d = fg(2, 2)
        q, r = divmod(n, d)

        self.assertTrue(q == fg(2, -1) == n // d)
        self.assertTrue(r == fg(-2, 0) == n % d)
        self.assertTrue((q, r) == fg.divmod(n, d))
        with self.assertRaises(TypeError):
            _ = n / d

    def testPow(self):
        fg = self.fg
        n = fg(3, 0)
        assert n**2 == fg(9, 0)

        n = fg(1, 1)
        assert n**2 == fg(0, 2)


if __name__ == '__main__':
    run_tests()
