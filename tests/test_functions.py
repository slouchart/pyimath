from unittest import TestCase
from unittest import main as run_tests
from random import randint


from pyimath.functions import *


class TestBinCoeff(TestCase):
    def test(self):
        r = bincoeff(5)
        self.assertListEqual(r, [1, 5, 10, 10, 5, 1])


class TestFactor(TestCase):
    def testPrime(self):
        r = factor(2137)
        self.assertListEqual(r, [(2137, 1,)])

    def testPowerOfTwo(self):
        r = factor(2**9)
        self.assertListEqual(r, [(2, 9,)])

    def testRandom(self):
        for _ in range(15):
            n = randint(200, 500)
            r = factor(n)
            self.assertEqual(n, mul_factor(r))

    def testPrimes(self):
        p = primes(100)
        self.assertEqual(len(p), 25)
        self.assertListEqual([47, 53, 59, 61, 67], p[14:19])
        self.assertEqual(p[-1], 97)
        self.assertEqual(p[0], 2)


if __name__ == '__main__':
    run_tests()
