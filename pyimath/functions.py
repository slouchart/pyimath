from random import randrange
from functools import reduce
from math import gcd as __internal_gcd
from typing import Union, List, Tuple, Iterator

__all__ = [
    'comb',
    'bincoeff',
    'primes',
    'factor',
    'gcd',
    'reduce_to_gcd',
    'power',
    'maybe_prime',
    'small_primes',
]


def bincoeff(n: int, k: int = None) -> Union[int, List[int]]:
    """Computes the binomial coefficients of `(1 + x)^n` or returns the k-th coefficient
    """
    if k is not None:
        return comb(n, k)
    else:
        result = []
        for i in range(0, n + 1):
            result.append(comb(n, i))
        return result


def comb(n: int, k: int) -> int:
    """Defines `C(n, k)` as the number of ways to pick `k` items among `n`
    """
    if k > n:
        raise ValueError

    result = 1
    for i in range(n - k + 1, n + 1):
        result *= i
    for i in range(2, k + 1):
        result /= i

    return int(result)


def factor(n: int) -> List[Tuple[int, int]]:
    """Computes the prime factorization of an integer the hard way
    """
    if n <= 1:
        raise ValueError

    factors = list()

    ml = 0
    p = 2
    while n % p == 0:
        n //= p
        ml += 1
    if ml > 0:
        factors.append((p, ml,))

    p = 3
    while p ** 2 <= n:
        ml = 0
        while n % p == 0:
            n //= p
            ml += 1
        if ml > 0:
            factors.append((p, ml,))
        p += 2

    if n > 2:
        factors.append((n, 1,))

    return factors


def gcd(a, b):
    """Computes the GCD of two operands recursively where
    `a` and `b` shall be of the same type or of compatible types.

    Mimics the behavior of `math.gcd` for zero operands.
    """
    if a == 0:
        return b
    if b == 0:
        return a

    if isinstance(a, int) and isinstance(b, int):
        _gcd = __internal_gcd
    else:
        def _gcd(g, r):
            if r == 0:
                return g
            else:
                return _gcd(r, g % r)

    return _gcd(a, b)


def maybe_prime(n: int, k: int = 3) -> bool:
    """Return `True` if `n` passes `k` rounds of the Miller-Rabin primality
    test (and is probably prime). Return `False` if `n` is proved to be
    composite.
    """
    if n < 2:
        return False
    for p in small_primes:
        if n < p * p:
            return True
        if n % p == 0:
            return False
    r, s = 0, n - 1
    while s % 2 == 0:
        r += 1
        s //= 2
    for _ in range(k):
        a = randrange(2, n - 1)
        x = pow(a, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def mul_factor(factors: List[Tuple[int, int]]) -> int:
    """Computes an integer whose prime factorization is given
    """
    n = 1
    for f in factors:
        n *= f[0] ** f[1]
    return n


def power(a, n: int):
    """Return `a` to the n-th power, uses exponentiation by squaring i.e square and multiply
    """
    # Shortcuts are evaluated here to avoid code duplication

    if a == 0:
        if n > 0:
            return 0  # 0^n = 0 for n > 0

    if n == 0:
        return 1  # a^0 = 1 (0^0 = 1 is a convention in number theory)

    if n == 1:
        return a  # a^1 = a

    if a == 1:
        return a  # 1^n = 1, n integer

    def sqr_mul(x, e):
        if e == 1:
            return x
        elif e % 2 == 0:
            return sqr_mul(x * x, e // 2)
        elif e % 2 != 0 and e > 2:
            return x * sqr_mul(x * x, (e - 1) // 2)

    return sqr_mul(a, n)


def primes(n_max: int = 100) -> List[int]:
    """Implements the Eratosthene's sieve
    """
    if n_max < 2:
        raise ValueError

    t = list(range(2, n_max + 1))
    for i in t:
        for j in (k for k in t if k > i):
            if j % i == 0:
                t.remove(j)

    return sorted(t)


def reduce_to_gcd(it: Iterator):
    """ Computes the GCD of all items from an `Iterator`
    """
    return reduce(gcd, it, next(it))


small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,)
"""A tuple that defines all primes up to 43"""
