Module pyimath.functions
========================

Variables
---------

`small_primes`
:   A tuple that defines all primes up to 43

Functions
---------

    
`bincoeff(n, k=None)`
:   Computes the binomial coefficients of `(1 + x)^n` or returns the k-th coefficient

    
`comb(n, k)`
:   Defines `C(n, k)` as the number of ways to pick `k` items among `n`

    
`factor(n)`
:   Computes the prime factorization of an integer the hard way

    
`gcd(a, b)`
:   Computes the GCD of two operands recursively where
    `a` and `b` shall be of the same type or of compatible types.
    
    Mimics the behavior of `math.gcd` for zero operands.

    
`maybe_prime(n, k=3)`
:   Return `True` if `n` passes `k` rounds of the Miller-Rabin primality
    test (and is probably prime). Return `False` if `n` is proved to be
    composite.

    
`power(a, n)`
:   Return `a` to the n-th power, uses exponentiation by squaring i.e square and multiply

    
`primes(n_max=100)`
:   Implements the Eratosthene's sieve

    
`reduce_to_gcd(it)`
:   Computes the GCD of all items from an `Iterator`