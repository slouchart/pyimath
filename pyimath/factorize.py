__all__ = [
    'factorize',
    'Factorization',
    'Factor',
]


from collections import namedtuple
from typing import Sequence, Tuple

from pyimath.polynomial import Polynomial
from pyimath.functions import gcd
from pyimath.annotations import BaseField, BaseNumber


def factorize(p: Polynomial) -> 'Factorization':
    """
    Main entry point of the module, provide an instance of the Factorization class as a place holder
    to call any factorization algorithms : `square_free`, `distinct_degree`, `equal_degree` or `cantor_zassenhaus`.

    Accepts an instance of `Polynomial` (parameter `p`) and returns a instance of `Factorization`
    """
    return Factorization(p.base_field, p)


class Factor(namedtuple('Factor', 'value, multiplicity, max_degree', defaults=(0,))):
    """Defines a `Factor` class that wraps a `Polynomial` along with a `multiplicity`.

    Works a `namedtuple`with the following fields:

    * `value`: an instance of `Polynomial`
    * `multiplicity`: an integer representing the multiplicity of this factor in a factorization
    * `max_degree`: for an reducible factor, represents the maximum degree of its irreducible factors
    """
    __slots__ = ()

    @property
    def is_irreducible(self) -> bool:
        """returns `True` if the `Polynomial` referenced by `value`is irreducible otherwise returns `False`"""
        return self.value.is_irreducible

    def __str__(self) -> str:
        if self.multiplicity > 1:
            s = f'({self.value})**{self.multiplicity}'
        else:
            s = f'({self.value})'
        if self.is_irreducible:
            s += ' IRREDUCIBLE'
        return s


class Factorization:
    """
    Service provider for polynomial factorization.

    The main service is the full factorization by the Cantor-Zassenhaus algorithm

    Never use it directly, subclass it if you want but in any case, please use `factorize(poly).<method>()`
    """

    def __init__(self, base_field: BaseField, poly: Polynomial):
        self.base_field = base_field
        self.poly = poly

    def factors_product(self, factors: Sequence[Factor]) -> Polynomial:
        """
        Converts a list of factors into a polynomial. The `factors` parameter must be an iterable of `Factor` instances.
        Returns a `Polynomial` that is the product of these factors along with their multiplicity
        """
        bf = self.base_field
        if len(factors) == 0:
            res = bf.polynomial(bf.zero)
        else:
            res = bf.polynomial(bf.one)
            for fct in factors:
                res *= (fct.value ** fct.multiplicity)

        return res

    def square_free(self) -> Tuple[Polynomial, Sequence[Factor]]:
        """Square free factorization of a monic polynomial over a finite field of prime characteristic.
           Returns the square free part and a partial or complete factorization"""

        f = self.poly.copy.make_monic()
        q = f.base_field.characteristic

        factors = []
        # step 2: identify all factors of multiplicity m not divisible by q
        if f.formal_derivative() != f.null:
            g = gcd(f, f.formal_derivative())
            if not g.is_abs_unit:
                w = f / g
                i = 1
                while not w.is_abs_unit:
                    y = w.gcd(g)
                    factor = w / y
                    if factor.degree > 0:
                        factors.append(Factor(factor.make_monic(), i, 0))

                    w, g = y, g / y
                    i += 1
            else:
                return f, factors
        else:
            g = f.copy

        if not g.is_abs_unit:
            # step 2: g is the product of all factors of multiplicity divisible by q
            g = ~g  # frobenius reciprocal

            # we try another sqf
            sqf, fct = factorize(g).square_free()
            factors += [Factor(sf.value.make_monic(), sf.multiplicity * q, 0) for sf in fct]
            factors += [Factor(sqf.make_monic(), q, 0)]

        # step 3: compute square free part
        sqf = f.copy
        for fact in factors[:]:
            if fact.multiplicity > 1:
                sqf /= fact.value ** fact.multiplicity
            else:
                factors.remove(fact)

        return sqf, factors

    def distinct_degree(self) -> Sequence[Factor]:
        """Distinct degree factorisation of a square free monic polynomial.
        Returns a list of Factors with a multiplicity of 1 and the maximum degree of any irreducible factor"""
        f = self.poly.copy.make_monic()
        factors = []

        fp = f.copy
        q = f.base_field.characteristic
        i = 1
        while fp.degree >= 2 * i:
            g = gcd(fp, f.monic(q ** i) - f.monic(1))
            if not g.make_monic().is_abs_unit:
                factors.append(Factor(g.make_monic(), 1, i))
            fp /= g
            i += 1

        if not fp.is_abs_unit and fp.degree > 0:
            factors.append(Factor(fp.make_monic(), 1, fp.degree))
        if len(factors) == 0:
            factors.append(Factor(f, 1, f.degree))

        return factors

    def equal_degree(self, nb_factors: int, max_degree: int) -> Sequence[Factor]:
        """Cantor–Zassenhaus factorization algorithm.
        Input: Over a finite field Fq of odd order q,
               p is a monic square free polynomial in Fq[x] of degree n = rd,
               which has r ≥ 2 irreducible factors each of degree d
        Output: The set of monic irreducible factors of f."""

        factors = []
        r, d = nb_factors, max_degree
        f = self.poly.copy.make_monic()
        q = self.poly.base_field.characteristic

        assert f.degree == r * d, f'Degree of {str(f)} mismatches input parameters {nb_factors}, {max_degree}'
        assert r >= 2

        max_retries = 2 * r * d
        found_factors, nb_retries = 0, 0
        while found_factors < r and nb_retries < max_retries:
            # pick a random polynomial a of degree < r * d
            a = f.base_field.random_polynomial(f.degree - 1)
            g = gcd(f, a).make_monic()
            if not g.is_unit and g not in factors:
                factors.append(g)
                found_factors += 1
                continue
            else:
                m = (q ** d - 1) // 2
                a_p = a ** m + a.unit
                if not a_p.is_null:
                    g = gcd(f, a_p).make_monic()
                    if not g.is_unit and g != f and g not in factors:
                        factors.append(g)
                        found_factors += 1
                        continue
            nb_retries += 1
        else:
            if nb_retries >= max_retries:
                raise RuntimeError(f'unable to find {r} degree {d} factors for {str(f)}')

        return [Factor(g, 1, 0) for g in factors]

    def cantor_zassenhaus(self) -> Tuple[Sequence[Factor], BaseNumber]:
        """Full factorisation algorithm for any polynomial over a finite field.
           Returns a tuple containing a constant term and a list of factors with their multiplicity"""

        irreducible_factors = []
        factors_to_consider = []
        # make monic polynomial
        if not self.poly.is_monic:
            f = self.poly.make_monic()
            constant_term = self.poly.leading
        else:
            f = self.poly.copy
            constant_term = self.poly.base_field.one

        # attempt a square free factorisation
        sqf, multiple_factors = factorize(f).square_free()
        # sqf may be irreducible or not
        # each factors may be irreducible or not
        if not sqf.is_unit:
            factors_to_consider.append(Factor(sqf, 1, 0))
        factors_to_consider += multiple_factors

        while len(factors_to_consider) > 0:  # we keep on factorizing while some reducible factors remain
            mfct = factors_to_consider.pop()
            if mfct.is_irreducible:
                irreducible_factors.append(mfct)
            elif mfct.max_degree == 0:
                subfactors = factorize(mfct.value).distinct_degree()
                for s in subfactors:
                    factors_to_consider.append(Factor(s.value, mfct.multiplicity, s.max_degree))
            else:
                d = int(mfct.max_degree)
                assert mfct.value.degree % d == 0
                r = mfct.value.degree // d
                factors_to_consider += factorize(mfct.value).equal_degree(r, d)

        return irreducible_factors, constant_term
