from typing import Iterator, Tuple, Union, Any, List, Dict
import random
import operator

from pyimath.functions import maybe_prime, power
from pyimath.polynomial import Polynomial, symbolic_polynomial


AdditiveGroup = List[int]
MultiplicativeGroup = Dict[Tuple[int, int], Union[int, List[int]]]


__all__ = ['PrimeField', 'PFElement']


class PrimeField:
    """Defines a prime field definition of characteristic P.

       The field elements are not represented by integer modulo P
       but rather by signed integers in the range `-(P-1)/2..(P-1)/2`
       """
    def __init__(self, prime: int):
        self.characteristic = prime
        self.additive_group = self._additive_group_representation(prime)
        self.multiplicative_group = self._multiplicative_group_representation(self.additive_group)

    def add(self, a: 'PFElement', b: 'PFElement') -> 'PFElement':
        """Returns the sum of two elements
        """
        return self(self._pf_add(a.value, b.value, self.additive_group))

    def additive_inverse(self, a: 'PFElement') -> 'PFElement':
        """Returns the additive inverse of an element
        """
        return self(self._pf_additive_inverse(a.value, self.additive_group))

    def div(self, a: 'PFElement', b: 'PFElement') -> 'PFElement':
        """Returns the quotient of two elements
        """
        return self(self._pf_div(a.value, b.value, self.multiplicative_group))

    def divmod(self, a: 'PFElement', b: 'PFElement') -> Tuple['PFElement', 'PFElement']:
        """Kept for interface purposes. Always returns `(a / b, self.zero)`
        """
        return self.div(a, b), self.mod(a, b)

    def element(self, n: Any) -> 'PFElement':
        """Casts an integer into an element of the field
        """
        if not int(n) in self.additive_group:
            raise ValueError(f'{n} does not belong to the additive group of {self}')
        return PFElement(self, int(n))

    def ext_mul(self, n: int, a: 'PFElement') -> 'PFElement':
        """Returns the n-th iterated sum of an element
        i.e  `a + ... + a, n times`"""
        return self(self._pf_ext_mul(n, a.value, self.additive_group))

    def floor_div(self, a: 'PFElement', b: 'PFElement') -> 'PFElement':
        """Same as `self.div`
        """
        return self.div(a, b)

    def frobenius_reciprocal(self, a: 'PFElement') -> 'PFElement':
        """Returns the Frobenius reciprocal of an element. As a remainder, the Frobenius automorphism is `a -> a^q`
        where `q` is the characteristic of the field. In prime fields, it is the identity.
        """
        assert a in self
        return a

    def generate_irreducible_polynomial(self, degree: int, max_retries: int = 15) -> Polynomial:
        """Returns an irreducible polynomial of a given degree over the base field. This algorithm is not deterministic
        and may raise exceptions. `degree` is the degree of the irreducible polynomial and `max_retries` sets
        the maximum number of attempts.
        """
        max_retries = max(degree // 2, max_retries)
        tries = 0
        retries = 0
        while retries <= max_retries:
            while tries <= degree:
                coefficients = [self.random_element() for _ in range(degree)]
                p = self.polynomial(coefficients)
                p += p.monic(degree)

                if p.is_irreducible:
                    return p
                tries += 1
            retries += 1
        err_msg = f'Could not find an irreducible polynomial of degree {degree} over {self} in {tries*retries} attempts'
        raise RuntimeError(err_msg)

    def linear_polynomial(self, e: 'PFElement') -> Polynomial:
        """Returns the polynomial `X + (-e)`
        """
        poly = self.polynomial(-e)
        poly += poly.monic(1)
        return poly

    def mul(self, a: 'PFElement', b: 'PFElement') -> 'PFElement':
        """Returns the product of two elements
        """
        return self(self._pf_mul(a.value, b.value, self.multiplicative_group))

    def multiplicative_inverse(self, a: 'PFElement') -> 'PFElement':
        """Returns the multiplicative inverse of an element
        """
        return self(self._pf_multiplicative_inverse(a.value, self.multiplicative_group))

    def mod(self, *_) -> 'PFElement':
        """Kept for interface purposes. Always return `self.zero`
        """
        return self.zero

    @property
    def neutral(self) -> 'PFElement':
        """Same as `self.one`
        """
        return self.one

    @property
    def null(self) -> 'PFElement':
        """Same as `self.zero`
        """
        return self.zero

    @property
    def one(self) -> 'PFElement':
        """Returns the neutral of the multiplicative group
        """
        return self(1)

    def parse_poly(self, expr: str) -> Polynomial:
        """Returns a polynomial from its symbolic expression
        """
        return symbolic_polynomial(expr, self)

    def polynomial(self, *args, indeterminate: str = 'X') -> Polynomial:
        """Returns a polynomial from the arguments
        """
        return Polynomial([self.element(c) for c in args], base_field=self, indeterminate=indeterminate)

    def pow(self, a: 'PFElement', n: int) -> 'PFElement':
        """Return the n-th power of an element
        """
        res = power(a, n)
        if not isinstance(res, PFElement):
            return self.element(res)
        else:
            return res

    def random_element(self) -> 'PFElement':
        """Returns an element of the field at random
        """
        return random.choice(list(iter(self)))

    def random_polynomial(self, degree: int) -> Polynomial:
        """Returns a random polynomial of a given degree
        """
        p = self.polynomial(*[self.random_element() for _ in range(0, degree)])
        p += p.monic(degree)
        return p

    @property
    def zero(self) -> 'PFElement':
        """Returns the neutral of the additive group
        """
        return self(0)

    def __call__(self, n: int) -> 'PFElement':
        """Syntactic sugar for self.element(e)
        """
        return self.element(n)

    def __contains__(self, item: Union['PFElement', int]) -> bool:
        if isinstance(item, PFElement):
            assert item.field == self
            return True
        elif isinstance(item, int):
            return self(item) in self.additive_group
        else:
            return False

    def __eq__(self, other: 'PrimeField') -> bool:
        assert isinstance(other, self.__class__)
        return other.characteristic == self.characteristic

    def __iter__(self) -> Iterator['PFElement']:
        """Allows iteration over all the elements
        """
        return (self(n) for n in self.additive_group)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.characteristic})'

    def __str__(self) -> str:
        return f'Prime field of characteristic {self.characteristic}'

    """LOW LEVEL FIELD OPERATIONS"""

    @staticmethod
    def _additive_group_representation(p: int) -> AdditiveGroup:
        assert p >= 2  # indeed, p must be a prime
        assert maybe_prime(p, 3)
        res = [0, 1]
        if p == 2:
            return res

        res = [-1] + res
        for elt in range(2, (p + 1) // 2):
            res.append(elt)
            res.insert(0, -elt)
        return res

    @staticmethod
    def _pf_add(a: int, b: int, gr: AdditiveGroup) -> int:
        assert a in gr and b in gr

        if a == 0:
            return b

        if b == 0:
            return a

        if a == -b:
            return 0

        p = len(gr)
        # special cas for F2, this field isn't symmetric
        if p == 2:
            if a == 1 and b == 1:
                return 0

        bound = p - 1
        pos = gr.index(a)
        offset = b
        if offset > 0:
            if pos + offset > bound:
                return gr[(pos + offset) - bound - 1]
            else:
                return gr[pos + offset]
        else:
            if pos + offset < 0:
                return gr[bound + (pos + offset) + 1]
            else:
                return gr[pos - abs(offset)]

    @staticmethod
    def _pf_additive_inverse(n: int, gr: AdditiveGroup) -> int:
        p = len(gr)
        if p == 2:
            return n
        else:
            return -n

    @staticmethod
    def _multiplicative_group_representation(gr: AdditiveGroup) -> MultiplicativeGroup:
        def uni_deg_one_poly(elt):
            return [elt - 1]

        def mul_uni_deg_one_poly(p1, p2, gr_add, gr_mul):
            return [PrimeField._pf_mul(p1[0], p2[0], gr_mul), PrimeField._pf_add(p1[0], p2[0], gr_add)]

        def eval_uni_poly_at_1(p0, gr_add):
            s = 1
            for elt in p0:
                s = PrimeField._pf_add(s, elt, gr_add)
            return s

        p = len(gr)  # characteristic of the field
        res = {(1, 1): 1}
        reciprocals = {k: None for k in gr if k != 0}
        reciprocals[1] = 1

        if p > 2:
            # step 0, initialize the matrix through the action of group {-1, 1} over the field
            for e in gr:
                if e != 0:
                    res[(1, e)] = res[(e, 1)] = e
                    res[(e, -1)] = res[(-1, e)] = -e
                    res[(1, -e)] = res[(-e, 1)] = -e
                    res[(-e, -1)] = res[(-1, -e)] = e
            res[(-1, -1)] = 1
            reciprocals[-1] = -1

            if p > 3:

                for e in range(2, (p + 1) // 2):
                    # step 1, compute squares
                    pn = uni_deg_one_poly(e)
                    v = eval_uni_poly_at_1(mul_uni_deg_one_poly(pn, pn, gr, res), gr)
                    res[(e, e)] = res[(-e, -e)] = v
                    res[(-e, e)] = res[(e, -e)] = -v
                    if v == 1:
                        reciprocals[e] = e
                    if -v == 1:
                        reciprocals[e] = -e
                        reciprocals[-e] = e

                    # step 2, compute lateral products
                    for f in range(e + 1, (p + 1) // 2):
                        pn2 = uni_deg_one_poly(f)
                        v = eval_uni_poly_at_1(mul_uni_deg_one_poly(pn2, pn, gr, res), gr)
                        res[(e, f)] = res[(f, e)] = v
                        res[(-e, f)] = res[(f, -e)] = -v
                        res[(e, -f)] = res[(-f, e)] = -v
                        res[(-e, -f)] = res[(-f, -e)] = v
                        if v == 1:
                            reciprocals[e] = f
                            reciprocals[f] = e
                            reciprocals[-e] = -f
                        if -v == 1:
                            reciprocals[-e] = f
                            reciprocals[f] = -e
                            reciprocals[-f] = e
                            reciprocals[e] = -f

        assert len(res) == (p - 1) ** 2

        res[(0, 0)] = reciprocals
        return res

    @staticmethod
    def _pf_mul(a: int, b: int, gr: MultiplicativeGroup) -> int:
        if a == 0 or b == 0:
            return 0

        assert (a, b) in gr, f'{(a, b)} not in {gr}'

        if a == 1:
            return b

        if b == 1:
            return a

        if a == -1:
            return -b

        if b == -1:
            return -a

        return gr[(a, b)]

    @staticmethod
    def _pf_ext_mul(n: int, a: int, gr: AdditiveGroup) -> int:
        if a == 0:
            return 0
        res = a
        for i in range(1, n):
            res = PrimeField._pf_add(res, a, gr)
        return res

    @staticmethod
    def _pf_multiplicative_inverse(a: int, gr: MultiplicativeGroup) -> int:
        if a == 0:
            raise ZeroDivisionError

        return gr[(0, 0)][a]

    @staticmethod
    def _pf_div(a: int, b: int, gr: MultiplicativeGroup) -> int:
        _1_b = PrimeField._pf_multiplicative_inverse(b, gr)
        return PrimeField._pf_mul(a, _1_b, gr)


class PFElement:
    """Represents a single element from a prime field by, basically, duck-typing an integer
    """
    def __init__(self, field: PrimeField, n: int):
        self.field = field
        self.value = n
        assert self.value in self.field.additive_group, f'{self.value} is not a field element'

    @property
    def is_one(self) -> bool:
        """Return `True` if the element is the multiplicative neutral of its field, `False` otherwise
        """
        return self.field.one == self

    @property
    def is_zero(self) -> bool:
        """Returns `True` if the element is the additive neutral of its field, `False` otherwise
        """
        return self.field.zero == self

    def __abs__(self) -> 'PFElement':
        return self.field(abs(self.value))

    def __add__(self, other: Any) -> 'PFElement':
        if isinstance(other, self.__class__):
            return self.field.add(self, other)
        elif isinstance(other, int):
            return self.field.add(self, self.field(other))
        else:
            return operator.add(other, self)

    def __divmod__(self, other: Any) -> Tuple['PFElement', 'PFElement']:
        return self.__truediv__(other), self.field.zero

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, self.__class__):
            return other.value == self.value
        elif isinstance(other, int):
            return other == self.value
        else:
            return other == self

    def __floordiv__(self, other: Any) -> 'PFElement':
        if isinstance(other, int):
            other = self.field(other)
        return self.field.div(self, other)

    def __format__(self, format_spec: str = '') -> str:
        try:
            return format(int(self.value), format_spec)
        except ValueError:
            # invalid format specifier
            return format(int(self.value))

    def __gt__(self, other: Any) -> bool:
        return int(self.value) > int(other)

    def __hash__(self) -> int:
        return self.value

    def __iadd__(self, other: Any) -> 'PFElement':
        self.value = int(self.__add__(other))
        return self

    def __ifloordiv__(self, other: Any) -> 'PFElement':
        if isinstance(other, int):
            other = self.field(other)
        self.value = int(self // other)
        return self

    def __int__(self) -> int:
        return self.value

    def __invert__(self) -> 'PFElement':
        return self.field.frobenius_reciprocal(self)

    def __isub__(self, other: Any) -> 'PFElement':
        if isinstance(other, int):
            other = self.field(other)
        self.value = int(self - other)
        return self

    def __itruediv__(self, other: Any) -> 'PFElement':
        return self.__ifloordiv__(other)

    def __imul__(self, other: Any) -> 'PFElement':
        self.value = int(self.__mul__(other))
        return self

    def __lt__(self, other: Any) -> bool:
        return int(self.value) < int(other)

    def __mod__(self, other: Any) -> 'PFElement':
        return self.field.zero  # element of a field always has a mult. inv.

    def __mul__(self, other: Any) -> 'PFElement':
        if isinstance(other, self.__class__):
            return self.field.mul(self, other)
        elif isinstance(other, (int, float)):
            return self.field.ext_mul(other, self)
        else:
            return operator.mul(other, self)

    def __neg__(self) -> 'PFElement':
        return self.field.additive_inverse(self)

    def __pos__(self) -> 'PFElement':
        return self.field(self.value)

    def __pow__(self, e: int, modulo=None) -> 'PFElement':
        return self.field.pow(self, e)

    def __radd__(self, other: Any) -> 'PFElement':
        return self.__add__(other)

    def __repr__(self) -> str:
        return repr(int(self))

    def __rmul__(self, other: Any) -> 'PFElement':
        return self.__mul__(other)

    def __rsub__(self, other: Any) -> 'PFElement':
        if isinstance(other, int):
            other = self.field(other)
        return self.field.add(other, -self)

    def __rfloordiv__(self, other: Any) -> 'PFElement':
        if isinstance(other, int):
            other = self.field(other)
        return self.field.div(other, self)

    def __rtruediv__(self, other: Any) -> 'PFElement':
        return self.__rfloordiv__(other)

    def __str__(self) -> str:
        return str(int(self))

    def __sub__(self, other: Any) -> 'PFElement':
        return self.__add__(-other)

    def __truediv__(self, other: Any) -> 'PFElement':
        return self.__floordiv__(other)
