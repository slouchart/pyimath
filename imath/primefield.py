from typing import Iterator, List, Dict, Tuple, Union
import random


import operator
from imath.functions import maybe_prime, power
from imath.base import IntegralDomain
from imath.polynomial import Polynomial

AdditiveGroup = List[int]
MultiplicativeGroup = Dict[Tuple[int, int], Union[int, List[int]]]


class PrimeField(IntegralDomain):
    """Prime field definition of characteristic P
       The representation of the field elements are not integer modulo P
       but rather relative integers in the range -(P-1)/2..(P-1)/2"""
    def __init__(self, prime: int):
        super().__init__(prime)
        self.additive_group = additive_group_representation(prime)
        self.multiplicative_group = multiplicative_group_representation(self.additive_group)

    def element(self, n):
        """Casts an integer into an element of the field"""
        if not int(n) in self.additive_group:
            raise ValueError(f'{n} does not belong to the additive group of {self}')
        return PFElement(self, int(n))

    @property
    def zero(self):
        """returns the neutral of the additive group"""
        return self(0)

    @property
    def null(self):
        return self.zero

    @property
    def one(self):
        """returns the neutral of the multiplicative group"""
        return self(1)

    @property
    def neutral(self):
        return self.one

    def add(self, a, b):
        """Addition in the field"""
        return self(_pf_add(a.value, b.value, self.additive_group))

    def additive_inverse(self, a):
        """Returns -a """
        return self(_pf_additive_inverse(a.value, self.additive_group))

    def mul(self, a, b):
        """Multiplication in the field"""
        return self(_pf_mul(a.value, b.value, self.multiplicative_group))

    def ext_mul(self, n: int, a):
        """Multiplication of a by integer n <=> a+...+a n times"""
        return self(_pf_ext_mul(n, a.value, self.additive_group))

    def multiplicative_inverse(self, a):
        """Returns 1/a"""
        return self(_pf_multiplicative_inverse(a.value, self.multiplicative_group))

    def floor_div(self, a, b):
        return self.div(a, b)

    def mod(self, a, b):
        return self.zero

    def divmod(self, a, b):
        return self.div(a, b), self.mod(a, b)

    def div(self, a, b):
        """Exact division in the field"""
        return self(_pf_div(a.value, b.value, self.multiplicative_group))

    def pow(self, a, e: int):
        """Returns a**e using rapid exponentiation"""
        return power(a, e, self.zero, self.one)

    def __eq__(self, other):
        assert isinstance(other, self.__class__)
        return other.characteristic == self.characteristic

    def __call__(self, n: int):
        """Syntactic sugar for self.element(e)"""
        return self.element(n)

    def __iter__(self) -> Iterator:
        """Allows iteration over all the elements"""
        return (self(n) for n in self.additive_group)

    def __contains__(self, item):
        if isinstance(item, PFElement):
            assert item.field == self
            return True
        elif isinstance(item, int):
            return self(item) in self.additive_group
        else:
            return False

    def __str__(self) -> str:
        return f'Prime field of characteristic {self.characteristic}'

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.characteristic})'

    def rand_element(self):
        return random.choice(list(iter(self)))

    def polynomial(self, *args, indeterminate='X') -> Polynomial:
        """Polynomial factory"""
        return Polynomial([self.element(c) for c in args], base_field=self, indeterminate=indeterminate)

    def random_polynomial(self, degree):
        p = self.polynomial(*[self.rand_element() for _ in range(0, degree)])
        p += p.monic(degree)
        return p

    def linear_polynomial(self, e) -> Polynomial:
        """Returns the polynomial X + (-e)"""
        poly = self.polynomial(-e)
        poly += poly.monic(1)
        return poly

    def generate_irreducible_polynomial(self, degree, max_retries=15):
        """"""
        max_retries = max(degree // 2, max_retries)
        tries = 0
        retries = 0
        while retries <= max_retries:
            while tries <= degree:
                coefficients = [self.rand_element() for _ in range(degree)]
                p = self.polynomial(coefficients)
                p += p.monic(degree)

                if p.is_irreducible:
                    return p
                tries += 1
            retries += 1
        err_msg = f'Could not find an irreducible polynomial of degree {degree} over {self} in {tries*retries} attempts'
        raise RuntimeError(err_msg)

    def p_th_root(self, a):
        """Taking the p-th root for a is finding an element b as b^p = a where p is the field characteristic
           Because the multiplicative group of the field is cyclic of order p-1 if such an element b exists
           b^(p-1) = 1 => b^p = b thus a = b"""
        assert a in self
        return a


class PFElement:
    """Generic class that duck-types an integer
       Represents a single element from a prime field"""
    def __init__(self, field: PrimeField, n: int):
        self.field = field
        self.value = n
        assert self.value in self.field.additive_group, f'{self.value} is not a field element'

    def __pos__(self):
        return self.field(self.value)

    def __neg__(self):
        return self.field.additive_inverse(self)

    def __add__(self, other):
        if isinstance(other, self.__class__):
            return self.field.add(self, other)
        elif isinstance(other, int):
            return self.field.add(self, self.field(other))
        else:
            return operator.add(other, self)

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        self.value = int(self.__add__(other))
        return self

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        if isinstance(other, int):
            return self.field(other).__sub__(self).value
        return other.__sub__(self)

    def __isub__(self, other):
        if isinstance(other, int):
            other = self.field(other)
        self.value = int(self - other)
        return self

    def __mul__(self, other):
        if isinstance(other, self.__class__):
            return self.field.mul(self, other)
        elif isinstance(other, (int, float)):
            return self.field.ext_mul(other, self)
        else:
            return operator.mul(other, self)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __imul__(self, other):
        self.value = int(self.__mul__(other))
        return self

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return other.value == self.value
        elif isinstance(other, int):
            return other == self.value
        else:
            return other == self

    def __lt__(self, other):
        return int(self.value) < int(other)

    def __gt__(self, other):
        return int(self.value) > int(other)

    def __int__(self):
        return self.value

    def __pow__(self, e: int, modulo=None):
        return self.field.pow(self, e)

    def __floordiv__(self, other):
        if isinstance(other, int):
            other = self.field(other)
        return self.field.div(self, other)

    def __truediv__(self, other):
        return self.__floordiv__(other)

    def __rfloordiv__(self, other):
        if isinstance(other, int):
            other = self.field(other)
        return self.field.div(other, self)

    def __rtruediv__(self, other):
        return self.__rfloordiv__(other)

    def __ifloordiv__(self, other):
        if isinstance(other, int):
            other = self.field(other)
        self.value = int(self // other)
        return self

    def __itruediv__(self, other):
        return self.__ifloordiv__(other)

    def __mod__(self, other):
        return self.field.zero  # element of a field always has a mult. inv.

    def __divmod__(self, other):
        return self.__truediv__(other), self.field.zero

    def __repr__(self):
        return repr(int(self))

    def __str__(self):
        return str(int(self))

    def __format__(self, format_spec=''):
        try:
            return format(int(self.value), format_spec)
        except ValueError:
            # invalid format specifier
            return format(int(self.value))

    def __abs__(self):
        return self.field(abs(self.value))

    def __hash__(self):
        return self.value

    @property
    def is_zero(self):
        return self.field.zero == self

    @property
    def is_one(self):
        return self.field.one == self


"""LOW LEVEL FIELD OPERATIONS"""

"""Note on the principles of operations

Usually addition and multiplication over Z/pZ
are defined from the addition and multiplication over Z modulo p
thus the classes of integers are 0, 1 ... p-1
Because we want to ease the arithmetic, we require our elements
to be in the range -(P-1)/2..(P+1)/2 for any P>2
For P=2, there is no change the operations are trivial anyway

So, how we do that?

For the addition, we use the group axioms and its cyclicity
i.e for a in the group 1+...+1 (p times) = 0
e.g for P=7, we map the elements 4, 5 and 6
respectively to -3, -2 and -1 because
1 + 1 + 1 + 1 = 4 and 1 + 1 + 1 = 3 
=> 1 + ... + 1 (7 times) = 0 => 4 = -3

For the multiplication, we need a table (a,b) = a*b
we construct this table using the following axioms
    a * 1 = 1 * a (note: this is enough for P = 2)
    a * -1 = -1 * a = -a
    This is enough for P = 3

when P > 3, we compute the squares of 2..(p-1)/2
using the distributivity of + over *
e.g 2 * 2 (F5) = (1 + 1)(1 + 1) = 1 + 1 + 1 + 1 = -1
then, the abelian group axioms gives the value of -a * -a
and -a * a
-a * -a = -1 * a * -1 * a  = a * a
-a * a = -1 * a * a = - a*a

while we are at it, once a*a has been computed
we can compute the products of a and all elements > a
a(a+1) = a*a + a
a(a+1+1) = a*a + a + a
and so on

The other products are computed using the commutativity
and the associativity of an abelian group

There is certainly a way to optimize this but it does not
seem to be trivial.

For instance, the only product we need to know for F5
is 2 * 2 = -2 * -2 = -(-2 * 2)

For F7, there are 2*2, 2*3 and 3*3
For F11, there are 2*2, 2*3, 2*4, 2*5, 3*3, 3*4, 3*5,
4*4, 4*5 and 5*5

The gain is the following:
the table takes O(p^2) space because its a (p-1) * (p-1)
matrix.
We note k = (p-1)/2 - 1 = (p-3)/2
The minimal number of products is then k(k+1)/2
The order is the same: O(p^2)
The only difference is a reduction of the overall size by 4
as p grows
"""


def additive_group_representation(p: int) -> AdditiveGroup:
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


def _pf_additive_inverse(n: int, gr: AdditiveGroup) -> int:
    p = len(gr)
    if p == 2:
        return n
    else:
        return -n


def multiplicative_group_representation(gr: AdditiveGroup) -> MultiplicativeGroup:
    def uni_deg_one_poly(elt):
        return [elt - 1]

    def mul_uni_deg_one_poly(p1, p2, gr_add, gr_mul):
        return [_pf_mul(p1[0], p2[0], gr_mul), _pf_add(p1[0], p2[0], gr_add)]

    def eval_uni_poly_at_1(p0, gr_add):
        s = 1
        for elt in p0:
            s = _pf_add(s, elt, gr_add)
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


def _pf_ext_mul(n: int, a: int, gr: AdditiveGroup) -> int:
    if a == 0:
        return 0
    res = a
    for i in range(1, n):
        res = _pf_add(res, a, gr)
    return res


def _pf_multiplicative_inverse(a: int, gr: MultiplicativeGroup) -> int:
    if a == 0:
        raise ZeroDivisionError

    return gr[(0, 0)][a]


def _pf_div(a: int, b: int, gr: MultiplicativeGroup) -> int:
    _1_b = _pf_multiplicative_inverse(b, gr)
    return _pf_mul(a, _1_b, gr)


def _pf_pow(*_) -> int:
    raise DeprecationWarning
