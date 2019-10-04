from itertools import product as cartesian_product
import operator

from typing import Optional, Sequence, MutableSequence, Iterator, Union, List, Tuple, Any, Collection


from pyimath.functions import maybe_prime
from pyimath.functions import power
from pyimath.polynomial import Polynomial
from pyimath.polyparse import symbolic_polynomial
from pyimath.primefield import PrimeField, PFElement

__all__ = ['Polynomial', 'FiniteField', 'PrimeField', 'finite_field', 'FFElement', 'PFElement']


BaseAtom = Union[int, PFElement]
Vector = Union[List[BaseAtom], MutableSequence[BaseAtom], Sequence[BaseAtom], Collection[BaseAtom]]


class FiniteField:
    """Represents a non-prime finite field"""
    def __init__(self,
                 prime: int,
                 dimension: int,
                 poly: Polynomial,
                 generator: Optional[Sequence[int]] = None,
                 root_symbol: Optional[str] = 'j'):
        """
        :param prime: the order of the prime field Fprime
        :param dimension: the dimension of the finite field seen as a vector space over its prime field
        :param poly: the minimal polynomial of the adjunct root j, must be irreducible over Fprime
        :param generator: an element of the field that generates the multiplicative group (tuple of int)
        :param root_symbol: for formatting elements
        """
        self._prime_field = PrimeField(prime)
        self.dimension = dimension
        self.base_polynomial = poly
        self.generator = generator
        self.root_symbol = root_symbol

        assert self.dimension >= 2
        assert self.base_polynomial.base_field == self.prime_field
        assert self.base_polynomial.degree == self.dimension
        assert self.base_polynomial.is_monic
        assert self.base_polynomial.is_irreducible

        self.root_powers = dict()
        self._compute_root_powers()

        self.generator_powers = dict()
        self.element_as_powers = dict()
        if self.generator is not None:
            self.generator = self.element(self.generator)
            self._check_generator_order()

        self._frobenius_map = self._compute_frobenius_map()

    def _compute_frobenius_map(self) -> MutableSequence['FFElement']:
        """The Frobenius automorphisme is defined by a -> a^p where p is the prime field characteristic
        This methods computes the inverse map of the automorphisme"""
        p_th_roots = []
        for b in self.basis:
            for e in self:
                if e ** self.characteristic == b:
                    p_th_roots.append(e)
        return p_th_roots

    def element_order(self, e: 'FFElement') ->int:
        """Returns the order of an element, that is the minimal integer k such as e^k = 1"""
        if self.has_valid_generator:
            return self.generator_powers[e]
        else:
            max_order = self.order - 1  # the order of the multiplicative group
            order = 1
            g = e
            while g != self.one and order < max_order:
                g *= e
                order += 1
            return order

    @property
    def frobenius_map(self) -> MutableSequence['FFElement']:
        """It's actually the linear map of the inverse function of Frobenius automorphism (see MISCELLANEOUS.md)"""
        return self._frobenius_map[:]

    def _check_generator_order(self):
        """computes the powers of the wanna-be generator and checks
        if the results form the multiplicative group. This method can be
        computer-intensive for large finite fields"""
        expected_order = self.order - 1  # the order of the multiplicative group
        order = 1
        g = self.generator
        powers = dict()
        powers[order] = g

        while g != self.one and order <= expected_order:
            g *= self.generator
            order += 1
            powers[order] = g

        if order != expected_order:
            raise ValueError(f'Element {self.generator} is not a generator for {self}')
        else:
            self.generator_powers = powers

        if self.has_valid_generator:
            for (o, e) in self.generator_powers.items():
                if e != self.zero:
                    self.element_as_powers[e] = o  # reverse map

    @property
    def has_valid_generator(self) -> bool:
        return len(self.generator_powers) > 0 and self.generator is not None

    @property
    def prime_field(self) -> PrimeField:
        return self._prime_field

    @property
    def basis(self) -> Iterator['FFElement']:
        return (self(0, 1) ** e for e in range(0, self.dimension))

    def _compute_root_powers(self):
        p = self.base_polynomial
        f = self.prime_field
        pr = p - p.monic(p.degree)
        p_neutral = f.polynomial(f.one)
        r = -pr
        e = self.dimension
        self.root_powers[e] = r
        while r != p_neutral and e < self.order - 1:
            r *= r.monic(1)
            e += 1
            if r.degree == p.degree:
                r -= p.monic(p.degree)
                r += -pr
            self.root_powers[e] = r

    def _safe_convert_vector(self, v: Vector) -> List[PFElement]:

        if isinstance(v, FFElement):
            assert v.field == self
            o = v.vector[:]
        else:
            o = []
            if not isinstance(v, (list, tuple)):
                v = [v]

            for c in v:
                if isinstance(c, int):
                    o.append(self.prime_field(c))
                else:
                    assert isinstance(c, PFElement)
                    o.append(c)

            for _ in range(self.dimension - len(v)):
                o.append(self.prime_field.zero)

        return o

    @property
    def order(self) -> int:
        """
        :return: The order of the field
        """
        return self.prime_field.characteristic ** self.dimension

    def element(self, v: Vector) -> 'FFElement':
        return FFElement(self, self._safe_convert_vector(v))

    @property
    def zero(self) -> 'FFElement':
        return self.element([self.prime_field.zero] * self.dimension)

    @property
    def one(self) -> 'FFElement':
        return self.element([self.prime_field.one] + [self.prime_field.zero] * (self.dimension - 1))

    @property
    def characteristic(self) -> int:
        return self._prime_field.characteristic

    @property
    def null(self) -> 'FFElement':
        return self.zero

    @property
    def neutral(self) -> 'FFElement':
        return self.one

    def __call__(self, *args) -> 'FFElement':
        return self.element([self.prime_field(c) for c in args])

    def __eq__(self, other: 'FiniteField') -> bool:
        assert isinstance(other, self.__class__)
        return other.prime_field == self.prime_field and other.dimension == self.dimension

    def __iter__(self) -> Iterator['FFElement']:
        """
        Iterates over all the elements
        :return: an iterator
        """
        iterables = [self.prime_field] * self.dimension
        return iter((self.element(list(v)) for v in cartesian_product(*iterables)))

    def __str__(self) -> str:
        return f'Finite field of order {self.prime_field.characteristic}^{self.dimension}'

    def __repr__(self) -> str:
        s = f'{self.__class__.__name__}('
        s += f'{self.prime_field.characteristic}, '
        s += f'{self.dimension}, '
        s += f'{repr(self.base_polynomial)}, '
        if self.has_valid_generator:
            s += f' generator={repr(self.generator.vector)})'
        else:
            s += ')'
        return s

    def add(self, a: 'FFElement', b: 'FFElement') -> 'FFElement':
        assert a.field == b.field == self
        s = [self.prime_field.zero] * self.dimension
        for i in range(0, self.dimension):
            s[i] = a[i] + b[i]
        return self.element(s)

    def additive_inverse(self, a: 'FFElement') -> 'FFElement':
        return self.element([-e for e in a.vector])

    def mul(self, a: 'FFElement', b: 'FFElement') -> 'FFElement':
        assert a.field == b.field == self

        if self.has_valid_generator:
            p1 = self.element_as_powers[a]
            p2 = self.element_as_powers[b]
            p = p1 + p2
            if p in self.generator_powers.keys():
                return self.generator_powers[p]
            else:
                assert p > self.order - 1
                return self.generator_powers[p - self.order + 1]

        else:
            p1 = self.polynomial_from_element(a)
            p2 = self.polynomial_from_element(b)
            p = (p1 * p2) % self.base_polynomial

            return self.element_from_polynomial(p)

    def ext_mul(self, n: int, a: 'FFElement') -> 'FFElement':
        assert a.field == self
        r = self.zero
        for _ in range(0, n):
            r += a
        return r

    def pow(self, a: 'FFElement', e: int) -> 'FFElement':
        assert a.field == self
        assert e >= 0
        res = power(a, e)
        if not isinstance(res, FFElement):
            return self(res)
        else:
            return res

    def multiplicative_inverse(self, a: 'FFElement') -> 'FFElement':

        if self.has_valid_generator:
            p = self.element_as_powers[a]
            return self.generator_powers[self.order - 1 - p]
        else:
            p_a = self.polynomial_from_element(a)
            p = self.base_polynomial

            if p_a == p or p_a == p_a.null:
                raise ZeroDivisionError

            neutral = p_a.unit
            t = p_a.null
            newt = neutral
            r = p
            newr = p_a
            while newr != newr.null:
                quotient = r / newr
                r, newr = newr, r - quotient * newr
                t, newt = newt, t - quotient * newt

            assert r.degree == 0
            # if r.degree > 0
            # Either p is not irreducible over the base_field or p_a is a multiple of p
            return self.element_from_polynomial((neutral / r) * t)

    def div(self, a: 'FFElement', b: 'FFElement') -> 'FFElement':
        if b.null:
            raise ZeroDivisionError
        if a.null:
            return self.zero

        return a * self.multiplicative_inverse(b)

    def floor_div(self, a: 'FFElement', b: 'FFElement') -> 'FFElement':
        return self.div(a, b)

    def mod(self, _, b: 'FFElement') -> 'FFElement':
        if b.null:
            raise ZeroDivisionError

        return self.zero

    def divmod(self, a: 'FFElement', b: 'FFElement') -> Tuple['FFElement', 'FFElement']:
        return self.div(a, b), self.mod(a, b)

    def format_element(self, e: 'FFElement', format_spec: str = '') -> str:

        if format_spec == 'short':
            s = ''
            if e.is_scalar:
                s = str(e.vector[0])
            else:
                for c in str(self.polynomial_from_element(e)):
                    if c != ' ':
                        s += c
                return f'({s})'
        else:
            s = str(self.polynomial_from_element(e))
        return s

    def polynomial(self, *args, indeterminate: Optional[str] = 'X') -> Polynomial:
        coefficients = []
        for c in args:
            if isinstance(c, (tuple, list)):
                coefficients.append(self.element(c))
            elif isinstance(c, FFElement):
                coefficients.append(c)
            else:
                # shall be an integer or, at least an atom (e.g. a PFElement)
                c = [c] + [self.prime_field.zero] * (self.dimension - 1)
                coefficients.append(self.element(c))

        return Polynomial(coefficients, base_field=self, indeterminate=indeterminate)

    def linear_polynomial(self, e: 'FFElement') -> Polynomial:
        p = self.polynomial(-e)
        p += p.monic(1)
        return p

    def check_irreducible_polynomial(self, p: Polynomial) -> bool:
        assert p.base_field == self
        return p.is_irreducible

    def frobenius_reciprocal(self, a: 'FFElement') -> 'FFElement':

        assert isinstance(a, FFElement)
        assert a in self
        r = self.zero
        for i, e in enumerate(a.vector):
            b = self._frobenius_map[i]
            r += self(e) * b
        return r

    def polynomial_from_element(self, e: 'FFElement') -> Polynomial:
        return Polynomial(e.vector, self.prime_field, indeterminate=self.root_symbol)

    def element_from_polynomial(self, p: Polynomial) -> 'FFElement':
        assert p.base_field == self.prime_field
        assert p.degree <= self.dimension - 1
        v = p.coefficients
        if p.degree < self.dimension - 1:
            v += [self.prime_field.zero] * (self.dimension - (1 + p.degree))
        return self.element(v)

    def find_generator(self, set_generator=False) -> 'FFElement':

        if self.has_valid_generator:
            g = self.generator
        elif maybe_prime(self.order - 1):
            # if the order of the multiplicative group is prime, any element but 1 and -1 is a generator
            p = self.base_polynomial.monic(self.dimension - 1)
            p += self.base_polynomial.monic(0)
            g = self.element_from_polynomial(p)
        else:
            # find it the hard way
            g = None
            for elt in iter(self):
                if elt != self.one and elt != self.zero:
                    o = self.element_order(elt)
                    if o == self.order - 1:
                        g = elt
                        break
            if g is None:
                raise ValueError(f'Unable to find a generator for {self}. Please check the field definition')
        if not self.has_valid_generator and set_generator:
            self.generator = g
            self._check_generator_order()
        return g

    def parse_poly(self, expr: str) -> Polynomial:
        """Returns a polynomial from its symbolic expression"""
        return symbolic_polynomial(expr, self)


class FFElement:
    def __init__(self, field: FiniteField, v: Vector):
        self.field = field
        self.vector = v

        assert isinstance(field, FiniteField)
        assert len(self.vector) == self.field.dimension

    @property
    def null(self) -> bool:
        return self == self.field.zero

    @property
    def is_scalar(self) -> bool:
        r = True
        for c in self.vector[1:]:
            if not c.is_zero:
                r = False
                break
        return r

    def __str__(self) -> str:
        return self.field.format_element(self)

    def __format__(self, format_spec: str) -> str:
        return self.field.format_element(self, format_spec)

    def __repr__(self) -> str:
        s = f'{self.__class__.__name__}('
        s += f'{repr(self.field)}, '
        s += str(self.vector)
        s += ')'
        return s

    def __getitem__(self, item: int) -> PFElement:
        return self.vector[item]

    def __add__(self, other: Any) -> 'FFElement':
        if isinstance(other, self.__class__):
            return self.field.add(self, other)
        elif isinstance(other, PFElement):
            return self.field.add(self, self.field(other))
        elif isinstance(other, int):
            return self.field.add(self, self.field(other))
        else:
            return operator.add(other, self)

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, int):
            other = self.field(self.field.prime_field(other))

        elif isinstance(other, PFElement):
            assert other.field == self.field.prime_field
            other = self.field(other)

        for i in range(0, self.field.dimension):
            if other[i] != self[i]:
                return False
        return True

    def __neg__(self) -> 'FFElement':
        return self.field.additive_inverse(self)

    def __sub__(self, other: Any) -> 'FFElement':
        return self.__add__(-other)

    def __mul__(self, other: Any) -> 'FFElement':
        if isinstance(other, FFElement):
            return self.field.mul(self, other)
        elif isinstance(other, int):
            return self.field.ext_mul(other, self)
        elif isinstance(other, PFElement):
            return self.field.mul(self, self.field(other))
        else:
            return operator.add(other, self)

    def __rmul__(self, other: Any) -> 'FFElement':
        return self.__mul__(other)

    def __pow__(self, e: int, modulo=None) -> 'FFElement':
        return self.field.pow(self, e)

    def __truediv__(self, other: 'FFElement') -> 'FFElement':
        return self.field.div(self, other)

    def __floordiv__(self, other: 'FFElement') -> 'FFElement':
        return self.__truediv__(other)

    def __mod__(self, other: Any) -> 'FFElement':
        return self.field.zero

    def __hash__(self) -> int:
        return hash(tuple(self.vector))

    def __len__(self) -> int:
        return len(self.vector)

    def __invert__(self) -> 'FFElement':
        return self.field.frobenius_reciprocal(self)


# Pre-computed finite fields
FINITE_FIELD = {
    'F2': PrimeField(2),
    'F3': PrimeField(3),
    'F4': FiniteField(2, 2, symbolic_polynomial('X^2+X+1', PrimeField(2)), generator=(0, 1)),
    'F5': PrimeField(5),
    'F7': PrimeField(7),
    'F8': FiniteField(2, 3, symbolic_polynomial('X^3+X+1', PrimeField(2)), generator=(0, 1, 0)),
    'F9': FiniteField(3, 2, symbolic_polynomial('X^2+1', PrimeField(3)), generator=(1, -1)),
    'F11': PrimeField(11),
    'F13': PrimeField(13),
    'F16': FiniteField(2, 4, symbolic_polynomial('X^4+X+1', PrimeField(2)), generator=(0, 1, 0, 0)),
    'F17': PrimeField(17),
    'F19': PrimeField(19),
    'F23': PrimeField(23),
    'F25': FiniteField(5, 2, symbolic_polynomial('X^2+2', PrimeField(5)), generator=(1, 1)),
    'F27': FiniteField(3, 3, symbolic_polynomial('X^3 - X - 1', PrimeField(3)), generator=(-1, -1, -1)),
}


def finite_field(order):
    if f'F{order}' in FINITE_FIELD.keys:
        return FINITE_FIELD[f'F{order}']
    else:
        raise ValueError(f'No pre-instantiated finite field of order {order}')
