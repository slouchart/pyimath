from itertools import product as cartesian_product
import operator

from typing import Optional, Sequence, MutableSequence, Iterator, Union, List, Tuple, Any, Collection


from pyimath.functions import maybe_prime
from pyimath.functions import power
from pyimath.polynomial import Polynomial, symbolic_polynomial
from pyimath.primefield import PrimeField, PFElement

__all__ = ['FiniteField', 'finite_field', 'FFElement']


BaseAtom = Union[int, PFElement]
Vector = Union[List[BaseAtom], MutableSequence[BaseAtom], Sequence[BaseAtom], Collection[BaseAtom]]


class FiniteField:
    """Represents a non-prime finite field where:
    """
    def __init__(self,
                 prime: int,
                 dimension: int,
                 ideal: Polynomial,
                 generator: Optional[Sequence[int]] = None,
                 root_symbol: Optional[str] = 'j'):
        """
        `prime` is the order of the base prime field,

        `dimension` is the dimension of the finite field seen as a vector space over its prime field,

        `ideal` is the minimal polynomial of the adjunct root j, must be irreducible over the prime field,

        `generator` must be an element of the field that generates the multiplicative group

        and `root_symbol` is a one-character string to represent the adjunct root
        """
        self._prime_field = PrimeField(prime)
        self.dimension = dimension
        self.base_polynomial = ideal
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

    def add(self, a: 'FFElement', b: 'FFElement') -> 'FFElement':
        """Adds two elements
        """
        assert a.field == b.field == self
        s = [self.prime_field.zero] * self.dimension
        for i in range(0, self.dimension):
            s[i] = a[i] + b[i]
        return self.element(s)

    def additive_inverse(self, a: 'FFElement') -> 'FFElement':
        """Returns the additive inverse of an element
        """
        return self.element([-e for e in a.vector])

    @property
    def basis(self) -> Iterator['FFElement']:
        """Returns a vector-space basis of field elements
        """
        return (self(0, 1) ** e for e in range(0, self.dimension))

    @property
    def characteristic(self) -> int:
        """Returns the characteristic of the field
        """
        return self._prime_field.characteristic

    def check_irreducible_polynomial(self, p: Polynomial) -> bool:
        """Returns `True` if the polynomial is irreducible, `False` otherwise
        """
        assert p.base_field == self
        return p.is_irreducible

    def div(self, a: 'FFElement', b: 'FFElement') -> 'FFElement':
        """Returns the true division of two elements
        """
        if b.null:
            raise ZeroDivisionError
        if a.null:
            return self.zero

        return a * self.multiplicative_inverse(b)

    def divmod(self, a: 'FFElement', b: 'FFElement') -> Tuple['FFElement', 'FFElement']:
        """Kept for interface purpose, returns always `(a / b, self.zero)`
        """
        return self.div(a, b), self.mod(a, b)

    def element(self, v: Union[Vector, 'FFElement']) -> 'FFElement':
        """Returns an instance of a `FFElement` from a vector
        """
        return FFElement(self, self._safe_convert_vector(v))

    def element_from_polynomial(self, p: Polynomial) -> 'FFElement':
        """Return an element whose basis components are the coefficients of a given polynomial over the field
        e.g. The polynomial `1 + X + X^2` would yield `(1+j+j^2)`
        """
        assert p.base_field == self.prime_field
        assert p.degree <= self.dimension - 1
        v = p.coefficients
        if p.degree < self.dimension - 1:
            v += [self.prime_field.zero] * (self.dimension - (1 + p.degree))
        return self.element(v)

    def element_order(self, e: 'FFElement') -> int:
        """Returns the order of an element, that is the minimal integer k such as e^k = 1
        """
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

    def ext_mul(self, n: int, a: 'FFElement') -> 'FFElement':
        """Returns the n-th iterated addition of an element with itself
        """
        assert a.field == self
        r = self.zero
        for _ in range(0, n):
            r += a
        return r

    def find_generator(self, set_generator=False) -> 'FFElement':
        """Returns a valid generator of the multiplicative group of the field
        """
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

    def format_element(self, e: 'FFElement', format_spec: str = '') -> str:
        """Returns a printable representation of an element
        """
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

    @property
    def frobenius_map(self) -> MutableSequence['FFElement']:
        """Returns what is actually the linear map of the inverse function of Frobenius automorphism
        (see MISCELLANEOUS.md)
        """
        return self._frobenius_map[:]

    def frobenius_reciprocal(self, a: 'FFElement') -> 'FFElement':
        """Returns the bijective inverse of an element through the Frobenius automorphism. This is
        equivalent of taking the p-th root of an element where p is the field characteristic
        """
        assert isinstance(a, FFElement)
        assert a in self
        r = self.zero
        for i, e in enumerate(a.vector):
            b = self._frobenius_map[i]
            r += self(e) * b
        return r

    def linear_polynomial(self, e: 'FFElement') -> Polynomial:
        """Returns the polynomial `X - e`
        """
        p = self.polynomial(-e)
        p += p.monic(1)
        return p

    @property
    def has_valid_generator(self) -> bool:
        """Returns `True` if the field was set with a valid generator i.e. an element whose order is equal
        to the order of the field
        """
        return len(self.generator_powers) > 0 and self.generator is not None

    def mod(self, _, b: 'FFElement') -> 'FFElement':
        """Kept for interface purpose, returns always `self.zero`
        """
        if b.null:
            raise ZeroDivisionError

        return self.zero

    def mul(self, a: 'FFElement', b: 'FFElement') -> 'FFElement':
        """Returns the product of two elements
        """
        assert a.field == b.field == self

        if self.has_valid_generator:

            if a.null or b.null:
                return self.null

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

    def multiplicative_inverse(self, a: 'FFElement') -> 'FFElement':
        """Returns the multiplicative inverse of an element
        """
        if self.has_valid_generator:
            p = self.element_as_powers[a]
            if self.order - 1 - p == 0:
                return self.generator_powers[self.order - 1]
            else:
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

    @property
    def neutral(self) -> 'FFElement':
        """Same as `self.one`, returns the multiplicative neutral element of the field
        """
        return self.one

    @property
    def null(self) -> 'FFElement':
        """Same as `self.zero`, returns the additive neutral element of the field
        """
        return self.zero

    def parse_poly(self, expr: str) -> Polynomial:
        """Returns a polynomial from its symbolic expression
        """
        return symbolic_polynomial(expr, self)

    def polynomial(self, *args, indeterminate: Optional[str] = 'X') -> Polynomial:
        """Returns a polynomial from a list of coefficients. The coefficients are converted to `FFElement`
        """
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

    def polynomial_from_element(self, e: 'FFElement') -> Polynomial:
        """Returns a polynomial whose coefficients are those of the element
        e.g. `(1+j+j^2)` would yield `1 + X + X^2`
        """
        return Polynomial(e.vector, self.prime_field, indeterminate=self.root_symbol)

    @property
    def one(self) -> 'FFElement':
        """Returns the multiplicative neutral element of the field
        """
        return self.element([self.prime_field.one] + [self.prime_field.zero] * (self.dimension - 1))

    @property
    def order(self) -> int:
        """Returns the order of the field
        """
        return self.prime_field.characteristic ** self.dimension

    def pow(self, a: 'FFElement', n: int) -> 'FFElement':
        """Returns the n-th power of an element
        """
        assert a.field == self
        assert n >= 0
        res = power(a, n)
        if not isinstance(res, FFElement):
            return self(res)
        else:
            return res

    @property
    def prime_field(self) -> PrimeField:
        """Returns the instance of the base prime field
        """
        return self._prime_field

    @property
    def zero(self) -> 'FFElement':
        """Returns the additive neutral element of the field
        """
        return self.element([self.prime_field.zero] * self.dimension)

    def __call__(self, *args) -> 'FFElement':
        """Creates an element from a list of arguments, syntactic sugar for `self.element()`
        """
        return self.element([self.prime_field(c) for c in args])

    def __eq__(self, other: 'FiniteField') -> bool:
        """Compares two finite fields
        """
        assert isinstance(other, self.__class__)
        return other.prime_field == self.prime_field and other.dimension == self.dimension

    def __iter__(self) -> Iterator['FFElement']:
        """Iterates over all the elements
        """
        iterables = [self.prime_field] * self.dimension
        return iter((self.element(list(v)) for v in cartesian_product(*iterables)))

    def __repr__(self) -> str:
        """Returns an evaluable representation of the field
        """
        s = f'{self.__class__.__name__}('
        s += f'{self.prime_field.characteristic}, '
        s += f'{self.dimension}, '
        s += f'{repr(self.base_polynomial)}, '
        if self.has_valid_generator:
            s += f' generator={repr(self.generator.vector)})'
        else:
            s += ')'
        return s

    def __str__(self) -> str:
        """Returns a printable representation of the finite field
        """
        return f'Finite field of order {self.prime_field.characteristic}^{self.dimension}'

    def _check_generator_order(self):
        """Computes the powers of the wanna-be generator and checks
        if the results form the multiplicative group. This method can be
        computer-intensive for large finite fields
        """
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

    def _compute_frobenius_map(self) -> MutableSequence['FFElement']:
        """The Frobenius automorphisme is defined by a -> a^p where p is the prime field characteristic.
        This methods computes the inverse map of the automorphism
        """
        p_th_roots = []
        for b in self.basis:
            for e in self:
                if e ** self.characteristic == b:
                    p_th_roots.append(e)
        return p_th_roots

    def _compute_root_powers(self):
        """Creates a list of elements j^p where j is the adjunct root and p is an integer
        in the range `self.dimension..self.order-1`
        """
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
        """Creates an element of the field from a variety of input values
        """
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


class FFElement:
    def __init__(self, field: FiniteField, v: Vector):
        """Represents an element of a finite field of non-prime order as a vector
        of elements of the base prime field
        """
        self.field = field
        self.vector = v

        assert isinstance(field, FiniteField)
        assert len(self.vector) == self.field.dimension

    @property
    def is_scalar(self) -> bool:
        """Returns `True` if the element belongs to the base prime field, `False` otherwise
        """
        r = True
        for c in self.vector[1:]:
            if not c.is_zero:
                r = False
                break
        return r

    @property
    def null(self) -> bool:
        """Returns `True` if the element is the additive neutral of the field
        """
        return self == self.field.zero

    def __add__(self, other: Any) -> 'FFElement':
        """Returns the sum of two elements
        """
        if isinstance(other, self.__class__):
            return self.field.add(self, other)
        elif isinstance(other, PFElement):
            return self.field.add(self, self.field(other))
        elif isinstance(other, int):
            return self.field.add(self, self.field(other))
        else:
            return operator.add(other, self)

    def __eq__(self, other: Any) -> bool:
        """Returns `True` if two elements are equal, `False` otherwise
        """
        if isinstance(other, int):
            other = self.field(self.field.prime_field(other))

        elif isinstance(other, PFElement):
            assert other.field == self.field.prime_field
            other = self.field(other)

        for i in range(0, self.field.dimension):
            if other[i] != self[i]:
                return False
        return True

    def __floordiv__(self, other: 'FFElement') -> 'FFElement':
        """Support for the `//` operator. Same as `self.__truediv__`
        """
        return self.__truediv__(other)

    def __format__(self, format_spec: str) -> str:
        """Returns a customizable printable representation of an element
        """
        return self.field.format_element(self, format_spec)

    def __getitem__(self, n: int) -> PFElement:
        """Returns the n-th basis component of an element
        """
        return self.vector[n]

    def __hash__(self) -> int:
        """Returns a hash for map `dict` purposes
        """
        return hash(tuple(self.vector))

    def __invert__(self) -> 'FFElement':
        """Support for the `~` operator. Same as `self.field.frobenius_reciprocal`
        """
        return self.field.frobenius_reciprocal(self)

    def __len__(self) -> int:
        """Returns the length of the element, that is the dimension of the field basis
        """
        return len(self.vector)

    def __mod__(self, other: Any) -> 'FFElement':
        """Returns always `self.field.zero`
        """
        return self.field.zero

    def __mul__(self, other: Any) -> 'FFElement':
        """Returns the product of two elements
        """
        if isinstance(other, FFElement):
            return self.field.mul(self, other)
        elif isinstance(other, int):
            return self.field.ext_mul(other, self)
        elif isinstance(other, PFElement):
            return self.field.mul(self, self.field(other))
        else:
            return operator.add(other, self)

    def __neg__(self) -> 'FFElement':
        """Returns the additive inverse of an element
        """
        return self.field.additive_inverse(self)

    def __pow__(self, n: int, modulo=None) -> 'FFElement':
        """Returns the n-th power of an element
        """
        return self.field.pow(self, n)

    def __repr__(self) -> str:
        """Returns an evaluable representation of an element
        """
        s = f'{self.__class__.__name__}('
        s += f'{repr(self.field)}, '
        s += str(self.vector)
        s += ')'
        return s

    def __rmul__(self, other: Any) -> 'FFElement':
        """Same as `self.__mul__`, used when adding non-`FFElement` instances to a `FFElement`
        """
        return self.__mul__(other)

    def __str__(self) -> str:
        """Returns a printable representation of an element
        """
        return self.field.format_element(self)

    def __sub__(self, other: Any) -> 'FFElement':
        """Returns the subtraction of two elements
        """
        return self.__add__(-other)

    def __truediv__(self, other: 'FFElement') -> 'FFElement':
        """Returns the quotient of two elements
        """
        if isinstance(other, FFElement):
            return self.field.div(self, other)
        elif isinstance(other, (int, PFElement)):
            return self.field.mul(self, self.field.multiplicative_inverse(self.field(other)))


# Pre-computed finite fields
class Record:
    """Record in the register of pre-instantiated finite field
    """
    def __init__(self, code: str):
        self.code = code
        self.instance = None

    def instantiate(self):
        if not self.instantiated:
            self.instance = eval(self.code)

        return self.instance

    @property
    def instantiated(self) -> bool:
        return self.instance is not None


FINITE_FIELD = {
    'F2': Record("PrimeField(2)"),
    'F3': Record("PrimeField(3)"),
    'F4': Record("FiniteField(2, 2, symbolic_polynomial('X^2+X+1', PrimeField(2)), generator=(0, 1))"),
    'F5': Record("PrimeField(5)"),
    'F7': Record("PrimeField(7)"),
    'F8': Record("FiniteField(2, 3, symbolic_polynomial('X^3+X+1', PrimeField(2)), generator=(0, 1))"),
    'F9': Record("FiniteField(3, 2, symbolic_polynomial('X^2+1', PrimeField(3)), generator=(1, 1))"),
    'F11': Record("PrimeField(11)"),
    'F13': Record("PrimeField(13)"),
    'F16': Record("FiniteField(2, 4, symbolic_polynomial('X^4+X+1', PrimeField(2)), generator=(1, 1, 0, 0))"),
    'F17': Record("PrimeField(17)"),
    'F19': Record("PrimeField(19)"),
    'F23': Record("PrimeField(23)"),
    'F25': Record("FiniteField(5, 2, symbolic_polynomial('X^2+2', PrimeField(5)), generator=(1, 1))"),
    'F27': Record("FiniteField(3, 3, symbolic_polynomial('X^3-X-1', PrimeField(3)), generator=(1, 0, 1))")
}


def finite_field(order):
    """Returns the finite field of the given order.
    Fields that may be instantiated this way are:

    * all prime fields of order less then 23
    * finite fields of non-prime order up to order 27. These finite fields have all a valid generator.

    """
    if f'F{order}' in FINITE_FIELD:
        return FINITE_FIELD[f'F{order}'].instantiate()
    else:
        raise ValueError(f'No pre-instantiated finite field of order {order}')
