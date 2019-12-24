from typing import Tuple, Iterator, Sequence, Union, Any


from pyimath.polynomial import Polynomial, symbolic_polynomial
from pyimath.integer import IntegerRing
from pyimath.functions import power


__all__ = ['GaussianIntegerRing', 'GaussInt']


class GaussianIntegerRing:
    """Set of numbers of the form `x + i*y` where `i^2 = -1`. Also defined as `Z[X]/<X^2+1>`
    """
    def __init__(self):
        self.root_symbol = 'i'
        self.base_polynomial = Polynomial([1, 0, 1], IntegerRing())  # X2 + 1 irreducible over Z

    def add(self, a: 'GaussInt', b: 'GaussInt') -> 'GaussInt':
        """Returns the sum of two gaussian integers
        """
        return self(a.x + b.x, a.y + b.y)

    def additive_inverse(self, a: 'GaussInt') -> 'GaussInt':
        """Returns the additive inverse of a gaussian integer
        """
        return self(-a.x, -a.y)

    @property
    def characteristic(self) -> int:
        """Returns the ring characteristic which is always zero.
        """
        return 0

    def divmod(self, a: 'GaussInt', b: 'GaussInt') -> Tuple['GaussInt', 'GaussInt']:
        """Returns the quotient and the remainder of the euclidean division of two gaussian integers
        """
        p_a = Polynomial([a.x, a.y], IntegerRing())
        p_b = Polynomial([b.x, b.y], IntegerRing())
        q, r = p_a.long_division_reversed(p_b)
        base = self.base_polynomial
        q %= base
        r %= base

        return self.element_from_polynomial(q), self.element_from_polynomial(r)

    @staticmethod
    def element(*args) -> 'GaussInt':
        """Returns a gaussian integer from a tuple of integers
        """
        return GaussInt(*args)

    def element_from_polynomial(self, p: Polynomial) -> 'GaussInt':
        """Returns the element `(a, b)` from the linear integer polynomial `a + bZ`
        """
        assert p.base_field == IntegerRing()
        assert p.degree < 2
        if p.degree == 0:
            return self(p[0], 0)
        else:
            return self(*tuple(p.coefficients))

    def ext_mul(self, n: int, a: 'GaussInt') -> 'GaussInt':
        """Returns the product of a gaussian integer and a regular integer
        """
        return self.mul(self(n, 0), a)

    def floor_div(self, a: 'GaussInt', b: 'GaussInt') -> 'GaussInt':
        """Returns the quotient of the euclidean division of two gaussian integers
        """
        q, _ = self.divmod(a, b)
        return q

    def format_element(self, e: 'GaussInt') -> str:
        """Returns a printable representation of a gaussian integer
        """
        return str(self.polynomial_from_element(e))

    def linear_polynomial(self, e: 'GaussInt') -> Polynomial:
        """Returns the polynomial `X - e`
        """
        return self.polynomial([e, self.one])

    def mod(self, a: 'GaussInt', b: 'GaussInt') -> 'GaussInt':
        """Returns the remainder of the division of two gaussian integers
        """
        _, r = self.divmod(a, b)
        return r

    def mul(self, a: 'GaussInt', b: 'GaussInt') -> 'GaussInt':
        """Returns the product of two gaussian integers
        """
        return self(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x)

    @property
    def neutral(self) -> 'GaussInt':
        """Same as `self.one`
        """
        return self.one

    @property
    def null(self) -> 'GaussInt':
        """Same as `self.zero`
        """
        return self.zero

    def polynomial(self, coefficients: Sequence, indeterminate='z') -> Polynomial:
        """Returns a polynomial from its coefficients
        """
        return Polynomial(coefficients, self, indeterminate)

    def polynomial_from_element(self, e: 'GaussInt') -> Polynomial:
        """Returns the polynomial `yZ + x` from the element `'(x, y)`
         """
        return Polynomial([e.x, e.y], IntegerRing(), indeterminate=self.root_symbol)

    @property
    def one(self) -> 'GaussInt':
        """Returns the multiplicative neutral
        """
        return GaussInt(1, 0)

    def parse(self, expr: str) -> 'GaussInt':
        """Returns a gaussian integer from its symbolic expression
        """
        p = symbolic_polynomial(expr, IntegerRing(), indeterminate=self.root_symbol)
        if p.degree > 1:
            raise ValueError(f'{expr} is not a valid gaussian integer')
        return self(p[0], p[1])

    @staticmethod
    def pow(a: 'GaussInt', n: int) -> 'GaussInt':
        """Returns the n-th power of a gaussian integer
        """
        return power(a, n)

    @property
    def zero(self) -> 'GaussInt':
        """Returns the additive neutral
        """
        return GaussInt(0, 0)

    def __call__(self, *args) -> 'GaussInt':
        """Returns a gausian integer from its arguments. Syntactic sugar for `self.element`
        """
        return self.element(*args)

    def __iter__(self) -> Iterator:
        """Kept for interface purposes
        """
        raise ValueError('Infinite sets are not iterable')

    def __repr__(self) -> str:
        """Returns an evaluable representation of the ring of gaussian integers
        """
        return f'{self.__class__.__name__}()'

    def __str__(self) -> str:
        """Returns a printable representation of the ring of gaussian integer
        """
        return "Field of Gauss's integers"


CompatibleType = Union[int, Tuple[int, int]]


class GaussInt:
    """Represents a Gaussian Integer
    """
    def __init__(self, x: int = 0, y: int = 0):
        self.vector = [x, y]
        self.ring = GaussianIntegerRing()

    @property
    def x(self) -> int:
        """Returns the integer component of a gaussian integer
        """
        return self.vector[0]

    @property
    def y(self) -> int:
        """Returns the complex component of a gaussian integer
        """
        return self.vector[1]

    @property
    def conjugate(self) -> 'GaussInt':
        """Returns the complex conjugate of a gaussian integer
        """
        return self.ring(self.x, -self.y)

    def __add__(self, other: Union['GaussInt', CompatibleType]) -> 'GaussInt':
        if isinstance(other, tuple):
            return self.ring.add(self, self.ring(*other))
        elif isinstance(other, int):
            return self.ring.add(self, self.ring(other, 0))
        else:
            return self.ring.add(self, other)

    def __radd__(self, other: Union['GaussInt', CompatibleType]) -> 'GaussInt':
        return self.__add__(other)

    def __sub__(self, other: Union['GaussInt', CompatibleType]) -> 'GaussInt':
        if isinstance(other, tuple):
            return self.ring.add(self, -self.ring(*other))
        elif isinstance(other, int):
            return self.ring.add(self, self.ring(other, 0))
        else:
            return self.ring.add(self, -other)

    def __rsub__(self, other: Union['GaussInt', CompatibleType]) -> 'GaussInt':
        if isinstance(other, tuple):
            return self.ring.add(self.ring(*other), -self)
        elif isinstance(other, int):
            return self.ring.add(self.ring(other, 0), -self)
        else:
            return self.ring.add(other, -self)

    def __neg__(self) -> 'GaussInt':
        return self.ring.additive_inverse(self)

    def __mul__(self, other: Union['GaussInt', CompatibleType]) -> 'GaussInt':
        if isinstance(other, tuple):
            return self.ring.mul(self, self.ring(*other))
        elif isinstance(other, int):
            return self.ring.mul(self, self.ring(other, 0))
        else:
            return self.ring.mul(self, other)

    def __rmul__(self, other: Union['GaussInt', CompatibleType]) -> 'GaussInt':
        return self.__mul__(other)

    def __pow__(self, e: int, modulo=None) -> 'GaussInt':
        return self.ring.pow(self, e)

    def __floordiv__(self, other: Union['GaussInt', CompatibleType]) -> 'GaussInt':
        if isinstance(other, tuple):
            return self.ring.floor_div(self, self.ring(*other))
        elif isinstance(other, int):
            return self.ring.floor_div(self, self.ring(other, 0))
        else:
            return self.ring.floor_div(self, other)

    def __truediv__(self, other: Any) -> None:
        raise TypeError('Gauss Integers do not define true division / but implement floor division //')

    def __mod__(self, other: Union['GaussInt', CompatibleType]) -> 'GaussInt':
        if isinstance(other, tuple):
            return self.ring.mod(self, self.ring(*other))
        elif isinstance(other, int):
            return self.ring.mod(self, self.ring(other, 0))
        else:
            return self.ring.mod(self, other)

    def __divmod__(self, other: Union['GaussInt', CompatibleType]) -> Tuple['GaussInt', 'GaussInt']:
        if isinstance(other, tuple):
            return self.ring.divmod(self, self.ring(*other))
        elif isinstance(other, int):
            return self.ring.divmod(self, self.ring(other, 0))
        else:
            return self.ring.divmod(self, other)

    def __eq__(self, other: Union['GaussInt', CompatibleType]) -> bool:
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y
        elif isinstance(other, tuple):
            return self.__eq__(self.ring(*other))
        elif isinstance(other, int):
            if self.y != 0:
                return False
            else:
                return other == self.x
        else:
            raise TypeError

    def __invert__(self) -> 'GaussInt':
        """Defines the complex conjugate"""
        return self.conjugate

    def __abs__(self) -> int:
        """Defines the multiplicative norm of algebraic integers"""
        return self.x ** 2 + self.y ** 2

    def __str__(self) -> str:
        return self.ring.format_element(self)

    def __repr__(self) -> str:
        s = f'{repr(self.ring)}('
        s += f'{self.x}, {self.y})'
        return s
