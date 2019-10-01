from imath.polynomial import Polynomial
from imath.polyparse import symbolic_polynomial
from imath.integer import IntegerRing
from imath.functions import power


__all__ = ['GaussianIntegerRing']


class GaussianIntegerRing:
    def __init__(self):
        self.root_symbol = 'i'
        self.base_polynomial = Polynomial([1, 0, 1], IntegerRing())  # X2 + 1 irreducible over Z

    @property
    def characteristic(self):
        return 0

    @staticmethod
    def element(*args):
        return GaussInt(*args)

    @property
    def zero(self):
        return GaussInt(0, 0)

    @property
    def one(self):
        return GaussInt(1, 0)

    @property
    def null(self):
        return self.zero

    @property
    def neutral(self):
        return self.one

    def add(self, a, b):
        return self(a.x + b.x, a.y + b.y)

    def additive_inverse(self, a):
        return self(-a.x, -a.y)

    def mul(self, a, b):
        return self(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x)

    def ext_mul(self, n: int, a):
        return self.mul(self(n, 0), a)

    def divmod(self, a, b):
        p_a = Polynomial([a.x, a.y], IntegerRing())
        p_b = Polynomial([b.x, b.y], IntegerRing())
        q, r = p_a.long_division_reversed(p_b)
        base = self.base_polynomial
        q %= base
        r %= base

        return self.element_from_polynomial(q), self.element_from_polynomial(r)

    def floor_div(self, a, b):
        q, _ = self.divmod(a, b)
        return q

    def mod(self, a, b):
        _, r = self.divmod(a, b)
        return r

    @staticmethod
    def pow(a, e: int):
        return power(a, e)

    def __call__(self, *args):
        return self.element(*args)

    def __iter__(self):
        raise ValueError('Infinite sets are not iterable')

    def __str__(self) -> str:
        return "Field of Gauss's integers"

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}()'

    def format_element(self, e) -> str:
        return str(self.polynomial_from_element(e))

    def polynomial(self, coefficients, indeterminate='z') -> Polynomial:
        return Polynomial(coefficients, self, indeterminate)

    def linear_polynomial(self, e) -> Polynomial:
        return self.polynomial([e, self.one])

    def polynomial_from_element(self, e) -> Polynomial:
        return Polynomial([e.x, e.y], IntegerRing(), indeterminate=self.root_symbol)

    def element_from_polynomial(self, p: Polynomial):
        assert p.base_field == IntegerRing()
        assert p.degree < 2
        if p.degree == 0:
            return self(p.coefficients[0], 0)
        else:
            return self(*tuple(p.coefficients))

    def parse(self, expr):
        p = symbolic_polynomial(expr, IntegerRing(), indeterminate=self.root_symbol)
        if p.degree > 1:
            raise ValueError(f'{expr} is not a valid gaussian integer')
        return self(p[0], p[1])


class GaussInt:
    def __init__(self, x: int = 0, y: int = 0):
        self.vector = [x, y]
        self.ring = GaussianIntegerRing()

    @property
    def x(self):
        return self.vector[0]

    @property
    def y(self):
        return self.vector[1]

    @property
    def conjugate(self):
        return self.ring(self.x, -self.y)

    def __add__(self, other):
        if isinstance(other, tuple):
            return self.ring.add(self, self.ring(*other))
        elif isinstance(other, int):
            return self.ring.add(self, self.ring(other, 0))
        else:
            return self.ring.add(self, other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, tuple):
            return self.ring.add(self, -self.ring(*other))
        elif isinstance(other, int):
            return self.ring.add(self, self.ring(other, 0))
        else:
            return self.ring.add(self, -other)

    def __rsub__(self, other):
        if isinstance(other, tuple):
            return self.ring.add(self.ring(*other), -self)
        elif isinstance(other, int):
            return self.ring.add(self.ring(other, 0), -self)
        else:
            return self.ring.add(other, -self)

    def __neg__(self):
        return self.ring.additive_inverse(self)

    def __mul__(self, other):
        if isinstance(other, tuple):
            return self.ring.mul(self, self.ring(*other))
        elif isinstance(other, int):
            return self.ring.mul(self, self.ring(other, 0))
        else:
            return self.ring.mul(self, other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, power, modulo=None):
        return self.ring.pow(self, power)

    def __floordiv__(self, other):
        if isinstance(other, tuple):
            return self.ring.floor_div(self, self.ring(*other))
        elif isinstance(other, int):
            return self.ring.floor_div(self, self.ring(other, 0))
        else:
            return self.ring.floor_div(self, other)

    def __truediv__(self, other):
        raise TypeError('Gauss Integers do not define true division / but implement floor division //')

    def __mod__(self, other):
        if isinstance(other, tuple):
            return self.ring.mod(self, self.ring(*other))
        elif isinstance(other, int):
            return self.ring.mod(self, self.ring(other, 0))
        else:
            return self.ring.mod(self, other)

    def __divmod__(self, other):
        if isinstance(other, tuple):
            return self.ring.divmod(self, self.ring(*other))
        elif isinstance(other, int):
            return self.ring.divmod(self, self.ring(other, 0))
        else:
            return self.ring.divmod(self, other)

    def __eq__(self, other) -> bool:
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

    def __invert__(self):
        """Defined the complex conjugate"""
        return self.conjugate

    def __abs__(self):
        """Defines the multiplicative norm of algebraic integers"""
        return self.x ** 2 + self.y ** 2

    def __str__(self):
        return self.ring.format_element(self)

    def __repr__(self):
        s = f'{repr(self.ring)}('
        s += f'{self.x}, {self.y})'
        return s
