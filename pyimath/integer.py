from typing import Tuple, Iterator, Optional
import operator

from pyimath.polynomial import Polynomial

__all__ = ['IntegerRing']


class IntegerRing:
    """
    Mimics a Field as defined in `primefield` or `finitefield` modules (and only mimics a field, Z is a ring!)

    Mainly used to define polynomials from `Z[X]`
    """

    def __init__(self):
        ...

    @staticmethod
    def add(a: int, b: int) -> int:
        """Returns the sum of two integers"""
        return operator.add(a, b)

    @staticmethod
    def additive_inverse(a: int) -> int:
        """Returns the additive inverse of an integer"""
        return operator.neg(a)

    @property
    def characteristic(self) -> int:
        """Returns the characteristic of Z which is always 0 by definition"""
        return 0

    @staticmethod
    def element(n: int) -> int:
        """Element constructor. Kept for interface compatibility"""
        return int(n)

    @staticmethod
    def ext_mul(n: int, a: int) -> int:
        """Kept for interface purposes. Simply returns `n*a`"""
        return operator.mul(n, a)

    @staticmethod
    def divmod(a: int, b: int) -> Tuple[int, int]:
        """Returns the quotient and the remainder of the integer division"""
        return divmod(a, b)

    @staticmethod
    def floor_div(a: int, b: int) -> int:
        """Return the quotient of the integer division"""
        return operator.floordiv(a, b)

    def linear_polynomial(self, e: int) -> Polynomial:
        """Returns a polynomial defined as `X - e`"""
        return self.polynomial(-e, 1)

    @staticmethod
    def mod(a: int, b: int) -> int:
        """Returns the remainder of the integer division"""
        return operator.mod(a, b)

    @staticmethod
    def mul(a: int, b: int) -> int:
        """Return the product of two integers"""
        return operator.mul(a, b)

    @property
    def neutral(self) -> int:
        """Same as `self.one`, returns the multiplicative neutral element of Z, 1"""
        return self.one

    @property
    def null(self) -> int:
        """Same as `self.zero`, returns the additive neutral element of Z, 0"""
        return self.zero

    @property
    def one(self) -> int:
        """Returns the multiplicative neutral element of Z"""
        return 1

    def polynomial(self, *args, indeterminate: Optional[str] = 'X') -> Polynomial:
        """Returns a polynomial built from a list of coefficients."""
        coefficients = [c for c in args]
        return Polynomial(coefficients, base_field=self, indeterminate=indeterminate)

    @staticmethod
    def pow(a: int, e: int) -> int:
        """Returns `a` to the `e`-th power """
        return a ** e

    @property
    def zero(self) -> int:
        """Returns the additive neutral element of Z"""
        return 0

    def __call__(self, n: int) -> int:
        return self.element(n)

    def __eq__(self, other: 'IntegerRing') -> bool:
        return isinstance(other, self.__class__)

    def __iter__(self) -> Iterator:
        raise NotImplementedError("Iteration over infinite number sets is undefined")

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}()'

    def __str__(self) -> str:
        return f'Ring of integers'
