from typing import Tuple, Iterator, Optional
import operator

from pyimath.polynomial import Polynomial

__all__ = ['IntegerRing']


class IntegerRing:
    """
    Mimic a Field as defined in `primefield` or `finitefield` modules
    used to define polynomials from `Z[X]`
    """

    def __init__(self):
        ...

    @staticmethod
    def add(a: int, b: int) -> int:
        return operator.add(a, b)

    @staticmethod
    def additive_inverse(a: int) -> int:
        return operator.neg(a)

    @property
    def characteristic(self) -> int:
        return 0

    @staticmethod
    def element(n: int) -> int:
        return int(n)

    @staticmethod
    def ext_mul(n: int, a: int) -> int:
        return operator.mul(n, a)

    @staticmethod
    def divmod(a: int, b: int) -> Tuple[int, int]:
        return divmod(a, b)

    @staticmethod
    def floor_div(a: int, b: int) -> int:
        return operator.floordiv(a, b)

    def linear_polynomial(self, e: int) -> Polynomial:
        return self.polynomial(-e, 1)

    @staticmethod
    def mod(a: int, b: int) -> int:
        return operator.mod(a, b)

    @staticmethod
    def mul(a: int, b: int) -> int:
        return operator.mul(a, b)

    @property
    def neutral(self) -> int:
        return self.one

    @property
    def null(self) -> int:
        return self.zero

    @property
    def one(self) -> int:
        return 1

    def polynomial(self, *args, indeterminate: Optional[str] = 'X') -> Polynomial:
        coefficients = [c for c in args]
        return Polynomial(coefficients, base_field=self, indeterminate=indeterminate)

    @staticmethod
    def pow(a: int, e: int) -> int:
        return a ** e

    @property
    def zero(self) -> int:
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
