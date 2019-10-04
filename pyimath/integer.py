from typing import Tuple, Iterator

from pyimath.polynomial import Polynomial

__all__ = ['IntegerRing']


class IntegerRing:
    """
    Mimic a Field as defined in primefield or finitefield modules
    used to define polynomials from Z[X]
    """

    def __init__(self):
        ...

    @staticmethod
    def element(n: int) -> int:
        return int(n)

    @property
    def characteristic(self):
        return 0

    @property
    def zero(self) -> int:
        return 0

    @property
    def one(self) -> int:
        return 1

    @property
    def null(self):
        return self.zero

    @property
    def neutral(self):
        return self.one

    @staticmethod
    def add(a: int, b: int) -> int:
        return a + b

    @staticmethod
    def additive_inverse(a: int) -> int:
        return -a

    @staticmethod
    def mul(a: int, b: int) -> int:
        return a * b

    @staticmethod
    def ext_mul(n: int, a: int) -> int:
        return n * a

    @staticmethod
    def floor_div(a: int, b: int) -> int:
        return a // b

    @staticmethod
    def mod(a: int, b: int) -> int:
        return a % b

    @staticmethod
    def divmod(a: int, b: int) -> Tuple[int, int]:
        return divmod(a, b)

    @staticmethod
    def pow(a: int, e: int) -> int:
        return a ** e

    def __eq__(self, other: 'IntegerRing') -> bool:
        return isinstance(other, self.__class__)

    def __call__(self, n: int) -> int:
        return self.element(n)

    def __iter__(self) -> Iterator:
        raise NotImplementedError("Iteration over infinite number sets is undefined")

    def __str__(self) -> str:
        return f'Ring of integers'

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}()'

    def polynomial(self, *args, indeterminate='X'):
        coefficients = [c for c in args]
        return Polynomial(coefficients, base_field=self, indeterminate=indeterminate)

    def linear_polynomial(self, e: int):
        return self.polynomial(-e, 1)
