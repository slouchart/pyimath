from typing import Tuple

from imath.base import IntegralDomain
from imath.polynomial import Polynomial


class Z(IntegralDomain):
    """It's actually a ring"""

    def __init__(self):
        super().__init__(0)

    def element(self, n: int) -> int:
        return int(n)

    @property
    def zero(self) -> int:
        return 0

    @property
    def one(self) -> int:
        return 1

    def add(self, a: int, b: int) -> int:
        return a + b

    def additive_inverse(self, a: int) -> int:
        return -a

    def mul(self, a: int, b: int) -> int:
        return a * b

    def ext_mul(self, n: int, a: int) -> int:
        return n * a

    def floor_div(self, a: int, b: int) -> int:
        return a // b

    def mod(self, a: int, b: int) -> int:
        return a % b

    def divmod(self, a: int, b: int) -> Tuple[int, int]:
        return divmod(a, b)

    def pow(self, a: int, e: int) -> int:
        return a ** e

    def __eq__(self, other):
        return isinstance(other, self.__class__)

    def __call__(self, n: int) -> int:
        return self.element(n)

    def __iter__(self):
        raise NotImplementedError("Iteration over infinite number sets is undefined")

    def __str__(self) -> str:
        return f'Pseudo-field (i.e ring) of integers'

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}()'

    def polynomial(self, *args, indeterminate='X'):
        coefficients = [c for c in args]
        return Polynomial(coefficients, base_field=self, indeterminate=indeterminate)

    def linear_polynomial(self, e: int):
        return self.polynomial(-e, 1)

    def p_th_root(self, a):
        raise NotImplementedError
