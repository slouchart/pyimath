

class IntegralDomain:
    """Defines an abstract class for an integral domain that is:
       a set which is an abelian group for + and a semi-group
       for * and where GCD are computable.
       A ring of polynomials is at least defined over an integral domain
       examples: Z, Z(i)"""

    def __init__(self, n: int):
        assert n >= 0
        self._characteristic = n

    @property
    def characteristic(self):
        return self._characteristic

    def element(self, n):
        raise NotImplementedError

    @property
    def zero(self):
        raise NotImplementedError

    @property
    def null(self):
        return self.zero

    @property
    def one(self):
        raise NotImplementedError

    @property
    def neutral(self):
        return self.one

    def add(self, a, b):
        raise NotImplementedError

    def additive_inverse(self, a):
        raise NotImplementedError

    def mul(self, a, b):
        raise NotImplementedError

    def ext_mul(self, n: int, a):
        raise NotImplementedError

    def floor_div(self, a, b):
        raise NotImplementedError

    def mod(self, a, b):
        raise NotImplementedError

    def divmod(self, a, b):
        raise NotImplementedError

    def pow(self, a, n: int):
        raise NotImplementedError

    def __eq__(self, other) -> bool:
        raise NotImplementedError

    def __call__(self, n):
        raise NotImplementedError

    def __str__(self) -> str:
        raise NotImplementedError

    def __repr__(self) -> str:
        raise NotImplementedError

    def polynomial(self, coefficients, indeterminate='X'):
        raise NotImplementedError

    def linear_polynomial(self, a):
        raise NotImplementedError

    def p_th_root(self, a):
        raise NotImplementedError
