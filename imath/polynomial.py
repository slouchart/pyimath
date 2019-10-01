from imath.functions import gcd, rgcd, power
from imath.polyparse import symbolic_polynomial


__all__ = ['Polynomial']


class Polynomial:
    """Represents a polynomial over a finite field or over the integers
    More generally, definition of polynomials over a ring is possible
    but not recommended
    The internal representation of the coefficients uses a dict that
    indexes the coefficients by their degree
    The usual definition of a polynomial as a sequence of numbers
    that are all zeroes from a certain index is used to initialize the polynomial
    on instantiation and can be retrieved through the coefficients property
    Moreover this class implements all the operations over a ring of polynomials

    base_field must represent an integral domain that is:
       a set which is an abelian group for + and a semi-group
       for * and where GCD are computable.
       A ring of polynomials is at least defined over an integral domain
       examples: Z, Z(i)
    """

    def __init__(self, coeffs, base_field, indeterminate='X'):
        """
        :param coeffs: an iterable of elements from the base field
        :param base_field: a NumberSet instance generally a finite field
        :param indeterminate: a single digit string used to format the polynomial
        """
        self.base_field = base_field
        assert hasattr(base_field, 'element') and hasattr(base_field, '__call__')

        self._coefficients = self._safe_convert_coefficients(self._remove_trailing_zeros(coeffs))
        self.indeterminate = indeterminate

    def _safe_convert_coefficients(self, seq):
        bf = self.base_field
        return dict({deg: bf.element(c) for deg, c in enumerate(seq) if c != bf.zero})

    def _remove_trailing_zeros(self, seq):
        if len(seq) == 0:
            return []
        revseq = list(reversed(seq))
        while revseq[0] == self.base_field.zero:
            revseq = revseq[1:]
            if len(revseq) == 0:
                return []

        return list(reversed(revseq))

    def set_term(self, deg, c):
        if c == self.base_field.zero:
            if deg in self._coefficients.keys():
                del self._coefficients[deg]
        else:
            self._coefficients[deg] = c

    @property
    def coefficients(self):
        deg = 0
        res = []
        while deg <= self.degree:
            res.append(self[deg])
            deg += 1
        return res

    @property
    def internal(self):
        return dict({deg: c for deg, c in self._coefficients.items()})

    @property
    def copy(self):
        res = self.null
        for deg, c in self.internal.items():
            res.set_term(deg, c)
        return res

    @property
    def null(self):
        return Polynomial([], base_field=self.base_field, indeterminate=self.indeterminate)

    def monic(self, degree: int = 1):
        res = self.null
        res.set_term(degree, self.base_field.one)
        return res

    @property
    def unit(self):
        return self.monic(0)

    @property
    def is_constant(self) -> bool:
        return self.degree == 0 and not self.is_null

    @property
    def is_monic(self) -> bool:
        return self._coefficients[self.degree] == self.base_field.one

    @property
    def is_null(self) -> bool:
        return len(self._coefficients.keys()) == 0

    @property
    def is_unit(self) -> bool:
        return self == self.monic(0)

    @property
    def is_abs_unit(self) -> bool:
        return self in (self.monic(0), -self.monic(0))

    @property
    def is_irreducible(self) -> bool:
        return self.check_irreducibility()

    @property
    def degree(self) -> int:
        if self.is_null:
            return 0  # rigorously, it should be -infinity
        else:
            return max(self._coefficients.keys())

    @property
    def constant(self):
        if not self.is_null:
            return self[0]
        else:
            return self.base_field.zero

    @property
    def leading(self):
        if not self.is_null:
            return self._coefficients[self.degree]
        else:
            return self.base_field.zero

    @property
    def trailing(self):
        if self.is_null:
            return self.base_field.zero
        else:
            return self._coefficients[self.valuation]

    @property
    def valuation(self) -> int:

        if self.is_null:
            raise ValueError('The valuation of the null polynomial is undefined')

        return min(self._coefficients.keys())

    def make_monic(self):
        """
        Attempts to divide the polynomial by its leading coefficient (for a field) or by the gcd of its
        coefficients (for a ring)
        :return: a monic polynomial or raises an error if the polynomial cannot be made monic
        """
        if not self.is_monic:
            if self.base_field.characteristic != 0:
                return self / self.leading
            else:
                g = rgcd(self.coefficients)
                if self.leading // g == self.base_field.one:
                    return self // self.leading
                else:
                    raise ValueError(f'Polynomial {self} over {self.base_field} cannot be made monic')

        else:
            return self.copy

    def __len__(self) -> int:
        return len([c for c in self._coefficients.values() if c != 0])

    def __getitem__(self, item):
        if item in self._coefficients.keys():
            return self._coefficients[item]
        else:
            return self.base_field.zero

    def __eq__(self, other) -> bool:
        if isinstance(other, Polynomial):
            if self.degree == other.degree:

                if self.is_null or other.is_null:
                    return self.is_null and other.is_null

                a = self.coefficients
                b = other.coefficients

                assert len(a) == len(b)
                return a == b
            else:
                return False
        else:
            other = self(other)
            return self == other

    def add(self, poly):

        a = self._coefficients
        s = poly.copy
        for deg_a, c_a in a.items():
            if deg_a in s.internal.keys():
                s.set_term(deg_a, s[deg_a] + c_a)
            else:
                s.set_term(deg_a, c_a)

        return s

    def add_constant(self, k):
        a = self._coefficients
        s = self.null

        if not self.is_null:
            for deg, c in a.items():
                if deg == 0:
                    s.set_term(deg, c + k)
                else:
                    s.set_term(deg, c)
        else:
            s.set_term(0, k)

        return s

    def __add__(self, other):
        if isinstance(other, Polynomial):
            return self.add(other)
        elif isinstance(other, list):
            return self.add(Polynomial(other, base_field=self.base_field))
        else:
            return self.add_constant(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        a = self._coefficients
        s = self.null
        for deg, c in a.items():
            s.set_term(deg, -c)
        return s

    def mul(self, poly):
        if poly.is_null or self.is_null:
            return self.null

        res = self.null
        for deg_p, c_p in poly.internal.items():
            a = self._coefficients
            for deg_a, c_a in a.items():
                deg = deg_a + deg_p
                if deg in res._coefficients.keys():
                    res.set_term(deg, res[deg] + c_p * c_a)
                else:
                    res.set_term(deg, c_p * c_a)

        return res

    def mul_constant(self, k):
        s = self.null
        if k != self.base_field.zero:
            for deg, c in self._coefficients.items():
                s.set_term(deg, k * c)
        return s

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            return self.mul(other)
        elif isinstance(other, list):
            return self.mul(Polynomial(other, base_field=self.base_field))
        else:
            return self.mul_constant(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def sub(self, poly):
        assert isinstance(poly, Polynomial)
        return self.add(-poly)

    def __sub__(self, other):
        if isinstance(other, (int, float,)):
            return self.add_constant(-other)
        elif isinstance(other, list):
            return self.sub(Polynomial(other, base_field=self.base_field))
        else:
            return self.sub(other)

    def __pow__(self, n, modulo=None):
        return self.pow(n)

    def pow(self, n):

        assert n >= 0
        if self.is_unit:
            return self.unit

        if n == 0:
            return self.unit

        if n == 1:
            return self.copy

        if self.is_null:
            return self.null

        return power(self.copy, n)

    def __hash__(self):
        return hash(tuple(self._coefficients.items()))

    def long_division(self, divisor):
        """Defines the long division according to decreasing degrees
        Be careful if the coefficients are from a ring"""

        quotient = self.null
        remainder = self.copy
        while remainder != self.null and remainder.degree >= divisor.degree:
            deg = remainder.degree - divisor.degree
            if remainder.leading % divisor.leading == self.base_field.zero:
                c = remainder.leading // divisor.leading
            else:
                raise ValueError(f'{remainder.leading} is not divisible by {divisor.leading}')

            if deg == 0:
                poly = self.unit.mul_constant(c)
            else:
                poly = self.monic(deg).mul_constant(c)
            quotient += poly
            remainder -= (divisor * poly)

        return quotient, remainder

    def long_division_reversed(self, divisor):
        """Defines the long division according to increasing degrees"""

        quotient = self.null
        remainder = self
        while remainder != self.null and remainder.valuation <= divisor.degree:

            deg = remainder.valuation - divisor.valuation
            if remainder.trailing % divisor.trailing == 0:
                c = self.base_field(remainder.trailing // divisor.trailing)
            else:
                raise ValueError(f'{remainder.trailing} is not divisible by {divisor.trailing}')
            if deg == 0:
                poly = self.unit.mul_constant(c)
            else:
                poly = self.monic(deg).mul_constant(c)
            quotient += poly
            remainder -= (divisor * poly)

        return quotient, remainder

    def __truediv__(self, other):
        return self.__floordiv__(other)

    def __floordiv__(self, other):
        if isinstance(other, self.__class__):
            return self.long_division(other)[0]
        else:
            other = self.base_field.element(other)
            return self.mul_constant(self.base_field.one / other)

    def __divmod__(self, other):
        return self.long_division(other)

    def __mod__(self, other):
        return self.long_division(other)[1]

    def __repr__(self) -> str:
        s = f'{repr(self.base_field)}.polynomial('
        s += f'{", ".join([repr(c) for c in self.coefficients])}, '
        s += f'indeterminate="{self.indeterminate}")'
        return s

    def _format_coefficient(self, c, display_plus_sign=False, raw=False):
        sf = ''
        if isinstance(c, int):
            if raw:
                sf = str(c)
            else:
                if display_plus_sign:
                    if c < 0:
                        sf += ' - '
                    else:
                        sf += ' + '
                    if abs(c) != 1:
                        sf += f'{abs(c)}'
                else:
                    if c < 0:
                        sf += '-'
                    if abs(c) != 1:
                        sf += f'{abs(c)}'
        else:
            # rely on the override of __format__
            sc = format(c, 'short')
            if raw:
                sf = sc
            else:
                if display_plus_sign:
                    if sc[0] == '-':
                        sf += ' - '
                    else:
                        sf += ' + '
                    # abs may not be defined but neg must be
                    if c != self.base_field.one and c != self.base_field.one.__neg__():
                        if sc[0] == '-':
                            sf += sc[1:]
                        else:
                            sf += sc
                else:
                    if sc[0] == '-':
                        sf += '-'
                    if c != self.base_field.one and c != self.base_field.one.__neg__():
                        if sc[0] == '-':
                            sf += sc[1:]
                        else:
                            sf += sc
        return sf

    def __str__(self) -> str:
        s = ''
        if self == self.null:
            return '0'

        for deg in sorted(self._coefficients.keys()):
            c = self[deg]
            if c != self.base_field.zero:
                if deg == 0:
                    s += self._format_coefficient(c, raw=True)
                else:
                    s += self._format_coefficient(c, len(s) > 0)
                    s += f'{self.indeterminate}'
                    if deg > 1:
                        s += f'^{deg}'
        return s

    def __call__(self, *args):
        """Syntactic sugar to create a polynomial from another one
        example: p = poly(1, 2, 1) -> p == 1 + 2X + X^2"""
        return Polynomial(list(args), base_field=self.base_field, indeterminate=self.indeterminate)

    def evaluate(self, value):
        value = self.base_field.element(value)
        # f(k) = f(x) % (x-k)
        # f(x) = q(x)(x - k) + r (deg(r) < d(x-k) = 1 => deg(r)=0
        # f(k) = r

        r = self % self.base_field.linear_polynomial(value)
        return r.trailing

    def formal_derivative(self):
        res = self.null
        for deg, c in self.internal.items():
            if deg > 0:
                res.set_term(deg - 1, self.base_field.ext_mul(deg, c))
        return res

    def gcd(self, p):
        return gcd(self, p)

    def check_irreducibility(self) -> bool:
        """Inspired by https://jeremykun.com/2014/03/13/programming-with-finite-fields/"""
        p = self.copy
        q = p.base_field.characteristic
        if q == 0:
            raise NotImplementedError(f'Cannot check polynomial irreducibility in {self.base_field}')

        x = p.monic(1)
        term = x.copy

        for _ in range(p.degree // 2):
            term = term ** q % p
            if not (term - x).is_null:
                if p.gcd(term - x).degree > 0:
                    return False
            else:
                return False

        return True

    def __invert__(self):
        return self.frobenius_reciprocal()

    def frobenius_reciprocal(self):
        """If this polynomial can be written as R^(p*m)
        where R is  polynomial, p the field characteristic and m an integer
        then taking the p-th root <=> returning R^m"""
        assert self.base_field.characteristic > 0
        if self.base_field.characteristic > 0:
            """if p is a power of a multiple of the field characteristic
            the function returns the base polynomial that is:
            if p == x^q where q % self.c == 0 it returns x^(q/self.c)"""
            assert hasattr(self.base_field, 'frobenius_reciprocal')
            p_th_root_func = self.base_field.frobenius_reciprocal

            if not self.formal_derivative().is_null:
                raise ValueError(f'The polynomial p is not a {self.base_field.characteristic} power')
            else:
                res = self.null

                for deg in range(0, self.degree + 1):
                    if deg == 0:
                        res += p_th_root_func(self.constant)
                    else:
                        if deg % self.base_field.characteristic == 0:
                            term = self.monic(deg // self.base_field.characteristic)
                            term *= p_th_root_func(self[deg])
                            res += term
                        else:
                            assert self[deg] == 0

                return res
        else:
            raise RuntimeError(f'{self.base_field} does not support taking the p-th root of a polynomial')

    @staticmethod
    def parse(expr, base_field, indeterminate='X'):
        return symbolic_polynomial(expr, base_field, indeterminate=indeterminate)
