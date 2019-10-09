from typing import Iterable, Optional, Dict, Collection, Tuple, Sequence, Iterator, Any
import operator
from collections import namedtuple
from enum import Enum
from re import finditer


from pyimath.annotations import BaseField, BaseNumber, Operand
from pyimath.functions import gcd, reduce_to_gcd, power


__all__ = ['Polynomial', 'symbolic_polynomial']


class Polynomial:
    """Represents a polynomial over a finite field or over the integers

    More generally, definition of polynomials over a ring is possible
    but not recommended

    The internal representation of the coefficients uses a `dict` that
    indexes the coefficients by their degree. The usual definition of a polynomial as a sequence of numbers
    that are all zeroes from a certain index is used to initialize the polynomial
    on instantiation and can be retrieved through the `coefficients` property

    Moreover this class implements all the operations over a ring of polynomials

    `base_field` must represent an integral domain that is:

    - a set which is an abelian group for `+` and a semi-group for `*` and where `GCD` are computable.

    - a ring of polynomials is at least defined over an integral domain

    examples: Z, Z(i)
    """

    def __init__(self, coeffs: Sequence[BaseNumber], base_field: BaseField, indeterminate: Optional[str] = 'X'):
        """`coeffs` is an iterable of elements from the base field, `base_field` an instance of what should generally be
        a finite field and `indeterminate` is a single digit string used to format the polynomial
        """
        self.base_field = base_field
        assert hasattr(base_field, 'element') and hasattr(base_field, '__call__')

        self._coefficients = self._safe_convert_coefficients(self._remove_trailing_zeros(coeffs))
        self.indeterminate = indeterminate

    def add(self, poly: 'Polynomial') -> 'Polynomial':
        """Addition in a ring of polynomials"""
        a = self._coefficients
        s = poly.copy
        for deg_a, c_a in a.items():
            if deg_a in s.internal.keys():
                s._set_term(deg_a, s[deg_a] + c_a)
            else:
                s._set_term(deg_a, c_a)

        return s

    def add_constant(self, k: BaseNumber) -> 'Polynomial':
        """External addition of a polynomial"""
        a = self._coefficients
        s = self.null

        if not self.is_null:
            for deg, c in a.items():
                if deg == 0:
                    s._set_term(deg, c + k)
                else:
                    s._set_term(deg, c)
        else:
            s._set_term(0, k)

        return s

    def check_irreducibility(self) -> bool:
        """Returns True if the polynomial is irreducible
        Inspired by https://jeremykun.com/2014/03/13/programming-with-finite-fields/"""
        p = self.copy
        q = p.base_field.characteristic
        if q == 0:
            raise NotImplementedError(f'Cannot check polynomial irreducibility in {self.base_field}')

        x = p.monic(1)
        term = x.copy

        for _ in range(p.degree // 2):
            term = term ** q % p
            if not (term - x).is_null:
                if gcd(p, term - x).degree > 0:
                    return False
            else:
                return False

        return True

    @property
    def coefficients(self) -> Collection:
        """Returns the coefficients as a list"""
        deg = 0
        res = []
        while deg <= self.degree:
            res.append(self[deg])
            deg += 1
        return res

    @property
    def constant(self) -> BaseNumber:
        """Returns the value of the coefficient of the term of degree 0"""
        if not self.is_null:
            return self[0]
        else:
            return self.base_field.zero

    @property
    def copy(self) -> 'Polynomial':
        """Returns a copy of itself"""
        res = self.null
        for deg, c in self.internal.items():
            res._set_term(deg, c)
        return res

    @property
    def degree(self) -> int:
        """Returns the degree of the polynomial"""
        if self.is_null:
            return 0  # rigorously, it should be -infinity
        else:
            return max(self._coefficients.keys())

    def evaluate(self, value: BaseNumber) -> BaseNumber:
        """Evaluate the polynomial for some value"""
        value = self.base_field.element(value)
        # f(k) = f(x) % (x-k)
        # f(x) = q(x)(x - k) + r (deg(r) < d(x-k) = 1 => deg(r)=0
        # f(k) = r

        r = self % self.base_field.linear_polynomial(value)
        return r.trailing

    def formal_derivative(self) -> 'Polynomial':
        """Computes and returns the formal derivative of a polynomial"""
        res = self.null
        for deg, c in self.internal.items():
            if deg > 0:
                res._set_term(deg - 1, self.base_field.ext_mul(deg, c))
        return res

    def gcd(self, p: 'Polynomial') -> 'Polynomial':
        """Returns the GCD of two polynomials"""
        return gcd(self, p)

    def frobenius_reciprocal(self) -> 'Polynomial':
        """If this polynomial can be written as R^(p*m)
        where R is a polynomial, p the field characteristic and m an integer
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

    @property
    def internal(self) -> Dict[int, BaseNumber]:
        """Returns the coefficients as a dict indexed by their degree
        Eschews null terms"""
        return dict({deg: c for deg, c in self._coefficients.items()})

    @property
    def is_abs_unit(self) -> bool:
        """Returns True if the polynomial is contant of constant term 1 or -1"""
        return self in (self.monic(0), -self.monic(0))

    @property
    def is_constant(self) -> bool:
        """Returns True if the polynomial is of degree zero and is not null"""
        return self.degree == 0 and not self.is_null

    @property
    def is_irreducible(self) -> bool:
        """Returns True if the polynomial is irreducible"""
        return self.check_irreducibility()

    @property
    def is_monic(self) -> bool:
        """Returns True if the leading coefficient is one"""
        return self._coefficients[self.degree] == self.base_field.one

    @property
    def is_null(self) -> bool:
        """Returns True if all coefficients are zero"""
        return len(self._coefficients.keys()) == 0

    @property
    def is_unit(self) -> bool:
        """Returns True if the polynomial is constant of constant term 1"""
        return self == self.unit

    @property
    def leading(self) -> BaseNumber:
        """Returns the value of the coefficient of the term of highest degree"""
        if not self.is_null:
            return self._coefficients[self.degree]
        else:
            return self.base_field.zero

    def long_division(self, divisor: 'Polynomial') -> Tuple['Polynomial', 'Polynomial']:
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

    def long_division_reversed(self, divisor: 'Polynomial') -> Tuple['Polynomial', 'Polynomial']:
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

    def make_monic(self) -> 'Polynomial':
        """
        Attempts to divide the polynomial by its leading coefficient (for a field) or by the gcd of its
        coefficients (for a ring)
        :return: a monic polynomial or raises an error if the polynomial cannot be made monic
        """
        if not self.is_monic:
            if self.base_field.characteristic != 0:
                return self / self.leading
            else:
                g = reduce_to_gcd(iter(self.coefficients))
                if self.leading // g == self.base_field.one:
                    return self // self.leading
                else:
                    raise ValueError(f'Polynomial {self} over {self.base_field} cannot be made monic')

        else:
            return self.copy

    def monic(self, degree: int = 1) -> 'Polynomial':
        """Returns a monic polynomial with a single term of a given degree"""
        res = self.null
        res._set_term(degree, self.base_field.one)
        return res

    def mul(self, poly: 'Polynomial') -> 'Polynomial':
        """Multiplication in a ring of polynomials"""
        if poly.is_null or self.is_null:
            return self.null

        res = self.null
        for deg_p, c_p in poly.internal.items():
            a = self._coefficients
            for deg_a, c_a in a.items():
                deg = deg_a + deg_p
                if deg in res._coefficients.keys():
                    res._set_term(deg, res[deg] + c_p * c_a)
                else:
                    res._set_term(deg, c_p * c_a)

        return res

    def mul_constant(self, k: BaseNumber) -> 'Polynomial':
        """External multiplication (vector space external product)"""
        s = self.null
        if k != self.base_field.zero:
            for deg, c in self._coefficients.items():
                s._set_term(deg, k * c)
        return s

    @property
    def null(self) -> 'Polynomial':
        """Returns the null polynomial"""
        return Polynomial([], base_field=self.base_field, indeterminate=self.indeterminate)

    @staticmethod
    def parse(expr: str, base_field: BaseField, indeterminate: Optional[str] = 'X') -> 'Polynomial':
        """Returns a polynomial from its algebraic expression
        """
        return symbolic_polynomial(expr, base_field, indeterminate=indeterminate)

    def pow(self, n: int) -> 'Polynomial':
        """Exponentiation of a polynomial"""
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

    def sub(self, poly: 'Polynomial') -> 'Polynomial':
        """Subtraction is indeed just an addition of an inverse"""
        assert isinstance(poly, Polynomial)
        return self.add(-poly)

    @property
    def trailing(self) -> BaseNumber:
        """Returns the value of the coefficient of the term of lowest degree"""
        if self.is_null:
            return self.base_field.zero
        else:
            return self._coefficients[self.valuation]

    @property
    def unit(self) -> 'Polynomial':
        """Return 1 as a polynomial of degree 0"""
        return self.monic(0)

    @property
    def valuation(self) -> int:
        """Returns the degree of the term of lowest degree"""
        if self.is_null:
            raise ValueError('The valuation of the null polynomial is undefined')

        return min(self._coefficients.keys())

    def __eq__(self, other: Operand) -> bool:
        """Term-wise comparison of two polynomials"""
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

    def __getitem__(self, degree: int) -> BaseNumber:
        """Returns the coefficient of the term of a given degree"""
        if degree in self._coefficients.keys():
            return self._coefficients[degree]
        else:
            return self.base_field.zero

    def __len__(self) -> int:
        """Returns the number of non zero terms"""
        return len([c for c in self._coefficients.values() if c != 0])

    def __add__(self, other: Operand) -> 'Polynomial':
        if isinstance(other, Polynomial):
            return self.add(other)
        elif isinstance(other, list):
            return self.add(Polynomial(other, base_field=self.base_field))
        else:
            return self.add_constant(other)

    def __radd__(self, other: Operand) -> 'Polynomial':
        return self.__add__(other)

    def __neg__(self) -> 'Polynomial':
        """Returns the inverse of a polynomial with respect to addition"""
        a = self._coefficients
        s = self.null
        for deg, c in a.items():
            s._set_term(deg, -c)
        return s

    def __mul__(self, other: Operand) -> 'Polynomial':
        if isinstance(other, Polynomial):
            return self.mul(other)
        elif isinstance(other, list):
            return self.mul(Polynomial(other, base_field=self.base_field))
        else:
            return self.mul_constant(other)

    def __rmul__(self, other: Operand) -> 'Polynomial':
        return self.__mul__(other)

    def __sub__(self, other: Operand) -> 'Polynomial':
        if isinstance(other, (int, float,)):
            return self.add_constant(-other)
        elif isinstance(other, list):
            return self.sub(Polynomial(other, base_field=self.base_field))
        else:
            return self.sub(other)

    def __pow__(self, n: int, modulo: 'Polynomial' = None) -> 'Polynomial':
        return self.pow(n)

    def __hash__(self) -> int:
        """Allows a polynomial to become a dictionary key"""
        return hash(tuple(self._coefficients.items()))

    def __truediv__(self, other: Operand) -> 'Polynomial':
        return self.__floordiv__(other)

    def __floordiv__(self, other: Operand) -> 'Polynomial':
        if isinstance(other, self.__class__):
            return self.long_division(other)[0]
        else:
            other = self.base_field.element(other)
            return self.mul_constant(self.base_field.one / other)

    def __divmod__(self, other: Operand) -> Tuple['Polynomial', 'Polynomial']:
        return self.long_division(other)

    def __mod__(self, other: Operand) -> 'Polynomial':
        return self.long_division(other)[1]

    def __repr__(self) -> str:
        s = f'{repr(self.base_field)}.polynomial('
        s += f'{", ".join([repr(c) for c in self.coefficients])}, '
        s += f'indeterminate="{self.indeterminate}")'
        return s

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

    def __call__(self, *args) -> 'Polynomial':
        """Syntactic sugar to create a polynomial from another one
        example: p = poly(1, 2, 1) -> p == 1 + 2X + X^2"""
        return Polynomial(list(args), base_field=self.base_field, indeterminate=self.indeterminate)

    def __invert__(self) -> 'Polynomial':
        """Return the Frobenius reciprocal with the operator ~"""
        return self.frobenius_reciprocal()

    # Gory Details (as usual)

    def _format_coefficient(self, c: BaseNumber, display_plus_sign: bool = False, raw: bool = False) -> str:
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

    def _remove_trailing_zeros(self, seq: Sequence) -> Collection:
        if len(seq) == 0:
            return []
        revseq = list(reversed(seq))
        while revseq[0] == self.base_field.zero:
            revseq = revseq[1:]
            if len(revseq) == 0:
                return []

        return list(reversed(revseq))

    def _safe_convert_coefficients(self, seq: Iterable) -> Dict[int, BaseNumber]:
        bf = self.base_field
        return dict({deg: bf.element(c) for deg, c in enumerate(seq) if c != bf.zero})

    def _set_term(self, deg: int, c: BaseNumber):
        if c == self.base_field.zero:
            if deg in self._coefficients.keys():
                del self._coefficients[deg]
        else:
            self._coefficients[deg] = c


def symbolic_polynomial(expression: str, base_field: BaseField, indeterminate: Optional[str] = 'X'):
    """Returns a polynomial from its algebraic expression
    :param expression: an algebraic expression in the indeterminate
    :param base_field: the field (or the ring) that coefficients are to be drawn from
    :param indeterminate: a single digit string in the range [a-zA-Z]
    :return: a polynomial"""
    return PolynomialParser.parse(expression, base_field, indeterminate=indeterminate)


"""Lexer et Parser code follows, should not be exported"""


class Lexer:
    INDETERMINATE = 'INDETERMINATE'
    INTEGER = 'INTEGER'
    OPERATOR = 'OPERATOR'
    EXPONENT = 'EXPONENT'
    SUBEXPR = 'SUBEXPR'
    EOF = 'EOF'
    IGNORE = 'IGNORE'
    MISMATCH = 'MISMATCH'

    class Token(namedtuple('Token', 'type, value, position')):
        __slots__ = ()

    def __init__(self, indeterminate: Optional[str] = 'X', root_symbol: Optional[str] = 'j'):
        self.indeterminate = indeterminate
        self.root_symbol = root_symbol

        self.symbols = [
            (Lexer.INDETERMINATE, r'[%s]' % self.indeterminate),
            (Lexer.INTEGER, r'[1-9][0-9]*'),
            (Lexer.OPERATOR, r'[+-]'),
            (Lexer.EXPONENT, r'[\^]'),
            (Lexer.SUBEXPR, r'\([^)]+\)'),
            (Lexer.EOF, r'$'),
            (Lexer.IGNORE, r'\s'),  # must stay before the last item
            (Lexer.MISMATCH, r'.')  # must stay the last item
        ]
        self.tokens_re = '|'.join([f'(?P<{tok}>{re})' for tok, re in self.symbols])

    def lex(self, expression: str) -> Iterator['Lexer.Token']:

        for tok in finditer(self.tokens_re, expression):
            token = Lexer.Token(tok.lastgroup, tok.group(), tok.start())

            if token.type == Lexer.IGNORE:
                continue
            elif token.type == Lexer.SUBEXPR:
                # remove left and right parentheses
                lexer = Lexer(indeterminate=self.root_symbol, root_symbol='')
                yield Lexer.Token(token.type, lexer.lex(token.value[1:-1], ), token.position)

            elif token.type == Lexer.INTEGER:
                yield Lexer.Token(token.type, int(token.value), token.position)

            elif token.type == Lexer.OPERATOR:
                if token.value == '+':
                    yield Lexer.Token(token.type, operator.add, token.position)
                elif token.value == '-':
                    yield Lexer.Token(token.type, operator.sub, token.position)

            else:
                yield token


class ParsingContext:
    def __init__(self, base_field: BaseField, indeterminate: str):
        self.base_field = base_field
        self.indeterminate = indeterminate

        self._stack = []
        self._stack.append(self.base_field.polynomial(self.base_field.zero, indeterminate=self.indeterminate))

    def accumulate_neutral(self, *_):
        self._stack.append(self.base_field.neutral)

    def accumulate_token(self, tok: Lexer.Token, *_):
        self._stack.append(tok.value)

    def accumulate_element(self, tok: Lexer.Token, *_):
        self._stack.append(self.base_field(tok.value))

    def accumulate_literal(self, _, v: Any):
        self._stack.append(v)

    def accumulate_subexpression(self, tok: Lexer.Token, *_):
        def convert_subexpr(subexpr):
            ctx = ParsingContext(self.base_field.prime_field, indeterminate=self.base_field.root_symbol)
            return self.base_field.element_from_polynomial(PolynomialParser.start(subexpr, ctx))

        self._stack.append(convert_subexpr(tok.value))

    def reduce(self, *_):
        try:
            (result, op, coefficient, degree), self._stack = self._stack[-4:], self._stack[:-4]
            self._stack.append(op(result, result.monic(degree).mul_constant(coefficient)))
        except BaseException as e:
            raise RuntimeError(e)

    def get_result(self) -> 'Polynomial':
        return self._stack.pop()


class PolynomialParser:

    @staticmethod
    def parse(expression: str, base_field: BaseField, indeterminate: Optional[str] = 'X'):
        """Main parsing utility"""
        if hasattr(base_field, 'root_symbol'):
            lexer = Lexer(indeterminate=indeterminate, root_symbol=base_field.root_symbol)
        else:
            lexer = Lexer(indeterminate=indeterminate, root_symbol='')

        ctx = ParsingContext(base_field, indeterminate)
        return PolynomialParser.start(lexer.lex(expression), ctx)

    @staticmethod
    def start(lexer: Iterator, context: ParsingContext) -> Polynomial:

        def format_syntax_error(s, t):
            return f'Syntax error at {t.position}: unexpected {t} in state {s}'

        state = PolynomialParser.States.starting

        for token in lexer:
            if (state, token.type) in PolynomialParser.transitions:
                transition = PolynomialParser.transitions[(state, token.type)]
                for op, args in transition.actions:
                    op(context, token, *args)
                state = transition.next_state
            else:
                raise SyntaxError(format_syntax_error(state, token))

        assert state == PolynomialParser.States.complete
        return context.get_result()

    """Internals"""
    Transition = namedtuple('Transition', 'actions next_state')

    class States(Enum):
        starting = 0
        coefficient = 1
        indeterminate = 2
        exponent = 3
        sub_expression = 4
        operator = 5
        complete = 6
        init = 7

    transitions = {
        (States.starting, Lexer.INTEGER):
            Transition([
                (ParsingContext.accumulate_literal, (operator.add,)),
                (ParsingContext.accumulate_element, (),)
            ], States.coefficient),

        (States.starting, Lexer.SUBEXPR):
            Transition([
                (ParsingContext.accumulate_literal, (operator.add,)),
                (ParsingContext.accumulate_subexpression, ())
            ], States.coefficient),

        (States.starting, Lexer.OPERATOR):
            Transition([
                (ParsingContext.accumulate_token, ())
            ], States.operator),

        (States.starting, Lexer.INDETERMINATE):
            Transition([
                (ParsingContext.accumulate_literal, (operator.add,)),
                (ParsingContext.accumulate_neutral, (),)
            ], States.indeterminate),

        (States.starting, Lexer.EOF):
            Transition([], States.complete),

        (States.coefficient, Lexer.INDETERMINATE):
            Transition([], States.indeterminate),

        (States.coefficient, Lexer.OPERATOR):
            Transition([
                (ParsingContext.accumulate_literal, (0,)),
                (ParsingContext.reduce, ()),
                (ParsingContext.accumulate_token, ()),
            ], States.operator),

        (States.coefficient, Lexer.EOF):
            Transition([
                (ParsingContext.accumulate_literal, (0,)),
                (ParsingContext.reduce, ())
            ], States.complete),

        (States.operator, Lexer.INTEGER):
            Transition([
                (ParsingContext.accumulate_element, ())
            ], States.coefficient),

        (States.operator, Lexer.SUBEXPR):
            Transition([
                (ParsingContext.accumulate_subexpression, ())
            ], States.coefficient),

        (States.operator, Lexer.INDETERMINATE):
            Transition([
                (ParsingContext.accumulate_neutral, ())
            ], States.indeterminate),

        (States.indeterminate, Lexer.EXPONENT):
            Transition([], States.exponent),

        (States.indeterminate, Lexer.OPERATOR):
            Transition([
                (ParsingContext.accumulate_literal, (1,)),
                (ParsingContext.reduce, ()),
                (ParsingContext.accumulate_token, ())
            ], States.operator),

        (States.indeterminate, Lexer.EOF):
            Transition([
                (ParsingContext.accumulate_literal, (1,)),
                (ParsingContext.reduce, ())
            ], States.complete),

        (States.exponent, Lexer.INTEGER):
            Transition([
                (ParsingContext.accumulate_token, ()),
                (ParsingContext.reduce, ())
            ], States.init),

        (States.init, Lexer.OPERATOR):
            Transition([
                (ParsingContext.accumulate_token, ())
            ], States.operator),

        (States.init, Lexer.EOF):
            Transition([], States.complete)
    }
