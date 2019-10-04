import operator
from collections import namedtuple
from enum import Enum
from re import finditer
from typing import Optional, Iterator, TypeVar


Polynomial = TypeVar('Polynomial')


def symbolic_polynomial(expression: str, base_field, indeterminate: Optional[str] = 'X'):
    """Returns a polynomial from its algebraic expression
    :param expression: an algebraic expression in the indeterminate
    :param base_field: the field (or the ring) that coefficients are to be drawn from
    :param indeterminate: a single digit string in the range [a-zA-Z]
    :return: a polynomial"""
    return PolynomialParser.parse(expression, base_field, indeterminate=indeterminate)


__all__ = ['symbolic_polynomial']


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
    def __init__(self, base_field, indeterminate: str):
        self.base_field = base_field
        self.indeterminate = indeterminate

        self._stack = []
        self._stack.append(self.base_field.polynomial(self.base_field.zero, indeterminate=self.indeterminate))

    def accumulate_neutral(self, *_):
        self._stack.append(self.base_field.neutral)

    def accumulate_token(self, tok, *_):
        self._stack.append(tok.value)

    def accumulate_element(self, tok, *_):
        self._stack.append(self.base_field(tok.value))

    def accumulate_literal(self, _, v):
        self._stack.append(v)

    def accumulate_subexpression(self, tok, *_):
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
    def parse(expression: str, base_field, indeterminate: Optional[str] = 'X'):
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
