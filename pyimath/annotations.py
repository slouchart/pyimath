from typing import Union, TypeVar

Polynomial = TypeVar('Polynomial')
BaseField = TypeVar('BaseField')
BaseNumber = TypeVar('BaseNumber')
Operand = Union[Polynomial, BaseNumber]
