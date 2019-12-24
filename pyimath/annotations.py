__all__ = [
    'Polynomial',
    'BaseField',
    'BaseNumber',
    'Operand',
]

from typing import Union, TypeVar

Polynomial = TypeVar('Polynomial')
"""User-defined annotation to hint signatures without importing `pyimath.polynomial`"""
BaseField = TypeVar('BaseField')
"""User-defined annotation to hint signatures without actually defining  a `BaseField` class"""
BaseNumber = TypeVar('BaseNumber')
"""User-defined annotation to hint signatures without actually defining  a `BaseNumber` class"""
Operand = Union[Polynomial, BaseNumber]
"""User-defined annotation to hint signatures with a generic `Operand` pseudo-type"""

