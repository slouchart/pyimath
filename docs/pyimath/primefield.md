Module pyimath.primefield
=========================

Classes
-------

`PFElement(field, n)`
:   Represents a single element from a prime field by, basically, duck-typing an integer

    ### Instance variables

    `is_one`
    :   Return `True` if the element is the multiplicative neutral of its field, `False` otherwise

    `is_zero`
    :   Returns `True` if the element is the additive neutral of its field, `False` otherwise

`PrimeField(prime)`
:   Defines a prime field definition of characteristic P.
    
    The field elements are not represented by integer modulo P
    but rather by signed integers in the range `-(P-1)/2..(P-1)/2`

    ### Instance variables

    `neutral`
    :   Same as `self.one`

    `null`
    :   Same as `self.zero`

    `one`
    :   Returns the neutral of the multiplicative group

    `zero`
    :   Returns the neutral of the additive group

    ### Methods

    `add(self, a, b)`
    :   Returns the sum of two elements

    `additive_inverse(self, a)`
    :   Returns the additive inverse of an element

    `div(self, a, b)`
    :   Returns the quotient of two elements

    `divmod(self, a, b)`
    :   Kept for interface purposes. Always returns `(a / b, self.zero)`

    `element(self, n)`
    :   Casts an integer into an element of the field

    `ext_mul(self, n, a)`
    :   Returns the n-th iterated sum of an element
        i.e  `a + ... + a, n times`

    `floor_div(self, a, b)`
    :   Same as `self.div`

    `frobenius_reciprocal(self, a)`
    :   Returns the Frobenius reciprocal of an element. As a remainder, the Frobenius automorphism is `a -> a^q`
        where `q` is the characteristic of the field. In prime fields, it is the identity.

    `generate_irreducible_polynomial(self, degree, max_retries=15)`
    :   Returns an irreducible polynomial of a given degree over the base field. This algorithm is not deterministic
        and may raise exceptions. `degree` is the degree of the irreducible polynomial and `max_retries` sets
        the maximum number of attempts.

    `linear_polynomial(self, e)`
    :   Returns the polynomial `X + (-e)`

    `mod(self, *_)`
    :   Kept for interface purposes. Always return `self.zero`

    `mul(self, a, b)`
    :   Returns the product of two elements

    `multiplicative_inverse(self, a)`
    :   Returns the multiplicative inverse of an element

    `parse_poly(self, expr)`
    :   Returns a polynomial from its symbolic expression

    `polynomial(self, *args, indeterminate='X')`
    :   Returns a polynomial from the arguments

    `pow(self, a, n)`
    :   Return the n-th power of an element

    `random_element(self)`
    :   Returns an element of the field at random

    `random_polynomial(self, degree)`
    :   Returns a random polynomial of a given degree