Module pyimath.gaussint
=======================

Classes
-------

`GaussInt(x=0, y=0)`
:   Represents a Gaussian Integer

    ### Instance variables

    `conjugate`
    :   Returns the complex conjugate of a gaussian integer

    `x`
    :   Returns the integer component of a gaussian integer

    `y`
    :   Returns the complex component of a gaussian integer

`GaussianIntegerRing()`
:   Set of numbers of the form `x + i*y` where `i^2 = -1`. Also defined as `Z[X]/<X^2+1>`

    ### Static methods

    `element(*args)`
    :   Returns a gaussian integer from a tuple of integers

    `pow(a, n)`
    :   Returns the n-th power of a gaussian integer

    ### Instance variables

    `characteristic`
    :   Returns the ring characteristic which is always zero.

    `neutral`
    :   Same as `self.one`

    `null`
    :   Same as `self.zero`

    `one`
    :   Returns the multiplicative neutral

    `zero`
    :   Returns the additive neutral

    ### Methods

    `add(self, a, b)`
    :   Returns the sum of two gaussian integers

    `additive_inverse(self, a)`
    :   Returns the additive inverse of a gaussian integer

    `divmod(self, a, b)`
    :   Returns the quotient and the remainder of the euclidean division of two gaussian integers

    `element_from_polynomial(self, p)`
    :   Returns the element `(a, b)` from the linear integer polynomial `a + bZ`

    `ext_mul(self, n, a)`
    :   Returns the product of a gaussian integer and a regular integer

    `floor_div(self, a, b)`
    :   Returns the quotient of the euclidean division of two gaussian integers

    `format_element(self, e)`
    :   Returns a printable representation of a gaussian integer

    `linear_polynomial(self, e)`
    :   Returns the polynomial `X - e`

    `mod(self, a, b)`
    :   Returns the remainder of the division of two gaussian integers

    `mul(self, a, b)`
    :   Returns the product of two gaussian integers

    `parse(self, expr)`
    :   Returns a gaussian integer from its symbolic expression

    `polynomial(self, coefficients, indeterminate='z')`
    :   Returns a polynomial from its coefficients

    `polynomial_from_element(self, e)`
    :   Returns the polynomial `yZ + x` from the element `'(x, y)`