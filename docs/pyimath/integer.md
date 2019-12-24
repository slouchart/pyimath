Module pyimath.integer
======================

Classes
-------

`IntegerRing()`
:   Mimics a Field as defined in `primefield` or `finitefield` modules (and only mimics a field, Z is a ring!)
    
    Mainly used to define polynomials from `Z[X]`

    ### Static methods

    `add(a, b)`
    :   Returns the sum of two integers

    `additive_inverse(a)`
    :   Returns the additive inverse of an integer

    `divmod(a, b)`
    :   Returns the quotient and the remainder of the integer division

    `element(n)`
    :   Element constructor. Kept for interface compatibility

    `ext_mul(n, a)`
    :   Kept for interface purposes. Simply returns `n*a`

    `floor_div(a, b)`
    :   Return the quotient of the integer division

    `mod(a, b)`
    :   Returns the remainder of the integer division

    `mul(a, b)`
    :   Return the product of two integers

    `pow(a, e)`
    :   Returns `a` to the `e`-th power

    ### Instance variables

    `characteristic`
    :   Returns the characteristic of Z which is always 0 by definition

    `neutral`
    :   Same as `self.one`, returns the multiplicative neutral element of Z, 1

    `null`
    :   Same as `self.zero`, returns the additive neutral element of Z, 0

    `one`
    :   Returns the multiplicative neutral element of Z

    `zero`
    :   Returns the additive neutral element of Z

    ### Methods

    `linear_polynomial(self, e)`
    :   Returns a polynomial defined as `X - e`

    `polynomial(self, *args, indeterminate='X')`
    :   Returns a polynomial built from a list of coefficients.