Module pyimath.finitefield
==========================

Functions
---------

    
`finite_field(order)`
:   Returns the finite field of the given order.
    Fields that may be instantiated this way are:
    
    * all prime fields of order less then 23
    * finite fields of non-prime order up to order 27. These finite fields have all a valid generator.

Classes
-------

`FFElement(field, v)`
:   Represents an element of a finite field of non-prime order as a vector
    of elements of the base prime field

    ### Instance variables

    `is_scalar`
    :   Returns `True` if the element belongs to the base prime field, `False` otherwise

    `null`
    :   Returns `True` if the element is the additive neutral of the field

`FiniteField(prime, dimension, ideal, generator=None, root_symbol='j')`
:   Represents a non-prime finite field where:
        
    
    `prime` is the order of the base prime field,
    
    `dimension` is the dimension of the finite field seen as a vector space over its prime field,
    
    `ideal` is the minimal polynomial of the adjunct root j, must be irreducible over the prime field,
    
    `generator` must be an element of the field that generates the multiplicative group
    
    and `root_symbol` is a one-character string to represent the adjunct root

    ### Instance variables

    `basis`
    :   Returns a vector-space basis of field elements

    `characteristic`
    :   Returns the characteristic of the field

    `frobenius_map`
    :   Returns what is actually the linear map of the inverse function of Frobenius automorphism
        (see MISCELLANEOUS.md)

    `has_valid_generator`
    :   Returns `True` if the field was set with a valid generator i.e. an element whose order is equal
        to the order of the field

    `neutral`
    :   Same as `self.one`, returns the multiplicative neutral element of the field

    `null`
    :   Same as `self.zero`, returns the additive neutral element of the field

    `one`
    :   Returns the multiplicative neutral element of the field

    `order`
    :   Returns the order of the field

    `prime_field`
    :   Returns the instance of the base prime field

    `zero`
    :   Returns the additive neutral element of the field

    ### Methods

    `add(self, a, b)`
    :   Adds two elements

    `additive_inverse(self, a)`
    :   Returns the additive inverse of an element

    `check_irreducible_polynomial(self, p)`
    :   Returns `True` if the polynomial is irreducible, `False` otherwise

    `div(self, a, b)`
    :   Returns the true division of two elements

    `divmod(self, a, b)`
    :   Kept for interface purpose, returns always `(a / b, self.zero)`

    `element(self, v)`
    :   Returns an instance of a `FFElement` from a vector

    `element_from_polynomial(self, p)`
    :   Return an element whose basis components are the coefficients of a given polynomial over the field
        e.g. The polynomial `1 + X + X^2` would yield `(1+j+j^2)`

    `element_order(self, e)`
    :   Returns the order of an element, that is the minimal integer k such as e^k = 1

    `ext_mul(self, n, a)`
    :   Returns the n-th iterated addition of an element with itself

    `find_generator(self, set_generator=False)`
    :   Returns a valid generator of the multiplicative group of the field

    `format_element(self, e, format_spec='')`
    :   Returns a printable representation of an element

    `frobenius_reciprocal(self, a)`
    :   Returns the bijective inverse of an element through the Frobenius automorphism. This is
        equivalent of taking the p-th root of an element where p is the field characteristic

    `linear_polynomial(self, e)`
    :   Returns the polynomial `X - e`

    `mod(self, _, b)`
    :   Kept for interface purpose, returns always `self.zero`

    `mul(self, a, b)`
    :   Returns the product of two elements

    `multiplicative_inverse(self, a)`
    :   Returns the multiplicative inverse of an element

    `parse_poly(self, expr)`
    :   Returns a polynomial from its symbolic expression

    `polynomial(self, *args, indeterminate='X')`
    :   Returns a polynomial from a list of coefficients. The coefficients are converted to `FFElement`

    `polynomial_from_element(self, e)`
    :   Returns a polynomial whose coefficients are those of the element
        e.g. `(1+j+j^2)` would yield `1 + X + X^2`

    `pow(self, a, n)`
    :   Returns the n-th power of an element