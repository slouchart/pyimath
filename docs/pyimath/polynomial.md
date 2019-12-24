Module pyimath.polynomial
=========================

Functions
---------

    
`symbolic_polynomial(expression, base_field, indeterminate='X')`
:   Returns a polynomial from its algebraic expression where:
    
     * `expression`is an algebraic expression in the `indeterminate`,
     * `base_field` is the field (or the ring) that coefficients are to be drawn from,
     * and `indeterminate` is a single digit string in the range [a-zA-Z], usually `'X'`.
    
    Returns an instance of `Polynomial`

Classes
-------

`Polynomial(coeffs, base_field, indeterminate='X')`
:   Represents a polynomial over a finite field or over the integers
    
    More generally, definition of polynomials over a ring is possible
    but not really recommended
    
    The internal representation of the coefficients uses a `dict` that
    indexes the coefficients by their degree. The usual definition of a polynomial as a sequence of numbers
    that are all zeroes from a certain index is used to initialize the polynomial
    on instantiation and can be retrieved through the `coefficients` property
    
    Moreover this class implements all the operations over a ring of polynomials
    
    `base_field` must represent an integral domain that is:
    
    - a set which is an abelian group for `+` and a semi-group for `*` and where `GCD` are computable.
    
    - a ring of polynomials is at least defined over an integral domain
    
    examples: Z, Z(i)
    
    `coeffs` is an iterable of elements from the base field, `base_field` an instance of what should generally be
    a finite field and `indeterminate` is a single digit string used to format the polynomial

    ### Static methods

    `parse(expr, base_field, indeterminate='X')`
    :   Returns a polynomial from its algebraic expression

    ### Instance variables

    `coefficients`
    :   Returns the coefficients as a list

    `constant`
    :   Returns the value of the coefficient of the term of degree 0

    `copy`
    :   Returns a copy of itself

    `degree`
    :   Returns the degree of the polynomial

    `internal`
    :   Returns the coefficients as a `dict` indexed by their degree.
        
        Eschews null terms

    `is_abs_unit`
    :   Returns `True` if the polynomial is a constant of constant term 1 or -1

    `is_constant`
    :   Returns `True` if the polynomial is of degree zero and is not null

    `is_irreducible`
    :   Returns `True` if the polynomial is irreducible

    `is_monic`
    :   Returns `True` if the leading coefficient is one

    `is_null`
    :   Returns `True` if all coefficients are zero

    `is_unit`
    :   Returns `True` if the polynomial is a constant of constant term 1

    `leading`
    :   Returns the value of the coefficient of the term of highest degree

    `null`
    :   Returns the null polynomial

    `trailing`
    :   Returns the value of the coefficient of the term of lowest degree

    `unit`
    :   Returns 1 as a polynomial of degree 0

    `valuation`
    :   Returns the degree of the term of lowest degree

    ### Methods

    `add(self, poly)`
    :   Returns the sum of two polynomials

    `add_constant(self, k)`
    :   External addition of a polynomial

    `check_irreducibility(self)`
    :   Returns True if the polynomial is irreducible
        
        Inspired by https://jeremykun.com/2014/03/13/programming-with-finite-fields/

    `evaluate(self, value)`
    :   Evaluate the polynomial for some value

    `formal_derivative(self)`
    :   Computes and returns the formal derivative of a polynomial

    `frobenius_reciprocal(self)`
    :   Returns a polynomial `R` if and only if  this polynomial can be written as `R^(p*m)`
        where `R` is a polynomial, `p` the field characteristic and `m` an integer.
        
        Equivalent of taking the p-th root of a polynomial over a finite field
        
        Do not use this method if the base field/ring has a characteristic of zero

    `gcd(self, p)`
    :   Returns the GCD of two polynomials

    `long_division(self, divisor)`
    :   Defines the long division according to decreasing degrees
        
        Be careful if the coefficients are from a ring

    `long_division_reversed(self, divisor)`
    :   Defines the long division according to increasing degrees

    `make_monic(self)`
    :   Attempts to divide the polynomial by its leading coefficient (for a field) or by the gcd of its
        coefficients (for a ring)
        
        Returns a monic polynomial or raises an error if the polynomial cannot be made monic

    `monic(self, degree=1)`
    :   Returns a monic polynomial with a single term of a given degree

    `mul(self, poly)`
    :   Multiplication in a ring of polynomials

    `mul_constant(self, k)`
    :   External multiplication (vector space external product) of a polynomial and a constant

    `pow(self, n)`
    :   Exponentiation of a polynomial

    `sub(self, poly)`
    :   Returns the difference between two polynomials