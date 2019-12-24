Module pyimath.factorize
========================

Functions
---------

    
`factorize(p)`
:   Main entry point of the module, provide an instance of the Factorization class as a place holder
    to call any factorization algorithms : `square_free`, `distinct_degree`, `equal_degree` or `cantor_zassenhaus`.
    
    Accepts an instance of `Polynomial` (parameter `p`) and returns a instance of `Factorization`

Classes
-------

`Factor(*args, **kwargs)`
:   Defines a `Factor` class that wraps a `Polynomial` along with a `multiplicity`.
    
    Works a `namedtuple`with the following fields:
    
    * `value`: an instance of `Polynomial`
    * `multiplicity`: an integer representing the multiplicity of this factor in a factorization
    * `max_degree`: for an reducible factor, represents the maximum degree of its irreducible factors

    ### Ancestors (in MRO)

    * builtins.tuple

    ### Instance variables

    `is_irreducible`
    :   returns `True` if the `Polynomial` referenced by `value`is irreducible otherwise returns `False`

`Factorization(base_field, poly)`
:   Service provider for polynomial factorization.
    
    The main service is the full factorization by the Cantor-Zassenhaus algorithm
    
    Never use it directly, subclass it if you want but in any case, please use `factorize(poly).<method>()`

    ### Methods

    `cantor_zassenhaus(self)`
    :   Full factorisation algorithm for any polynomial over a finite field.
        Returns a tuple containing a constant term and a list of factors with their multiplicity

    `distinct_degree(self)`
    :   Distinct degree factorisation of a square free monic polynomial.
        Returns a list of Factors with a multiplicity of 1 and the maximum degree of any irreducible factor

    `equal_degree(self, nb_factors, max_degree)`
    :   Cantor–Zassenhaus factorization algorithm.
        Input: Over a finite field Fq of odd order q,
               p is a monic square free polynomial in Fq[x] of degree n = rd,
               which has r ≥ 2 irreducible factors each of degree d
        Output: The set of monic irreducible factors of f.

    `factors_product(self, factors)`
    :   Converts a list of factors into a polynomial. The `factors` parameter must be an iterable of `Factor` instances.
        Returns a `Polynomial` that is the product of these factors along with their multiplicity

    `square_free(self)`
    :   Square free factorization of a monic polynomial over a finite field of prime characteristic.
        Returns the square free part and a partial or complete factorization