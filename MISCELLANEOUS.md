# Note on the principles of prime field operations

Finite fields of prime characteristic (hence prime fields) are isomorphic to Z/pZ, the field of integers modulo p.
Usually addition and multiplication over Z/pZ are defined from the addition and multiplication over Z modulo p thus the classes of integers are 0, 1 ... p-1
Because we want to ease the arithmetic, we require our elements to be in the range -(p-1)/2..(p+1)/2 for any p > 2
For p = 2, there is no change, the operations are trivial anyway

So, how we do that?

## Addition in a cyclic group
For the addition, we use the group axioms and its cyclicity
i.e for a in the group 1+...+1 (p times) = 0
e.g for P=7, we map the elements 4, 5 and 6 respectively to -3, -2 and -1 because 1 + 1 + 1 + 1 = 4 and 1 + 1 + 1 = 3 
which implies 1 + ... + 1 (7 times) = 0 => 4 = -3

## Multiplication in a cyclic group
For the multiplication, we need a table (a,b) = a*b
We thus construct this table using the following axioms:

    a * 1 = 1 * a (note: this is enough for P = 2)
    a * -1 = -1 * a = -a
    
This is enough for P = 3

When P > 3, we compute the squares of 2..(p-1)/2 using the distributivity of + over *
e.g 2 * 2 (F5) = (1 + 1)(1 + 1) = 1 + 1 + 1 + 1 = -1
Then, the abelian group axioms gives the value of: 
    
    -a * -a and -a * a
    -a * -a = -1 * a * -1 * a  = a * a
    -a * a = -1 * a * a = - a*a

While we are at it, once a*a has been computed, we can compute the products of a and all elements > a
    
    a(a+1) = a*a + a
    a(a+1+1) = a*a + a + a
    and so on

The other products are computed using the commutativity and the associativity of an abelian group
There is certainly a way to optimize this but it does not seem to be so trivial.

For instance, the only product we need to know for F5 is 

    2 * 2 = -2 * -2 = -(-2 * 2)

For F7, there are : 

    2*2, 2*3 and 3*3

For F11, there are : 

    2*2, 2*3, 2*4, 2*5, 3*3, 3*4, 3*5, 4*4, 4*5 and 5*5

The gain would be the following:
- the table takes O(p^2) space because its a (p-1) * (p-1) matrix.
- if we note k = (p-1)/2 - 1 = (p-3)/2 then
    - the minimal number of products is k(k+1)/2
    - the order is the same: O(p^2)
    
The only difference is a reduction of the overall size by 4 as p grows


# Frobenius automorphism in a prime field
Taking the p-th root for a is finding an element b such as b^p = a where p is the field characteristic.
Because the multiplicative group of the field is cyclic of order p-1 if such an element b exists:

    b^(p-1) = 1 
    => b^p = b => a = b
    
In a prime field, the Frobenius automorphism is the identity

# Frobenius automorphism in a non-prime finite field
Taking the p-th root of an element a of a finite field of order p^d is finding an element b such as b^p = a. 
Contrary to the case of prime fields the answer is not that straightforward.

First recall the Freshman's Dream identity in finite fields of characteristic p:

    (x + y)^p = x^p + y^p
    => frobenius((x + y)^p) = (x + y) = frobenius(x^p + y^p)
    => frobenius(x^p + y^p) = frobenius(x^p) + frobenius(y^p) = x + y
    => frobenius is a linear map
    
Then, as any element of the field is a vector, let's write it a_0 + a_1j + ...
where the a_i are elements of the prime field and the basis is {1, j, ... j^d-1}
where j is the root of the polynomial that defines the field as a quotient

Thus, 

    frobenius(a) = frobenius(a_0 + a_1j + ... a_d-1 j^d-1)
    = frobenius(a_0) + frobenius(a_1j) + ... frobenius(a_d-1 j^d-1)
    = a_0 + a_1 frobenius(j) + ... a_d-1 frobenius(j^d-1)
    
This implies that we only need to pre-compute the p-th root of the basis elements which is the purpose of the function FiniteField._compute_frobenius_map

# Factorization algorithms

The imath module provides a full factorization algorithm for polynomials over finite fields. Although it is possible to define polynomials over Z or any extension of Z (gaussian integers for instance), factorization is only possible over finite fields.
The Factorization class embeds 3 methods of factorization: square free, distinct degree and equal degree. The full factorization is performed by the cantor_zassenhaus method.

## Square free factorization

This method accepts a monic polynomial (make it monic beforehand) and returns the square free part (that can be 1) and a list of all factors with multiplicity > 1.
The algorithm is based on 2 principles:

1) A polynomial with factors of multiplicity > 1 and its formal derivative share at least one common factor. We identify this factor by computing the GCD repeatedly.
2) If a polynomial over Fq has a factor of multiplicity r | q then its formal derivative is zero, we use the inverse Frobenius map to compute this factor.

## Distinct degree factorization

This method accepts a monic square free polynomial and returns a list of distinct degree factors. Each factor is also tagged with the degree of any equal degree factor. For instance, if a factor of degree 4 is the product of two distinct factors of degree 2, this factor will have its max_degree property set to 2. This value is then used by the equal degree factorization algorithm that completes the full factorization.
The algorithm exploits the following property:

    The polynomial X^q - X over Fq is the product of all irreducible polynomials of degree 1.
    
For instance, over F2 there are only 2 distinct irreducible polynomials of degree 1 : X and X+1 whose product is X^2 + X.
The same stand for F3. Its irreducible polynomials of degree 1 are : X, X+1, X-1 whose product is X^3-X.
Moreover, the property above can be extended to any arbitrary degree:

    The polynomial X^(q^n) - X over Fq is the product of all irreducible polynomials of degree <= n.
    
The algorithme computes repeatedly the GCD of the polynomial to factorize with the polynomial X^i - X to find the relevant factors.

## Equal degree factorization

This is the Cantor-Zassenhaus algorithm (see https://en.wikipedia.org/wiki/Cantor%E2%80%93Zassenhaus_algorithm for a thorough mathematical description of the algorithm).
The input is a monic, square free polynomial known to have k factors of equal degree d. It outputs a list of the factors. Be aware that this algorithm is non-deterministic and may not always provide an answer after several retries. In this case a runtime exception is raised. 
    
