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
    
This implies that we only need to pre-compute the p_th root of the basis elements which is the purpose of the function FiniteField._compute_frobenius_map
