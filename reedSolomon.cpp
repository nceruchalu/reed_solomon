/* 
 * -----------------------------------------------------------------------------
 * -----                           reedSolomon.cpp                         -----
 * -----                          REED-SOLOMON CODES                       -----
 * -----------------------------------------------------------------------------
 *
 * File Description:
 *   This is the implementation file for `reedSolomon` encoder/decoder class
 *
 * Assumptions:
 *   - Polynomials are written in vector form, or really an array where 
 *     zero-based index i represents the coefficient of x^i.
 *   - The factors read in the array slots are the Galois field elements in 
 *     decimal form, i.e. alpha^i values, where alpha is the primitive element
 *     of a primitive polynomial
 *
 * References:
 *   - http://downloads.bbc.co.uk/rd/pubs/whp/whp-pdf-files/WHP031.pdf
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision- encoder works
 *   Jun 01, 2011    Nnoduka Eruchalu    Added comments and wrote decoder
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */

#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;
#include "primitives.h"     // primitive elements
#include "reedSolomon.h"    // class declaration


///////////////////////////////////////////////////////////////////////////
// Initializers
///////////////////////////////////////////////////////////////////////////
/*
 * gen_gf()
 * Description: 
 *   Generate the Galois field and update representative arrays `alpha_to` and
 *   `index_of`
 *   
 *   ## Background on galois field
 *   - The Galois field consists of a set of elements (numbers). The elements 
 *     are based on a primitive element, denoted `alpha`, and take the values:
 *       0, alpha^0, alpha^1, alpha^2, ..., alpha^(N-1)
 *     to form a set of 2^m elements, where N = 2^m -1. The field is then known
 *     as GF(2^m)
 *  
 *   - In addition to the powers of alpha form, each field element can also be
 *     represented by a polynomial expresion of the form:
 *       a_(m-1)*x^(m-1) + ... + a_1*x + a_0
 *     where the coefficients a_(m-1) to a_0 take the values 0 or 1. Thus we can
 *     describe a field element using the binary number a_(m-1),...,a_1,a_0 and
 *     the 2^m field elements correspond to the 2^m combinations of the m-bit 
 *     number.
 *   
 *   - For example, in the Galois field with 16 elements (known as GF(16), m=4),
 *     the polynomial representation is:
 *       a_3*x^3 + a_2*x^2 + a_1*x^1 + a_0*x^0
 *     with a_3,a_2,a_1,a_0 corresponding to binary numbers 0000 to 1111.
 *
 *
 *   ## Constructing the Galois field
 *   - All non-zero elements of the Galois field can be constructed by using the
 *     fact that the primitive element alpha is a root of the primitive 
 *     polynomial, so that 
 *       p(alpha) = 0
 *   
 *   - Thus for GF(16) with the primitive polynomial [p(x) = x^4 + x + 1] in 
 *     primitives.h, we can say
 *       alpha^4 + alpha + 1 = 0
 *     or
 *       alpha^4 = alpha + 1 (+ and - are the same in a GF; see refernce)
 *   
 *    - Multiplying by alpha at each stahe and using [alpha+1] to substitute for
 *      alpha^4 and adding the resulting terms can be used to obtain the
 *      complete field as shown in the table below:
 *      +------------+---------------------------+-------------+--------------+
 *      | index form | polynomial form           | binary form | decimal form |
 *      +------------+---------------------------+-------------+--------------+
 *      | 0          | 0                         | 0000        | 0            |
 *      | alpha^0    | 1                         | 0001        | 1            |
 *      | alpha^1    | alpha                     | 0010        | 2            |
 *      | alpha^2    | alpha^2                   | 0100        | 4            |
 *      | alpha^3    | alpha^3                   | 1000        | 8            |
 *      | alpha^4    | alpha+1                   | 0011        | 3            |
 *      | alpha^5    | alpha^2 + alpha           | 0110        | 6            |
 *      | alpha^6    | alpha^3 + alpha^2         | 1100        | 12           |
 *      | alpha^7    | alpha^3 + alpha + 1       | 1011        | 11           |
 *      | alpha^8    | alpha^2 + 1               | 0101        | 5            |
 *      | alpha^9    | alpha^3 + alpha           | 1010        | 10           |
 *      | alpha^10   | alpha^2 + alpha + 1       | 0111        | 7            |
 *      | alpha^11   | alpha^3 + alpha^2 + alpha | 1110        | 14           |
 *      | alpha^12   | alpha^3 +alpha^2 +alpha +1| 1111        | 15           |
 *      | alpha^13   | alpha^3 + alpha^2 + 1     | 1101        | 13           |
 *      | alpha^14   | alpha^3 + 1               | 1001        | 9            |
 *      +------------+---------------------------+-------------+--------------+
 *      If we continue this process beyond alpha^14, we see alpha^15 = alpha^0,
 *      alpha^16 = alpha^1,.... so it's a repeating sequence.
 * 
 *
 * Assumptions:
 *   alpha_to and index_of arrays have had their memory allocated
 *  
 * Arguments:
 *   None
 *
 * Return:
 *   None
 * 
 * Operation:
 *   - This Galois field is represented by the arrays alpha_to and index_of [See
 *     `reedSolomon.h` for more details]. The goal is to populate the array 
 *     elements with the binary/decimal representation of the polynomial form of
 *     alpha^0 to alpha^(2^m -2). 
 *   - These bounds of 0 to (2^m -2) were chosen because:
 *     + We dont worry about GF element 0; that is always represented as 0
 *     + alpha^(2^m -1) = alpha^0, alpha(2^m) = alpha^1, so at this point it's 
 *       a repeating sequence.
 *   - At the end of it all, it is important to note that for no value i is
 *     alpha^i == 0, so index_of[0] isnt a valid field. Indicate this with -1.
 *
 * References:
 *   - http://downloads.bbc.co.uk/rd/pubs/whp/whp-pdf-files/WHP031.pdf
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::gen_gf()
{   
    int i, mask;
    
    mask = 1 ;
    alpha_to[m] = 0 ;
    
    for (i=0; i<m; i++)
    { 
        alpha_to[i] = mask ;
        index_of[alpha_to[i]] = i;
        
        if (p_x[i]!=0)
            alpha_to[m] ^= mask;
        
        mask <<= 1;
    }
    index_of[alpha_to[m]] = m;
    mask >>= 1 ;
    
    for (i=m+1; i<n; i++)
    {
        if (alpha_to[i-1] >= mask)
            alpha_to[i] = alpha_to[m] ^ ((alpha_to[i-1]^mask)<<1) ;
        else
            alpha_to[i] = alpha_to[i-1]<<1 ;
       
        index_of[alpha_to[i]] = i ;
    }
    
    index_of[0] = -1 ;
  
}

/*
 * gen_prim_poly()
 * Description: 
 *   Generate the primitive polynomial by setting the primitive array, p_x 
 *   [or p(x) really] with an appropriate primitive polynomial of order m.
 *   This info is stored in string form in the array primitive[] (defined in
 *   primitives.h)
 *
 *  ## Background on primitive polynomial
 *  - A Primitive polynomial is a polynomial of degree m which is irreducible. 
 *     It forms part of the process of multiplying two field elements together. 
 *     For a Galois field of a particular size, there is sometimes a choice of
 *     suitable polynomials. Using a different primitive polynomial from that
 *     specified will produce incorrect results.
 *
 *   - As an example, for GF(16), the polynomial [p(x) = x^4 + x + 1] is
 *     irreducible and is thus a suitable primitive polynomial.
 *
 * Assumptions:
 *   - the primitive[] array has the primitive array of order m at index (m-1)
 *   - p(x) has already been allocated the appropriate memory size
 *  
 * Arguments:
 *   None
 *
 * Return:
 *   None
 *
 * Operation:
 *   - Set all coefficients of p(x) to 0.
 *   - read in the primitive polynomial string from primitives.h and update the
 *     non-zero coefficients based on that string.
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::gen_prim_poly() // create primitive polynomial p(x)
{
    string prim_poly_str = primitive[m-1];   // get string version of p(x)
    string token;                            // each character in string

    for(int i=0; i < size_p; i++)            // set all coeffs. of p(x) to 0
        p_x[i] = 0;
    
    stringstream ss(prim_poly_str);          // string contains only powers of
    while(getline(ss, token, ' '))           // x with coefficient = 1. Get them
        p_x[atoi(token.c_str())] = 1;        // and update p(x) accordingly.
}

/*
 * gen_g_poly()
 * Description: 
 *   Create the generator polynomial which is really just
 *   g(x) = (x + a)(x + a^2)...(x+a^2t)
 *   this makes g(x) a polynomial of order 2t (size_g is 2t+1 coz of constant) 
 *
 *   ## Background on (code) generator polynomial
 *   - An (n, k) Reed-Solomon code is constructed by forming the code generator
 *     polynomial g(x), consisting of n-k=2t factors, the roots of which are
 *     consecutive elements of the Galois field. Choosing consecutive elements
 *     ensures that the distance properties of the code are maximised.
 *     Thus the code generator polynomial takes the form:
 *       g(x) = (x - alpha)(x - alpha^2)...(x - alpha^2t)
 *     We are working with mod 2, so - and + are interchangeable:
 *       g(x) = (x + alpha)(x + alpha^2)...(x + alpha^2t)
 *
 *   - This makes g(x) of order 2t:
 *       g(x) = g_0 +  g_1*x + g_2*x^2 + ... + g_(2t-1)*x^(2t-1) + x^2t
 *     so the array representation of g(x) will need to account for 2t coeffs
 *     + 1 constant, so size of g(x) array ,`size_g`, is 2t + 1
 *  
 *
 * Assumptions:
 *   g(x) has already been allocated the appropriate memory size of 2t + 1
 *  
 * Arguments:
 *   None
 *
 * Return:
 *   None
 * 
 * Operation:
 *   - clear out g(x)
 *   - Start off with g(x) = (x+alpha)
 *   - Now consecutively multiply g(x) with (x+alpha^i) where i goes from 2 to
 *     2t (Note 2t == size_g-1)
 *     + So g(x) = x*g(x) + B*g(x), where B = alpha^i
 *     + x*g(x) and B*g(x) will be computed separately.
 *     + x*g(x): 
 *       ++ will never have a constant term so set coefficient of x^0 to 0
 *       ++ all other coefficients in g(x) are bumped up to one high power of x
 *     + B*g(x), where B=alpha^i:
 *       ++ All non-zero coefficients in g(x) are multiplied using the fact that
 *          we are multiplying powers of alpha. So each coefficient is converted
 *          to its appropriate alpha^j. The coefficient scaling then becomes:
 *            alpha^[index_of(coefficient)] * alpha^i
 *              == alpha^((index_of(coefficient) + i) % n)
 *     + Finally perform x*g(x) + B*g(x) keeping in mind that addition is binary
 *       XOR
 *   
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::gen_g_poly() // create the generator polynomial g(x)
{
    int i, j;
    int *xg_x = new int[size_g]; 
    int *Bg_x = new int[size_g];

    // clear out g(x) by initializing all slots to 0
    for(i =0; i < size_g; i++)
        g_x[i] = 0;

    // start off with g(x) = (x + alpha)
    g_x[0] = alpha_to[1]; 
    g_x[1] = 1;
    
    
    // consecutively multiply g(x) with (x+alpha^i) where i goes from 2 to 2t
    // take advantage of fact that 2t == size_g-1
    for(i = 2; i < size_g; i++)
    {   
        // g(x)  = g(x)*(x+B) = x*g(x) + B*g(x),  where B = alpha^i
        
        // first compute x*g(x) 
        xg_x[0] = 0;                       // xg(x) doesn't have a constant term
        for(int k = size_g-1; k >= 1; k--) // all other coefficients are bumped
            xg_x[k] = g_x[k-1];            // up one level
        
        // then compute B*g(x)
        for(int k = 0; k < size_g; k++)   
            Bg_x[k] = (g_x[k] ?  alpha_to[(index_of[g_x[k]] + i)%n] : 0);
       
        // finally add both results: g(x) = x*g(x) + B*g(x) using binary XOR
        // between corresponding coefficients.
        for(j = size_g - 1; j >= 0; j--) 
            g_x[j] = xg_x[j] ^ Bg_x[j];  
    }
    
    // deallocate temporary memory
    delete[] xg_x;
    delete[] Bg_x;
}


/*
 * update_params()
 * Description: 
 *     Initalizae the instance variables and arrays representing multiple
 *     polynomials and the Galois Field(2^m) elements.
 *     
 *     - This takes into account the assumptions of gen_gf(), gen_prim_poly(),
 *       gen_g_poly() and as such initializes p(x), g(x), alpha_to[], index_of[]
 *
 *     - This should be called by the constructors only, as it allocates memory
 *       without attempting to deallocate previously allocated memory. It set
 *
 * Assumptions
 *     -Assumes m and t have already been setup
 *  
 * Arguments:
 *     None
 *
 * Return:
 *     None
 *
 * Operation:
 *   - the constructor only passes the m and t values. so n and k have to be
 *     computed using their formulas: (n= 2^m -1) and (k = n -2t)
 *   - start object in a state of no transmission error
 *   - determine polynomial array sizes for g(x), p(x), s(x) by adding 1 to the
 *     polynomial's order. The +1 accounts for the constant term
 *   - When done allocating memory call functions to update p(x),
 *     generate the galois field by updating alpha_to[] and index_of[], then
 *     finally creates g(x)
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::update_params()
{
    n = pow(2,m) - 1;         // compute n and k
    k = n- 2*t;
    detect_error = false;     // start in state of no transmission error
    
                              // polynomial's size = its order +1 (for constant)
    size_g = 2*t + 1;         // g(x) is order 2t; see gen_g_poly();
    size_p = m +1;            // p(x) is order m; see gen_prim_poly()
    size_s = 2*t;             // s(x) is order (2t-1): see get_syndromes()
    g_x = new int[size_g];
    p_x = new int[size_p];
    s_x = new int[size_s];
    m_x = new int[n];         // should be of size k, but n makes algo easier.
    c_x = new int[n];
    rc_x = new int[n];
    dc_x = new int[n];
    alpha_to = new int[n+1];  // allocate memory for galois field representation
    index_of = new int[n+1];  // arrays
    
    gen_prim_poly();   // generate primitive polynomial
    gen_gf();          // update alpha_to and index_to -- they represent the GF
    gen_g_poly();      // now can create generator polynomial
}


/*
 * reedSolomon() -- default constructor
 * Description: 
 *   Default constructor that initializes object with m=4 and t=3
 *
 * Assumptions:
 *   None
 *  
 * Arguments:
 *   None
 *
 * Return:
 *   None
 *
 * Operation:
 *   Create the reedSolomon object with parameters m=4, t=3 then call
 *   update_params() to setup other instance variables
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
reedSolomon::reedSolomon() // default constructor
{
    m = 4;
    t = 3;
    update_params();
}


/*
 * reedSolomon(int mm, int tt) -- more useful constructor with 2 arguments
 * Description: 
 *   Constructor that initializes object with m and t value arguments.
 *
 * Assumptions:
 *   None
 *  
 * Arguments:
 *   mm - integer value representing the reedSolomon object's m value
 *   tt - integer value representing the reedSolomon object's t value
 *
 * Return:
 *   None
 *
 * Operation:
 *   Create the reedSolomon object with parameters m=mm, t=tt then call
 *   update_params() to setup other instance variables.
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
reedSolomon::reedSolomon(int mm, int tt) // more useful constructor
{
    m = mm;
    t = tt;
    update_params();
}


/*
 * ~reedSolomon() -- destructor
 * Description: 
 *   Delete any memory allocated by the constructor
 *
 * Assumptions:
 *   None
 *  
 * Arguments:
 *   None
 *
 * Return:
 *   None
 *
 * Operation:
 *   Deallocate memory of instance's arrays.
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
reedSolomon::~reedSolomon()  // destructor
{
    delete [] g_x;
    delete [] p_x;
    delete [] m_x;
    delete [] c_x;
    delete [] rc_x;
    delete [] dc_x;
    delete [] alpha_to;
    delete [] index_of;
}



///////////////////////////////////////////////////////////////////////////
// Galios field arithmetic
///////////////////////////////////////////////////////////////////////////
/*
 * doAdd(int * & cc, int *m1, int *m2, int size_n)  
 * Description: 
 *   cc(x) = m1(x) + m2(x)
 *  
 * Arguments:
 *   - cc = cc(x); 
 *   - m1 = m1(x);
 *   - m2 = m2(x)
 *   - size_n = size of polynomial arrays == (order of cc(x) +1)
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   - All polynomial arrays have already been setup and allocated memory.
 *   - all polynomials have at least size_n elements. Not error-checking this!
 *
 * Operation:
 *   Bitwise XOR of corresponding terms as that's what Galois field addition is.
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::doAdd(int * & cc, int *m1, int *m2, int size_n) 
{
    for(int i=0; i<size_n; i++)
        cc[i] = m1[i] ^ m2[i];
}


/*
 * doSub(int * &m1, int *m2, int size_n)
 * Description: 
 *   m1(x) -= m2(x)   or   m1(x) += m2(x)
 *  
 * Arguments:
 *   - m1 = m1(x); 
 *   - m2 = m2(x)
 *   - size_n = size of polynomial arrays == (order of m1(x) +1)
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   - All polynomial arrays have already been setup and allocated memory.
 *   - all polynomials have at least size_n elements. Not error-checking this!
 *
 * Operation:
 *   Bitwise XOR of corresponding terms per Galois field subtraction rules.
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::doSub(int * &m1, int * & m2, int size_n)
{
    for(int i=0; i<size_n; i++)  // m1(x)-= m2(x) <===> m1(x) += m2(x)
        m1[i] ^= m2[i];
}


/*
 * mul(int * & res1, int *g1, int fact, int w, int size_res1, int size_g1)
 * Description: 
 *   res1(x) = g1(x) * (alpha^fact) * x^w
 *  
 * Arguments:
 *   - res1      = res1(x); g1 = g1(x)
 *   - fact      = index form of GF element (power of alpha) to multiply g1(x)
 *   - w         = power of x to multiple g(x) by
 *   - size_res1 = size of array res1
 *   - size_g1   = order of g1(x) +1
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   - All polynomial arrays have already been setup and allocated memory.
 *   - size_res1 >= size_g1+w    NO ERROR CHECKING
 *
 * Operation:
 *   - When this operation is done, the order of res1(x) will be increased by w
 *     to size_g1 + w -1. So zero out all coefficients for variables with powers
 *     higher than that: [size_g1+w, size_res1-1] 
 *   - Since there is a multiplication by x^w, each term's power of x will be at
 *     least w, so zero out all coefficients for variables of lower powers
 *   - For all other terms with powers of x in the range [w, size_g1+w-1]
 *     perform multiplacation by adding powers of the index forms of the GF
 *     elements and keeping it modulo (2^m -1)
 *
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::mul(int * & res1, int *g1, int fact, int w, int size_res1,
                      int size_g1)
{
    int i;
    // zero out coefficients of terms with higher than possible powers of x
    for(i = size_res1-1; i >= size_g1+w; i--)
        res1[i] = 0;
    
    //  
    for(i = size_g1-1+w; i >= w; i--)
        res1[i]= (g1[i-w] ? alpha_to[(index_of[g1[i-w]] + fact)%n]: 0);
    
    // zero out coefficients of terms with lower than possible powers of x
    for (i = w-1; i >= 0; i--)
        res1[i] = 0;
}


/*
 * toPower(int * & res, int * mm, int size_mm, int w) 
 * Description:
 *   res(x) = m(x) * x^w
 *  
 * Arguments:
 *   - res = res(x)
 *   - mm = mm(x)
 *   - size_mm = size of array mm = order of mm(x) + 1 -- remember constant term
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   - res(x) and mm(x) have already been setup and allocated memory.
 *   - w >= 0
 *
 * Operation:
 *   - move coefficients to terms of higher power by a factor of w
 *   - Minimum power of res(x) will be w, so zero all coefficients for terms of
 *     lower powers
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments and stopped
 *                                       allocating result's memory
 */
void reedSolomon::toPower(int * & res, int * mm, int size_mm, int w)
{
    int i;
    
    // move coefficients to terms of higher powers by a factor of w
    for(i = size_mm-1; i >= w; i--)
        res[i] = mm[i-w];
    
    // zero coefficients for terms with too low powers
    for (i = w-1; i >= 0; i--)
        res[i] = 0;
}


/*
 * doMul(int * & res_x, int * a_x, int * b_x,  int size_res, int size_a, 
 *       int size_b)
 *
 * Description:
 *   res(x) = a(x)*b(x)
 *  
 * Arguments:
 *   - res_x = res(x); a_x = a(x); b_x = b(x)
 *   - size_a = size of array a_x
 *   - size_b = size of array b_x
 *   - size_res = size of array res_x
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   - all arrays setup and have memory allocated
 *   - size_res >= order_of(a(x)) + order_of(b(x))+1
 *
 * Operation:
 *   - get orders of a(x) and b(x). Note that this isn't simply (size_a-1) and 
 *     (size_b-1) as those values will tell max. possible orders not actuals.
 *   - initialize res(x) by clearing out all coefficiients. This will now be
 *     used as an accumulator
 *   - now do actual multiplication by looping through a(x)'s terms.
 *     + for each non-zero term, get the value:
 *       temp_res(x) = b(x)* [term's coefficient] * x^[term's power]
 *     + accumulate this temporary result into res(x)
 *
 * Revision History
 *   Jun 01, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::doMul(int * & res_x, int * a_x, int * b_x,  int size_res,
                        int size_a, int size_b)
{
    // get orders of a(x) and b(x)
    int order_a = get_order(a_x, size_a);
    int order_b = get_order(b_x, size_b);
    
    // clear out result polynomial to prepare it for accumulation
    for(int i = 0; i < size_res; i++) res_x[i] = 0;
    // create temporary result polynomial
    int *temp_res = new int[size_res];
    
    for(int i=0; i <= order_a; i++) // looping through a(x)'s terms
    {
        if(a_x[i])
        {
            // multiply this non-zero term with b(x) and add to res(x)
            int fact = index_of[a_x[i]];
            mul(temp_res, b_x, fact, i, size_res, order_b+1);
            doSub(res_x, temp_res, size_res);  // same thing as an add!
        }
    }
    
    // deallocate temporary result
    delete [] temp_res;
}


/*
 * doDiv(int * & quot_x, int * & rmd_x, int * a_x, int * b_x, int size)
 *
 * Description:
 *   quot(x) = a(x) / b(x)
 *   rmd(x)  = a(x) % b(x)
 *  
 * Arguments:
 *   - quot_x = quot(x); 
 *   - a_x = a(x): dividend (numerator)
 *   - b_x = b(x): divisor (denominator)
 *   - size = size of all arrays
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   - all arrays are setup and have memory allocated
 *   - all arrays are the same size
 *
 * Operation:
 *   - initialize quotient and remainder, quot(x) = 0, rmd(x) = a(x)
 *   - get highest powers (actual orders) of a(x) and b(x).
 *   - if order[a(x)] > order[b(x)] then done, else proceed
 *   - while order[b(x)] stil <= order[rmd(x)] keep performing this loop:
 *     + get coefficients of highest order terms of rmd(x) and b(x). Get these
 *       in GF index forms.
 *     + Divide both coefficients, by subtracting the powers (of alpha).
 *       ++ To ensure this GF index form stays valid, add `n` and take mod. `n`
 *       ++ save this resulting value's index as `fact`
 *     + update the appropriate term of quot(x) with the decimal form of
 *       alpha^`fact`
 *       ++ the appropriate term here is that whose variable's power is:
 *          order[rmd(x)] - order[b(x)]]
 *     + subtract the product of the added quotient term and b(x) from remainder
 *       ++ temp(x) = b(x) * alpha^`fact` * x^(order[rmd(x)] - order[b(x)]])
 *       ++ rmd(x) -= temp(x)
 
 *
 * Revision History
 *   Jun 01, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::doDiv(int * & quot_x, int * & rmd_x, int * a_x, int * b_x, 
                        int size)
{
    for(int i = 0; i < size; i++)  // initialize quotient and remainder;
    {
        quot_x[i] = 0;
        rmd_x[i] = a_x[i];
    }
    
    int order_b = get_order(b_x, size);
    int order_a = get_order(a_x, size);
    
    if(order_a < order_b)  // if order(a(x)) < order(b(x)), done
        return;
    
    int coeff_b = index_of[b_x[order_b]]; // coefficient of highest order term
                                          // of b_x = alpha^coeff_b
    int * temp = new int[size];
    while(order_b <= get_order(rmd_x, size)) // loop through dividend
    {
        // get coefficients of highest order terms of rmd(x) and b(x)
        // and divide both coefficients
        int w = get_order(rmd_x, size); // clearly mod[w] != 0
        int fact = index_of[rmd_x[w]];  // so the value here is a^fact
        fact += n - coeff_b;           
        fact %= n;    
        
        // add result of coefficient division as coefficient of a new quotient 
        // term. This new term's variable's power will be:
        //   order[rmd(x)] - order[b(x)]]
        quot_x[w-order_b] = alpha_to[fact];
        
        // subtract the product of the added quot. term and b(x) from remainder
        mul(temp, b_x, fact, w-order_b, size, order_b+1);
        doSub(rmd_x, temp, size); 
    }

    delete [] temp;
}



///////////////////////////////////////////////////////////////////////////
// Miscellaneous routines
///////////////////////////////////////////////////////////////////////////
/*
 * print_poly(int *pp_x, int size_pp)
 * Description: 
 *   print the polynomial in the Galois Field GF(2^m)
 *
 * Assumptions:
 *   - for each polynomial, zero-based index i represents the factor
 *     corresponding to x^i
 *   - the coefficients stored in the polynomial array slots are the decimal
 *     forms of the Galois field elements (alpha^j) values
 *  
 * Arguments:
 *   pp_x = integer array representing polynomial to be printed pp(x)
 *   size_pp = size of pp(x) == (order of pp(x) +1)
 *
 * Return:
 *   None
 *
 * Operation:
 *   Loop through all coefficients of the polynomial:
 *   - Only consider terms with non-zero coefficients.
 *   - Get index form of GF element, save that as `pow`
 *   - Print all coefficients in a^pow form, except special case of pow==1. 
 *     Print that as "a"
 *   - Print "x^(power of variables) + " for all variables.
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::print_poly(int *pp_x, int size_pp)
{
    for(int i = size_pp-1; i>= 0; i-- )    // loop through all terms of pp(x)
    {
        if(pp_x[i])                        // only consider non-zero coeffs.
        {
            int pow = index_of[pp_x[i]];   // this coeff. is really alpha^pow
            if(pow == 0) cout << "1";      // so print it as such, but note
            else if(pow == 1)  cout << "a";// that "a" is neater than "a^1"
            else if (pow > 1) cout << "a^" << pow;

            if(i > 0) cout << "x^" << i<<" + ";// print x for non-constant terms
        }
    }
    cout << endl;
}


/*
 * print_params()
 * Description: 
 *   print the properties of the reed Solomon object, including its
 *   many polynomials and the Galois Field info
 *
 * Assumptions:
 *   None
 *  
 * Arguments:
 *   None
 *
 * Return:
 *   None
 *
 * Operation:
 *   - print out m,t,n,k values
 *   - print out alpha_to and index_of values in tabular form
 *   - print polynomials
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::print_params()
{
    // print out m,t,n,k values
    cout <<  "m: " << m
         << "\tt: " << t
         << "\tn: " << n
         << "\tk: " << k
         << endl;
    
    // print out table resulting from alpha_to[]
    cout << "field elements for GF(" << (n+1);
    cout << "): index form -> decimal form" << endl;
    cout << "index form\tdecimal form" << endl;
    for(int i=0; i<n; i++)
      cout << "a^" << i <<"\t\t" << alpha_to[i] << endl;
    
    // print out table resulting from index_of[]
    cout << endl;
    cout << "field elements for GF(" << (n+1);
    cout << "): decimal form -> index form" << endl;
    cout << "decimal form\tindex form" << endl;
    for(int i=1; i<=n; i++)
      cout << i <<"\t\t" << "a^" << index_of[i] << endl;
    
    // print polynomials
    cout << "p(x):  ";  print_poly(p_x, size_p);
    cout << "g(x):  ";  print_poly(g_x, size_g);
    cout << "m(x):  ";  print_poly(m_x, n);
    cout << "c(x):  ";  print_poly(c_x, n);
    cout << "rc(x): ";  print_poly(rc_x, n);
    cout << "s(x):  ";  print_poly(s_x, size_s);
    cout << "dc(x):  "; print_poly(dc_x, n);
}


/*
 * get_order(int *pp_x, int size_pp)
 * Description: 
 *   Get the order of the given polynomial pp(x). Note that the order
 *   of a polynomial is the highest power of x with a non-zero coefficient
 *  
 * Arguments:
 *   pp_x = integer array representing polynomial of interest
 *   size_pp = size of pp(x) == (order of pp(x) +1)
 *
 * Return:
 *   (int) order of polynomial
 *
 * Assumptions:
 *   - Polynomials are written in vector form, or really an array where 
 *     zero-based index i represents the coefficient of x^i.
 *
 * Operation:
 *   loop through polynomial's terms starting from highest possible power term. 
 *    Stop on coming across a non-zero coefficient and return the index. 
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
int reedSolomon::get_order(int *pp_x, int size_pp)
{
    int i;
    for(i = size_pp-1; i>= 0; i--)
    {
        if(pp_x[i] != 0) break;
    }
    return i;
}


/*
 * copy_arr(int * & dst,int * src, int size)
 * Description:
 *   dst[] = src; Copy contents of src INTO dst
 *  
 * Arguments:
 *   dst = destination array
 *   src = source array
 *   size = size of both arrays
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   Both arrays are of same size
 *
 * Operation:
 *   Loop through both arrays and assign elements in dst to the same value as
 *   the corresponding elements in src.
 *
 * Revision History
 *   Jun 01, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::copy_arr(int * & dst,int * src, int size)
{
    for(int i =0; i< size; i++)
        dst[i] = src[i];
}



///////////////////////////////////////////////////////////////////////////
// Reed-Solomon encoding/decoding methods
///////////////////////////////////////////////////////////////////////////
/*
 * gen_rand_msg()
 * Description:
 *   Randomly populates message vector m(x) of order k-1 so that:
 *   a random message vector = (m_0, m_1, ..., m_(k-1)) is represented as a
 *   message polynomial, m(x) = m_0 + m_1*x + ... + m_(k-1)*x^(k-1)
 *
 * Arguments:
 *   None
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   - m(x) is an array of size n, but only first k (n>k) slots are used for
 *    polynomial m(x)
 * 
 * Operation:
 *   - loop through first k slots and randomly generate decimal form GF(2^m)
 *     elements
 *   - zero out
 *  
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::gen_rand_msg()
{
    // populate first k slots with random GF(2^m) elements
    int i, mod = pow(2,m);
    for(i =0; i<k; i++)
        m_x[i] = rand() % mod;
    // zero-out unused slots.
    for(i = k; i < n; i++)
        m_x[i] = 0;
}


/*
 * encode()
 * Description:
 *   Do a reed solomon decoding of message vector m(x) using this formula:
 *   c(x) = m(x) * x^(n-k)   +   [m(x) * x^(n-k)] % g(x)
 *   
 *   Clear out detect error flag
 *   
 *   Note r(x) = [m(x) * x^(n-k)] % g(x) ensures that encoded message 
 *   polynomial, c(x) is always divisible by the generator polynomial, g(x),
 *   without a remainder.
 *
 * Arguments:
 *   None
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   None
 *
 * Operation:
 *   - m_toPow = m(x)*x^(n-k) [using toPower()]
 *   - m_withMod = [m(x) * x^(n-k)] % g(x) = [m_toPow % g(x)]
 *     ++ This requires use of doDiv()
 *     ++ As such there will be need for temp array to hold unused quotient
 *     ++ Also g(x)'s polynomial array of size `size_g`, but doDiv() requires
 *        all arrays be the same size. For this reason g(x) will be copied into
 *        a larger array of size n. The unused (higher-order) polynomial terms
 *        will have their coefficients cleared.
 *   - c(x) = m_toPow(x) + m_withMod(x)
 *   - since codeword was just created, no error can be detected yet
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::encode()
{
    int *m_withMod = new int[n];     // m_withMod = [m(x) * x^(n-k)] % g(x)
    int *m_toPow = new int[n];       // m_toPow = m(x)*x^(n-k)
    int *m_withDiv = new int[n];     // temp array to hold doDiv()'s quotient
    int *mg_x = new int[n];          // need a version of g(x) of size n > k
    
    toPower(m_toPow, m_x, n, n-k);   
    
    // generae new version of g(x) as required by doDiv()
    for(int i =0; i<size_g; i++) 
        mg_x[i] = g_x[i];
    for(int i=size_g; i<n; i++)
        mg_x[i] = 0;
    
    doDiv(m_withDiv, m_withMod, m_toPow, mg_x, n);
    
    // finally create code word.
    doAdd(c_x, m_toPow, m_withMod, n);
    
    // no error detected yet
    detect_error = false;
    
    delete [] m_toPow;               // cleanup
    delete [] m_withMod;
    delete [] m_withDiv;
    delete [] mg_x;
}


/*
 * sim_channel()
 * Description:
 *   Pass encoded codeword c(x) through the channel and generate the
 *   received codeword rc(x) with at most t errors
 *  
 * Arguments:
 *   None
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   rc(x) has already been allocated memory, and c(x) is setup
 *
 * Operation:
 *   - first set rc(x) = c(x)
 *   - then randomly create up to t errors. Note that the goal isn't to always
 *     create exactly t errors, but to make t attempts at creating errors. It
 *     could happen that the error occurs on the same term multiple times.
 *
 * Revision History
 *   Jun 01, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::sim_channel()
{
    for(int i=0; i<n;i++)  // first copy entire codeword
        rc_x[i] = c_x[i];
    
    for(int i=0; i<t; i++) // then create up to t errors
        rc_x[rand()%n] ^= alpha_to[rand()%n];
}

/*
 * sim_channel(double Ps)
 * Description:
 *    pass encoded codeword c(x) through the channel and generate the
 *    received codeword rc(x) using the following channel model:
 *      The codeword elements (c_0, c_1, ..., c_(n-1)) are transmitted over a 
 *      channel where symbol errors occur independently with probability Ps, 
 *      i.e. each received symbol rc_i = c_i + e_i, where rc_i, c_i, e_i [error]
 *      are IN GF(2^m) and:
 *                       |1 - Ps      if a == 0
 *      Prob(e_k = a) = -|
 *                       |Ps/(2^m -1) if a != 0
 *  
 * Arguments:
 *   Ps - symbol error probability
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   rc(x) has already been allocated memory, and c(x) is setup
 *
 * Operation:
 *   - create each term of rc(x) from c(x) based on the following algo:
 *     + if a randomly generated percentage < Ps then:
 *        r_i = c_i + e_k, where e_k is a random element in GF(2^m)
 *     + else:
 *        r_i = c_i
 *
 * Revision History
 *   Jun 02, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::sim_channel(double Ps)
{
    for(int i=0; i<n; i++) // create errors with probability Ps
    {       
        if (((double) rand()/RAND_MAX) < Ps)
            rc_x[i] = c_x[i]^alpha_to[rand()%n];
        else
            rc_x[i] = c_x[i];
    }
}


/*
 * get_syndromes()
 * Description:
 *   Get syndrome polynomial from received codeword rc(x)
 *   then set the detect_error flag if s(x) != 0
 *
 *   # Background on calculating syndromes
 *   - As explained in encode(), c(x) is always perfectly divisble by g(x).
 *     Therefore the first step in the decoding process is to divide the
 *     received polynomial by each of the factors (x+ alpha^i) of the generator
 *     polynomial, g(x). This produces a quotient and a remainder that is:
 *       rc(x)/(x + alpha^i) = Q_i(x) + S_i/(x + alpha^i) for 1 <= 1 <= 2t
 *
 *   - The emainders S_i resulting from these divisions are known as syndromes
 *     and can be written as S_1,...,S_2t
 *  
 *   - Rearranging the equation above produces:
 *       S_i = Q_i(x)*(x + alpha^i) + rc(x)
 *
 *   - When x = alpha^i this reduces to:
 *       S_i = rc(alpha^i)
 *           = rc_0 + rc_1*alpha^i + ... + rc_(n-1)*[(alpha^i)^(n-1)]
 *       where the coefficients rc_0,rc_1,...rc_(n-1) are the symbols of the
 *       received code word, rc(x).
 *
 *   - This means that each of the syndrom values can be obtained by 
 *     substituting x = alpha^i in the received polynomial, rc(x), as the
 *     alternative to the division of rc(x) by (x + alpha^i) to form a remainder
 *  
 *   - Substituting in the equation: rc(x) = c(x) + Error(x)
 *       rc(alpha^i) = c(alpha^i) + E(alpha^i)
 *     in which c(alpha^i) = 0 because (x+alpha^i) is a factor of g(x), which
 *     is also a factor of c(x). So:
 *       rc(alpha^i) = E(alpha^i) = Si.
 *
 *   - This means that the syndrome values are only dependent on the error
 *     pattern and are not affected by the data values. Also, when no errors
 *     have occured, all the syndrome values are zero.
 *
 *   - The syndrome polynomial is used by the Euclidean algorithm for finding
 *     the coefficients of the error locator polynomial:
 *       s(x) = S_1 + S_2*x + ... + S_2t*x^(2t-1)
 *     Thus s(x) is of order (2t-1) and so its array representation's size,
 *     `size_s`, must be 2t.
 *
 * Arguments:
 *   None
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   - s(x) has already been allocated memory
 *   - rc(x) is setup
 *
 * Operation:
 *   - Generate all syndrome values for use in populating the syndrome
 *     polynomial, s(x) by looping over 1 <= i <= 2t (== size_s)
 *     + S_i = rc(alpha^i), which is implemented by:
 *       ++ looping over all n terms of rc(x) and converting the decimal
 *          form coefficients into index form GF(2^m) elements.
 *          The variables of these the corresponding terms will be converted to
 *          x^(power of term) == alpha^i*(power of term).
 *          Now that both coefficients and variables are in GF(2^m) index form,
 *          they can be easily multiplied and cummulatively added into S_i.
 *     + sc_(i-1) = S_i, to fit into how we define polynomial array reprs.
 *   - If any syndrome value is non-zero set detect_error to true.
 *
 * References:
 *   - http://downloads.bbc.co.uk/rd/pubs/whp/whp-pdf-files/WHP031.pdf 
 *
 * Revision History
 *   Jun 01, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::get_syndromes()
{
    for(int i=1; i <= size_s; i++)  // i is index into alpha^i values
    {                                // so its +1 of s_x index
        s_x[i-1] = rc_x[0];
        
        for(int j=1; j<n; j++)
        {
            if(rc_x[j])  //cummulative additions when there is a factor
                s_x[i-1] ^= alpha_to[(index_of[rc_x[j]]+i*j)%n];
        }
        
        if(s_x[i-1]) detect_error = true;
    }
}


/*
 * euclid(int *a_x, int *b_x, int * & r_xi, int * &t_xi)
 *
 * Description:
 *   Run the extended Euclidean Algorithm for fiding the highest common factor
 *   of two numbers. Use this for finding the coefficients of the error
 *   locator polynomial, lambda(x), with the special end condition of:
 *     order(r_xi) < t;
 * 
 * Arguments:
 *   - a_x = a(x); 
 *   - b_x = b(x);
 *   - r_xi = r(x)_i;  
 *   - t_xi = t(x)_i
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   -a(x) and b(x) are already setup
 *   -all arrays are setup and have memory allocated
 *   -all arrays except b(x) are the same size of (2t+1)
 *   -b(x) is of size 2t
 *
 * Operation:
 *  Below is the Euclidean algorithm in pseudo-code:
 *  Initialization: r_(-1)(x) = a(x);  r_0(x) = b(x); 
 *                  s_(-1)(x) = 1;     s_0(x) = 0
 *                  t_(-1)(x) = 0;     t_0(x) = 1;
 *  while(deg(r_i(x)) >= t)
 *  {
 *     // compute quotient[q_i(x)] and remainder [r_i(x)]
 *     q_i(x) = r_(i-2)(x) / r_(i-1)(x)     
 *     r_i(x) = r_(i-2)(x) % r_(i-1)(x) = r_(i-2)(x) - q_i(x)*r_(i-1)(x)
 *     s_i(x) = s_(i-2)(x) - q_i(x)*s_(i-1)(x) 
 *     t_i(x) = t_(i-2)(x) - q_i(x)*t_(i-1)(x) 
 *   }
 *
 * Revision History
 *   Jun 1, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::euclid(int *a_x, int *b_x, int * & r_xi, int * &t_xi)
{
    int size = 2*t+1; // remember a(x) is order 2t;
    int *s_xi = new int[size];
    int *q_xi = new int[size];   
    int *r_xio = new int[size];
    int *s_xio = new int[size];
    int *t_xio = new int[size];
    
    int *temp = new int[size];
    
    // initialize r_xi-- yes i know the size of b_x is (size-1)-- r0 = b
    for(int i=0; i < size-1; i++)
        r_xi[i] = b_x[i];
    r_xi[size-1] = 0;

    for(int i=0; i< size; i++) // rn1 = a, tn1 = 0, s0 = 0
    {
        r_xio[i] = a_x[i];
        s_xio[i] = s_xi[i] = t_xio[i] = t_xi[i]= 0;
    }
    s_xio[0] = t_xi[0] = 1; // t0 = 1, sn1 = 1
   
    int i = 0;
       
    while(get_order(r_xi, size) >= t)
    {
        i++;
        copy_arr(temp, r_xi, size); // about to trash r_xi so save it
        doDiv(q_xi, r_xi, r_xio, temp, size); // update qxi and r_xi
        copy_arr(r_xio, temp, size);          // update r_xio to last r_xi
        
        copy_arr(temp, s_xi, size); // about to trash s_xi so save it
        doMul(s_xi, q_xi, temp,  size, size, size); 
        doSub(s_xi, s_xio, size);            // update s_xi
        copy_arr(s_xio, temp, size);         // update s_xio to last s_xi
        
        copy_arr(temp, t_xi, size); // about to trash t_xi so save it
        doMul(t_xi, q_xi, temp,  size, size, size); 
        doSub(t_xi, t_xio, size);            // update t_xi
        copy_arr(t_xio, temp, size);         // update t_xio to last t_xi
    }
        
    delete [] s_xi;                          // cleanup
    delete [] q_xi;
    delete [] r_xio;
    delete [] s_xio;
    delete [] t_xio;
    delete [] temp;
}


/*
 * chien(int *lambda, int * & roots)
 * Description:
 *   Find the roots of the error locator polynomial, lambda(x) using the Chien 
 *   search, which is an exhaustive search over all the elements in the Galois
 *   field.
 *
 *   - These roots are found by trial and error in which all possible values of 
 *     the roots (the GF(2^m) values alpha^i, 0 <= i <= n-1) are substituted 
 *     into the error locator polynomial. 
 *
 *   - If the lambda(x) expression evaluates to zero, upon substitution, then
 *     the value of x is a root, and identifies the error position at 
 *     (alpha^i)^-1 = alpha^(n-i).
 *     + This can be interpreted as an error at location (n-i)
 *
 *   - Note that if lambda(x) is written in the form:
 *      lambda(x) = X_1(x + X_1^-1) * X_2(x + X_2^-1)
 *     + then the roots are x=X_1^-1, X_2^-2, ...
 *     + the inverses of these roots are the error locations/positions
 *
 * Arguments:
 *   - lambda = lambda(x): error locator polynomial
 *   - roots  = array of index forms of possible roots, where alpha^i is only a
 *              root if root[i] == 1          
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   - all arrays are setup and have memory allocated
 *   - the roots array is of size n
 *   - lambda array is of size 2*t + 1
 *   - lambda(x)'s coefficients have been identified before calling this func.
 *
 * Operation:
 *   - loop through all possible roots (alpha^i, 0 <= i < n):
 *     + for each root plug it into lambda[j] and see if it evaluates 0
 *     + if it does, indicate that this index form of the field element 
 *       (alpha^i) is a root, by setting roots[i] = 1. Otherwise, set roots[i]=0
 *
 * Revision History
 *   Jun 01, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::chien(int *lambda, int * & roots)
{
    for(int i=0; i < n; i++)  // loop through all possible roots
    {   int lbd = 0;                // for each root plug it into lambda(x)
        for(int j=0; j< 2*t+1; j++) //
        {
            if(lambda[j])
                lbd ^= alpha_to[ (index_of[lambda[j]]+i*j)%n ];
        }
        if(lbd == 0) roots[i] = 1;   // record if this is a root
        else roots[i] = 0;
    }    
}


/*
 * forney(int *lambda, int *omega, int *err_loc)
 * Description:
 *   Use Forney's algorithm and calculate the error values at the error
 *   locations (as identified by a Chien search)
 *   Using these error values, fix the received codeword to get the decoded
 *   codeword, dc(x)
 *
 *   - Forney's algorithm states that the error values for a Reed-solomon code
 *     are computed by:
 *     e_ik = omega(X_k^-1)/derivative[lambda(X_k^-1)] where X_k^-1 is a root of
 *            lambda(x) and i is the index of the error location
 *
 * Assumptions:
 *   - all arrays are setup and have memory allocated
 *   - the err_loc array is of size n, and is set at indices where the
 *     rc_x has errors
 *   - lambda and omega arrays are of size 2*t + 1
 *  
 * Arguments:
 *   - lambda = lambda(x);  error locator polynomial
 *   - omega = omega(x); error magnitude/value/evaluator polynomial
 *   - err_loc = array that states if a given polynomial coeff.
 *               is an error location
 *
 * Return:
 *   None
 *
 * Operation:
 *   - Get derivative of lambda(x) = lambda(x)'
 *   - Using this determine decoded codeword, dc(x) with algorithm:
 *     + FOR all ccodeword indexed by i FROM 0 to (n-1)
 *         if (at error location) dc_i = rc_i + e
 *         else                   dc_i = rc_i
 *
 * Revision History
 *   Jun 01, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::forney(int *lambda, int *omega, int *err_loc)
{
    int size = 2*t+1;
    
    // lambdap(x) = derivative[lambda(x)] = lambda(x)'
    int * lambdap = new int[size];  
    for(int i=0; i < size-1; i++)
        lambdap[i] = ((i+1)%2)*lambda[i+1];
    lambdap[size-1] = 0;
    
    // setup decoded word and fix it
    for(int i=0; i < n; i++)
    {
        if(err_loc[i])
        {
            int num = 0;   // numerator
            int denum = 0; // denominator
            for(int j=0; j< size; j++)
            {
                if(omega[j])
                    num ^= alpha_to[ (index_of[omega[j]]+ ((n-i)%n)*j)%n ];
                if(lambdap[j])
                    denum ^= alpha_to[ (index_of[lambdap[j]]+ ((n-i)%n)*j)%n ];
            }
            dc_x[i] = rc_x[i]^alpha_to[(index_of[num] - index_of[denum] + n)%n];
        }
        else // err_loc[i] == 0  // rc(x) is correct at this index
            dc_x[i] = rc_x[i];
    }
    
    delete [] lambdap;
}


/*
 * decode()
 * Description:
 *   a Reed-Solomon decoder.
 *  
 *   ## Background:
 *   - The error locator polynomial:
 *     lambda(x) = 1 + lambda_1*x + ... + lambda_(v-1)*x^(v-1) + lambda_v*x^v
 *
 *   - For reach error there is a corresponding root X_j^-1 that makes lambda(x)
 *     equal to 0.
 *
 *   - The error magnitutude polynomial
 *     omega(x) = omega_(v-1)*x^(v-1) + ... + omega_1*x + omega_0;
 *     + also called error value or error evaluator polynomial
 * 
 *   - omega(x) = [s(x)*lambda(x)] mod x^2t, where s(x) is syndrome polynomial
 *     any terms of degree x^2t or higher in the product are ignored
 *
 *   - This representation of omega in terms of s(x) and lambda(x) is the key
 *     equation.
 *
 * Arguments:
 *   None
 *
 * Return:
 *   None
 *
 * Assumptions:
 *   None
 *
 * Operation:
 *   - Get the Syndrome polynomial.
 *   - Use the Euclidean Algorithm to get the Error Locator polynomial,
 *     lambda(x) and the needed omega(x)
 *   - Use Chien Search to get the roots of lambda(x) and then get the
 *     error locations from these roots (just the inverse).
 *   - Use Forney's Algorithm to determine the error values at error locations
 *     then add that to the received codeword, to get the decoded codeword 
 *
 * Revision History
 *   Jun 01, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void reedSolomon::decode()
{
    get_syndromes();  // first get syndromes
    if(!detect_error) // no error detected
    {                 // then dc(x) = rc(x)
        for(int i=0; i<n; i++)
            dc_x[i] = c_x[i];
        return;
    }
    
    // there is an error so get error locator polynomial
    int size = 2*t +1;
    int *r_xi = new int[size];
    int *t_xi = new int[size];
    int *a_x = new int[size];
    a_x[2*t] = 1;//a(x) = x^2t
    for(int i=0; i < 2*t; i++)
        a_x[i] = 0;
    
    // use euclidean algorithm
    euclid(a_x, s_x, r_xi, t_xi);  // results in weird scaling by alpha
    // omega(x) = r(x)_i and lambda(x) = t(x)_i
    // now scale lambda and omega by dividing by t(0)
    int pow = index_of[t_xi[0]];
    for(int i = 0; i < size; i++)
    {
        t_xi[i] =  (t_xi[i] ? alpha_to[(index_of[t_xi[i]]+n-pow)%n] : 0);
        r_xi[i] =  (r_xi[i] ? alpha_to[(index_of[r_xi[i]]+n-pow)%n] : 0);
    }
    
    int *lbd_roots = new int[n]; // lambda's roots are a^i when lbd_roots[i]!= 0
    chien(t_xi, lbd_roots);
    
    int *err_loc = new int[n];  // error_locators are inverse of roots of lambda
    for(int i = 0; i<n; i++)    // initialize error locators array
        err_loc[i] = 0;
    for(int i = 0; i<n; i++)    // then setup the values... note that if
        if (lbd_roots[i]) err_loc[(n-i)%n] = 1; // err_loc[i] != 0, then index i
                                                // is in error in rc_x [rc(x)]
    
    forney(t_xi, r_xi, err_loc); // run forney's algorithm -- also updates dc(x)

    delete [] r_xi;            // cleanup
    delete [] t_xi;
    delete [] a_x;
    delete [] lbd_roots;
    delete [] err_loc;
}

/*
 * compare()
 * Description:
 *   Check the results of this Reed-Solomon Encoder/Channel/Decoder
 *
 * Arguments:
 *    None
 *
 * Return:
 *   Success Boolean:
 *   - true indicates decoding worked. 
 *   - false indicates decoding failed
 * 
 * Assumptions:
 *   None
 *  
 * Operation:
 *   - start by assuming correct decoding
 *   - loop through all corresponding terms of both the transmitted codeword
 *     and the decoded codeword and compare them along the way.
 *     + if any pair arent equal flag an error in the decoder
 *   
 * Revision History
 *   Jun 01, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
bool reedSolomon::compare()
{
    bool correctly_decoded = true;    // start by assuming correct decoding
    for(int i = 0; i< n; i++ )        // if at any position, a mismatch occurs
        if(c_x[i] != dc_x[i]) correctly_decoded = false; // record an incorrect
    
    return correctly_decoded;
}
