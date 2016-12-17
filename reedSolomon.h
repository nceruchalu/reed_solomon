/* 
 * -----------------------------------------------------------------------------
 * -----                            reedSolomon.h                          -----
 * -----                          REED-SOLOMON CODES                       -----
 * -----------------------------------------------------------------------------
 *
 * File Description:
 *   This is the header file for `reedSolomon` encoder/decoder class
 *
 * Assumptions:
 *   - Polynomials are written in vector form, or really an array where 
 *     zero-based index i represents the coefficient of x^i.
 *   - The factors read in the array slots are the Galois field elements in 
 *     decimal form, i.e. alpha^i values, where alpha is the primitive element
 *     of a primitive polynomial
 *
 * Revision History
 *   May 31, 2011    Nnoduka Eruchalu    Initial Revision- encoder works
 *   Jun 01, 2011    Nnoduka Eruchalu    added decoder
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */

#ifndef REED_SOLOMON_H_
#define REED_SOLOMON_H_

class reedSolomon
{
private:
    ///////////////////////////////////////////////////////////////////////////
    // Instance Variables
    ///////////////////////////////////////////////////////////////////////////
    int m;             // messages and codewords are over finite GF(2^m) field
    int t;             // number of error symbols that can be corrected
                       // i.e. this is a t-error correcting code
    int n;             // codeword length (in symbols): should be 2^m -1;
    int k;             // number of message symbols: n-k = 2t ==> k = n-2t
    
    bool detect_error; // any error detected in received vector?
    
    /* For a given Galois field, alpha is the primitive element of the chosen
     * primitive polynomial. The entire GF is represented by the arrays
     * `alpha_to` and `index_of`. 
     *
     * For a given m and GF(2^m):
     * - alpha_to: each element of the array is the decimal form of 
     *             a GF element. Indices into this array are the corresponding 
     *             index forms, such that 
                   alpha^i == alpha_to[i].
     * - index_of: is the inverse of `alpha_to` where each element of the array 
     *             is the index form of of a GF element. Indices to this array
     *             are the corresponding decimal forms, such that:
     *             i == index_of[j] and alpha^(index_of[j]) = j
     *   
     * - Example:  For m=4: alpha_to[3] == 8  and index_of[8] == 3
     *             For m=4: alpha_to[9] == 10 and index_of[10] == 9
     */ 
    int *alpha_to;    // decimal form of polynomial corresponding to alpha^i
    int *index_of;    // inverse of alpha_to; each entry cont

    int *m_x;         // m(x): message polynomial
    int *g_x;         // g(x): generator polynomial
    int *p_x;         // p(x) primitive polynomial
    int *c_x;         // c(x): encoded codeword
    int *rc_x;        // received code_word: c(x) after going through channel
    int *s_x;         // s(x): syndromes polynomial
    int *dc_x;        // dc(x): decoded codeword
    
    int size_g;       // size of g(x);
    int size_p;       // size of p(x);
    int size_s;       // size of s(x)
    // size_m = k; size_c = n;  // for m(x), c(x) respectively
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Initializers: Initialize instance variables
    ///////////////////////////////////////////////////////////////////////////
    void gen_gf();          // update alpha_to and index_of
    void gen_prim_poly();   // get primitive polynomial, p(x), for given m
    void gen_g_poly();      // create generator polynomial, g(x)
    void update_params();   // initialize instance variables
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Galios field arithmetic
    ///////////////////////////////////////////////////////////////////////////
    //c(x) = m1(x) + m2(x)
    void doAdd(int * & cc, int *m1, int *m2, int size_n); 
    
    // m1(x) -= m2(x)   or   m1(x) += m2(x)
    void doSub(int * & m1, int * & m2, int size_n);      
    
    // res1(x) = g(x) *(alpha^fact)* x^w
    void mul(int * & res1, int *g1, int fact, int w, int size_res1, 
             int size_g1);                       
    
    // res(x) = mm(x) * x^w
    void toPower(int * & res, int * mm, int size_mm, int w);
    
    // res(x) = a(x)*b(x)
    void doMul(int * & res_x, int * a_x, int * b_x,  int size_res,
               int size_a, int size_b);
    
    // quot(x) = a(x)/b(x)   and rmd(x) = a(x) % b(x)
    void doDiv(int * & quot_x, int * & rmd_x, int * a_x, int * b_x,
               int size);
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Miscellaneous routines
    ///////////////////////////////////////////////////////////////////////////
    void print_poly(int *p_x, int size_p); // print polynomial for humans
    int get_order(int *p_x, int size_p);   // get order of polynomial p(x)
    
    void copy_arr(int * & dst,int * src, int size); // dst[] = src[]
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Reed-Solomon encoding/decoding methods
    ///////////////////////////////////////////////////////////////////////////
    // populate s(x)
    void get_syndromes(); 
    
    // get the coefficients of the error locator polynomial
    void euclid(int *a_x, int *b_x, int * & r_xi, int * & t_xi);
    
    // find the roots of the error locator polynomial
    void chien(int *lambda, int * & roots);
    
    // calculate error values at the error locations (as identified by chien())
    void forney(int *lambda, int *omega, int *err_loc);
    
public:
    ///////////////////////////////////////////////////////////////////////////
    // Initializers: Constructors & Desctructor
    ///////////////////////////////////////////////////////////////////////////
    reedSolomon();               // default constructor
    reedSolomon(int mm, int tt); // constructor that takes object's m,t values
    ~reedSolomon();              // destructor
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Miscellaneous routines
    ///////////////////////////////////////////////////////////////////////////
    void print_params();  // print the many parameters and polynomials
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Reed-Solomon encoding/decoding methods
    ///////////////////////////////////////////////////////////////////////////
    void gen_rand_msg();         // create random message m(x)
    void encode();               // encode created message, m(x), to get c(x)
    void sim_channel();          // generate rc(x) with at most t errors
    void sim_channel(double Ps); // generate rc(x) with Ps symbol error prob.
    void decode();               // decode received vector rc(x);
    bool compare();              // compare c(x) and dc(x)


	///////////////////////////////////////////////////////////////////////////
	// Reed-Solomon set/get method, length of msg must be n
	// By yuanyuanxiang, yuan_yuanxiang@163.com
	///////////////////////////////////////////////////////////////////////////
	void get_original_msg(int *msg);
	void get_received_msg(int *msg);
	void set_msgfor_encode(int *msg);	// set message m(x)
	void get_encoded_msg(int *msg);		// get message c(x)
	void set_msgfor_decode(int *msg);	// set message rc(x)
	void get_decoded_msg(int *msg);		// get message dc(x)
	void rs_encode(int *msg);			// rs_encode by yuanyuanxiang
	void rs_decode(int *msg);			// rs_decode by yuanyuanxiang
};

#endif                                                     /* REED_SOLOMON_H_ */
