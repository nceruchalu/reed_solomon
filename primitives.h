/* 
 * -----------------------------------------------------------------------------
 * -----                             primitives.h                          -----
 * -----                          REED-SOLOMON CODES                       -----
 * -----------------------------------------------------------------------------
 *
 * File Description:
 *   This is the list of finite field generator polynomials or primitive 
 *   polynomials, p(x) (mod 2) up to order NUM_PRIMITIVES
 *
 *   - A Primitive polynomial is a polynomial of degree m which is irreducible. 
 *     It forms part of the process of multiplying two field elements together. 
 *     For a Galois field of a particular size, there is sometimes a choice of
 *     suitable polynomials. Using a different primitive polynomial from that
 *     specified will produce incorrect results.
 *
 *   - As an example, for GF(16), the polynomial [p(x) = x^4 + x + 1] is
 *     irreducible and is thus a suitable primitive polynomial.
 *
 *
 * Operation:
 *   To get the primitive polynomial of order m from this table:
 *   - read the string at index (m-1), i.e. primitive[m-1]
 *   - translate the string apropriately by noting that only the degrees of
 *     the separate terms in primitive[m-1] are given. Examples:
 *     + "4 1 0" stands for: x^4 + x + 1
 *
 *
 * References:
 *   - http://downloads.bbc.co.uk/rd/pubs/whp/whp-pdf-files/WHP031.pdf
 *
 * Revision History
 *   May 30, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */

#define NUM_PRIMITIVES 100
static string primitive[NUM_PRIMITIVES] = {
  "1 0",
  "2 1 0",
  "3 1 0",
  "4 1 0",
  "5 2 0",
  "6 1 0",
  "7 1 0",
  "8 4 3 2 0",
  "9 4 0",
  "10 3 0",
  "11 2 0",
  "12 6 4 1 0",
  "13 4 3 1 0",
  "14 5 3 1 0",
  "15 1 0",
  "16 5 3 2 0",
  "17 3 0",
  "18 5 2 1 0",
  "19 5 2 1 0",
  "20 3 0",
  "21 2 0",
  "22 1 0",
  "23 5 0",
  "24 4 3 1 0",
  "25 3 0",
  "26 6 2 1 0",
  "27 5 2 1 0",
  "28 3 0",
  "29 2 0",
  "30 6 4 1 0",
  "31 3 0",
  "32 7 5 3 2 1 0",
  "33 6 4 1 0",
  "34 7 6 5 2 1 0",
  "35 2 0",
  "36 6 5 4 2 1 0",
  "37 5 4 3 2 1 0",
  "38 6 5 1 0",
  "39 4 0",
  "40 5 4 3 0",
  "41 3 0",
  "42 5 4 3 2 1 0",
  "43 6 4 3 0",
  "44 6 5 2 0",
  "45 4 3 1 0",
  "46 8 5 3 2 1 0",
  "47 5 0",
  "48 7 5 4 2 1 0",
  "49 6 5 4 0",
  "50 4 3 2 0",
  "51 6 3 1 0",
  "52 3 0",
  "53 6 2 1 0",
  "54 6 5 4 3 2 0",
  "55 6 2 1 0",
  "56 7 4 2 0",
  "57 5 3 2 0",
  "58 6 5 1 0",
  "59 6 5 4 3 1 0",
  "60 1 0",
  "61 5 2 1 0",
  "62 6 5 3 0",
  "63 1 0",
  "64 4 3 1 0",
  "65 4 3 1 0",
  "66 8 6 5 3 2 0",
  "67 5 2 1 0",
  "68 7 5 1 0",
  "69 6 5 2 0",
  "70 5 3 1 0",
  "71 5 3 1 0",
  "72 6 4 3 2 1 0",
  "73 4 3 2 0",
  "74 7 4 3 0",
  "75 6 3 1 0",
  "76 5 4 2 0",
  "77 6 5 2 0",
  "78 7 2 1 0",
  "79 4 3 2 0",
  "80 7 5 3 2 1 0",
  "81 4 0",
  "82 8 7 6 4 1 0",
  "83 7 4 2 0",
  "84 8 7 5 3 1 0",
  "85 8 2 1 0",
  "86 6 5 2 0",
  "87 7 5 1 0",
  "88 8 5 4 3 1 0",
  "89 6 5 3 0",
  "90 5 3 2 0",
  "91 7 6 5 3 2 0",
  "92 6 5 2 0",
  "93 2 0",
  "94 6 5 1 0",
  "95 6 5 4 2 1 0",
  "96 7 6 4 3 2 0",
  "97 6 0",
  "98 7 4 3 2 1 0",
  "99 7 5 4 0",
  "100 8 7 2 0"
};
