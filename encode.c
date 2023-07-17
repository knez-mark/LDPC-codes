/* ENCODE.C - Encode message blocks. */

/* Copyright (c) 1995-2012 by Radford M. Neal.
 *
 * Permission is granted for anyone to copy, use, modify, and distribute
 * these programs and accompanying documents for any purpose, provided
 * this copyright notice is retained and prominently displayed, and note
 * is made of any changes made to these programs.  These programs and
 * documents are distributed without any warranty, express or implied.
 * As the programs were written for research purposes only, they have not
 * been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own
 * risk.
 */

#include <stdio.h>

#include "mod2dense.h"
#include "enc.h"

#define LDPC_N 192
#define LDPC_K 122

static mod2dense G [1];		/* Dense or mixed representation of generator matrix,
			   if type=='d' or type=='m' */

static mod2word * G_col [LDPC_K*((LDPC_N - LDPC_K+mod2_wordsize-1) >> mod2_wordsize_shift)];
static mod2word * u_col [1], * v_col [1];

static int cols [] = {0x00, 0x01, 0x02, 0x03, 0x04, 
                      0x05, 0x06, 0x07, 0x08, 0x09, 
                      0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 
                      0x0f, 0x10, 0x11, 0x12, 0x13, 
                      0x14, 0x15, 0x16, 0x17, 0x18, 
                      0x19, 0x1a, 0x1b, 0x1c, 0x1d, 
                      0x1e, 0x1f, 0x20, 0x21, 0x22, 
                      0x23, 0x24, 0x25, 0x26, 0x27, 
                      0x28, 0x29, 0x2a, 0x2b, 0x2c, 
                      0x2d, 0x2e, 0x2f, 0x30, 0x31, 
                      0x32, 0x33, 0x34, 0x35, 0x36, 
                      0x37, 0x38, 0x39, 0x3a, 0x3b, 
                      0x3c, 0x3d, 0x3e, 0x3f, 0x40, 
                      0x41, 0x42, 0x43, 0x44, 0x45, 
                      0x46, 0x47, 0x48, 0x49, 0x4a, 
                      0x4b, 0x4c, 0x4d, 0x4e, 0x4f, 
                      0x50, 0x51, 0x52, 0x53, 0x54, 
                      0x55, 0x56, 0x57, 0x58, 0x59, 
                      0x5a, 0x5b, 0x5c, 0x5d, 0x5e, 
                      0x5f, 0x60, 0x61, 0x62, 0x63, 
                      0x64, 0x65, 0x66, 0x67, 0x68, 
                      0x69, 0x6a, 0x6b, 0x6c, 0x6d, 
                      0x6e, 0x6f, 0x70, 0x71, 0x72, 
                      0x73, 0x74, 0x75, 0x76, 0x77, 
                      0x78, 0x79, 0x7a, 0x7b, 0x7c, 
                      0x7d, 0x7e, 0x7f, 0x80, 0x81, 
                      0x82, 0x83, 0x84, 0x85, 0x86, 
                      0x87, 0x88, 0x89, 0x8a, 0x8b, 
                      0x8c, 0x8d, 0x8e, 0x8f, 0x90, 
                      0x91, 0x92, 0x93, 0x94, 0x95, 
                      0x96, 0x97, 0x98, 0x99, 0x9a, 
                      0x9b, 0x9c, 0x9d, 0x9e, 0x9f, 
                      0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 
                      0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 
                      0xaa, 0xab, 0xac, 0xad, 0xae, 
                      0xaf, 0xb0, 0xb1, 0xb2, 0xb3, 
                      0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 
                      0xb9, 0xba, 0xbb, 0xbc, 0xbd, 
                      0xbe, 0xbf};

// Generator matrix for LDPC
static mod2word G_bits [] = {0xa01800f8, 0x0002001c, 0x00000028, 0xa2520078, 0x1000a00c, 
                           0x00000028, 0x02c00800, 0x00000450, 0x00000000, 0x1694fd49, 
                           0xd081e056, 0x00000026, 0xc1ca8387, 0x1000a134, 0x00000028, 
                           0xc0869b38, 0x00400004, 0x00000028, 0x0014000d, 0x5001e004, 
                           0x00000020, 0x104a02c0, 0x1080a012, 0x00000000, 0x00000000, 
                           0x00000000, 0x00000008, 0x1006c379, 0x00080094, 0x00000020, 
                           0x12a69d3e, 0x8800a846, 0x00000001, 0x81d402f7, 0x0000e034, 
                           0x00000038, 0x8a867475, 0x00000614, 0x00000028, 0x03864243, 
                           0x80000074, 0x00000028, 0x8140586c, 0x02000024, 0x00000028, 
                           0x90c00078, 0x00800116, 0x00000028, 0x215856fa, 0x0802003c, 
                           0x00000000, 0xb0c05880, 0x00800016, 0x00000028, 0x824a16fa, 
                           0x0800a400, 0x00000029, 0x0252582e, 0xd201e040, 0x00000002, 
                           0xa0de4883, 0x0002001c, 0x00000028, 0x82544075, 0x5001e004, 
                           0x0000002a, 0x849eff89, 0x00014014, 0x0000002c, 0x1a38a5b3, 
                           0xd401e842, 0x00000002, 0x2015186f, 0x1000a018, 0x00000000, 
                           0x9235ff25, 0xd021e056, 0x0000002a, 0x865ee5c9, 0xd0014044, 
                           0x0000002f, 0x09c0660a, 0x01040325, 0x00000028, 0x10405808, 
                           0x05040201, 0x00000000, 0x8818bd00, 0x01044281, 0x00000010, 
                           0x18c03d76, 0x01040211, 0x00000028, 0x99006672, 0x01040211, 
                           0x00000000, 0xb8543470, 0x01060209, 0x00000000, 0x18de7c05, 
                           0x01044215, 0x00000038, 0x830c04b9, 0x8800a160, 0x00000029, 
                           0x215858f8, 0x00000138, 0x00000000, 0x186afe8d, 0x20000a02, 
                           0x00000000, 0x1000827a, 0x20080080, 0x00000000, 0x00061d7d, 
                           0x28000004, 0x00000000, 0x29d867f2, 0x20000338, 0x00000000, 
                           0x88062508, 0x2000a204, 0x00000028, 0x08983d06, 0x20004200, 
                           0x00000010, 0x90cac4f8, 0x0800a080, 0x00000029, 0x13c00000, 
                           0x00000050, 0x00000000, 0x84c0e704, 0x00000000, 0x00000000, 
                           0x89406616, 0x04000124, 0x00000028, 0x001818f8, 0x00000000, 
                           0x00000000, 0x00405800, 0x00080000, 0x00000000, 0xc08083f9, 
                           0x00101000, 0x00000008, 0xc246c13a, 0x00101044, 0x00000028, 
                           0x80921e8c, 0x0800e000, 0x00000039, 0x305e9fff, 0x00000088, 
                           0x00000000, 0x80801a7c, 0x00400000, 0x00000028, 0xc006c56c, 
                           0x02101014, 0x00000028, 0x307f9593, 0x0022000a, 0x00000000, 
                           0x100adb87, 0x1008a090, 0x00000000, 0x1100067a, 0x00800002, 
                           0x00000000, 0x4001856f, 0x00101000, 0x00000000, 0xb0461a13, 
                           0x02040001, 0x00000000, 0x8087c369, 0x00000014, 0x00000028, 
                           0x9235c326, 0xd021e056, 0x0000002a, 0x08003d4a, 0x20000200, 
                           0x00000000, 0x8080830a, 0x00000000, 0x00000000, 0x80805a01, 
                           0x0000a010, 0x00000008, 0x105e1af4, 0x00014015, 0x00000028, 
                           0x06c1bf10, 0x00000050, 0x00000000, 0xc094db7c, 0x00002014, 
                           0x00000028, 0x12a1e568, 0x00200052, 0x00000000, 0x10544309, 
                           0x5001e005, 0x00000028, 0x81d24674, 0x0801e000, 0x00000029, 
                           0x91861803, 0x00000001, 0x00000000, 0x8040d912, 0x02080081, 
                           0x00000000, 0x50469b84, 0x00000015, 0x00000028, 0x81c61a07, 
                           0x00040014, 0x00000028, 0x809e03f1, 0x5001c004, 0x00000028, 
                           0x26d2bb48, 0xd801e044, 0x00000006, 0x94e67c79, 0x00000016, 
                           0x00000028, 0x00000000, 0x00800000, 0x00000000, 0x94e7246d, 
                           0x00000016, 0x00000028, 0x90e08104, 0x00000016, 0x00000028, 
                           0x8646bd44, 0xd001a054, 0x0000002e, 0x809e4289, 0x00414004, 
                           0x00000028, 0x986cb7fe, 0x04000816, 0x00000028, 0x504cdbfc, 
                           0x0004b005, 0x00000028, 0x5040e701, 0x00801002, 0x00000000, 
                           0x641e7684, 0x00021008, 0x00000000, 0x90c61202, 0x00000014, 
                           0x00000028, 0x00800216, 0x02100000, 0x00000000, 0xc086c77a, 
                           0x08001014, 0x00000028, 0x0a40360a, 0x00400610, 0x00000000, 
                           0x90fe99f2, 0x00016816, 0x00000028, 0x0c808106, 0x00000000, 
                           0x00000000, 0x18808106, 0x00080080, 0x00000000, 0x8800bd00, 
                           0x00040081, 0x00000000, 0x22c61083, 0x00000400, 0x00000000, 
                           0xb092d970, 0x5009e084, 0x00000028, 0x809e4083, 0x04000014, 
                           0x00000028, 0x8c14990b, 0x5001e004, 0x0000002e, 0x80874216, 
                           0x00202814, 0x00000028, 0x001a18f5, 0x00000000, 0x00000000, 
                           0x02c05a04, 0x00001000, 0x00000000, 0xd0c4007d, 0x00101094, 
                           0x00000028, 0x82450a6e, 0x00000414, 0x00000028, 0x9221836a, 
                           0x00200042, 0x00000028, 0x0022c107, 0x00000010, 0x00000000, 
                           0x90e6c179, 0x00400004, 0x00000028, 0x008a1b90, 0x02002000, 
                           0x00000000, 0x18f9e5e8, 0x00204202, 0x00000000, 0x002ad9f2, 
                           0x1000a010, 0x00000000, 0x400ad9fd, 0x5001a010, 0x00000000, 
                           0x4000d9fd, 0x00501010, 0x00000000, 0x9180647e, 0x00040031, 
                           0x00000000, 0x8ade248d, 0x04004014, 0x00000028, 0x12b2c14b, 
                           0x80004852, 0x00000000, 0x06c0ffbe, 0xd0018040, 0x00000016, 
                           0x200c1899, 0x5001e018, 0x00000000, 0x00269d6b, 0x08000004, 
                           0x00000000, 0x90c0d942, 0x00040001, 0x00000000, 0xa0800862, 
                           0x0003401c, 0x00000028, 0x089864fe, 0x04400010, 0x00000000, 
                           0x001e1277, 0x00004004, 0x00000010, 0x22cc50fb, 0x0000a410, 
                           0x00000000};

static mod2word u_bits [(LDPC_K + mod2_wordsize - 1) >> mod2_wordsize_shift];
static mod2word v_bits [(LDPC_N - LDPC_K + mod2_wordsize - 1) >> mod2_wordsize_shift];

static char sblk [] = {0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 1, 
                       0, 0, 0, 0, 0, 0, 1, 0, 
                       0, 0, 0, 0, 0, 0, 1, 1, 
                       0, 0, 0, 0, 0, 1, 0, 0, 
                       0, 0, 0, 0, 0, 1, 0, 1, 
                       0, 0, 0, 0, 0, 1, 1, 0, 
                       0, 0, 0, 0, 0, 1, 1, 1, 
                       0, 0, 0, 0, 1, 0, 0, 0, 
                       0, 0, 0, 0, 1, 0, 0, 1, 
                       0, 0, 0, 0, 1, 0, 1, 0, 
                       0, 0, 0, 0, 1, 0, 1, 1, 
                       0, 0, 0, 0, 1, 1, 0, 0, 
                       0, 0, 0, 0, 1, 1, 0, 1, 
                       0, 0, 0, 0, 1, 1, 1, 0, 
                       0, 0};

static char cblk [LDPC_N];

/* MAIN PROGRAM. */

int main ()
{
  static mod2dense u [1], v [1];

  int i, n;

  /* Allocate space for generator matrix. */

  G->n_rows = LDPC_N - LDPC_K;
  G->n_cols = LDPC_K;
  G->n_words = (LDPC_N - LDPC_K + mod2_wordsize-1) >> mod2_wordsize_shift;

  G->col = G_col;
  G->bits = G_bits;

  for (i = 0; i<G->n_cols; i++)
  { 
    G->col[i] = G->bits + i*G->n_words;
  }

  /* Allocate needed space. */

  u->n_rows = LDPC_K;
  u->n_cols = 1;
  u->n_words = (LDPC_K + mod2_wordsize - 1) >> mod2_wordsize_shift;
  u->col = u_col;
  u->bits = u_bits;
  u->col[0] = u->bits;

  v->n_rows = LDPC_N - LDPC_K;
  v->n_cols = 1;
  v->n_words = (LDPC_N - LDPC_K + mod2_wordsize - 1) >> mod2_wordsize_shift;
  v->col = v_col;
  v->bits = v_bits;
  v->col[0] = v->bits;

  /* Encode successive blocks. */
  for (n = 0; ; n++)
  {
    /* Compute encoded block. */

    dense_encode (sblk, cblk, G, u, v, cols);

    printf ("encoded: ");
    for (i = 0; i < LDPC_N; i++) {
      if (i % 16 == 0) {
        printf ("\n");
      }
      printf ("%d, ", cblk[i]);

    }
    printf("\n");

    //In this case, break after encoding only once
    break;
  }

  fprintf(stderr,
    "Encoded %d blocks, source block size %d, encoded block size %d\n", n+1, LDPC_K, LDPC_N);

  return 0;
}