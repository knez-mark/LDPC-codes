/* DECODE.C - Decode blocks of received data. */

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

#include "mod2sparse.h"
#include "check.h"
#include "dec.h"

#define LDPC_N 192
#define LDPC_M 56
#define ERROR_PROB 0.02

/* MAIN PROGRAM. */

int main
( int argc,
  char **argv
)
{
  int M, N;
  mod2sparse *H;

  static char dblk [LDPC_N], pchk [LDPC_M];
  static double lratio [LDPC_N];
  static double bitpr [LDPC_N];

  //Six bits altered
  int bsc_data [] = {1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 
                     0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 
                     1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 
                     0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 
                     0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 
                     0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 
                     0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 
                     0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 
                     0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 
                     0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0};

  unsigned iters;		/* Unsigned because can be huge for enum */
  double chngd;	/* Double because can be fraction if lratio==1*/

  int valid, block_no;

  int i, j;

  /* Read parity check file. */

  H = mod2sparse_read_H();

  M = mod2sparse_rows(H);
  N = mod2sparse_cols(H);

  /* Read received blocks, decode, and write decoded blocks. */

  for (block_no = 0; block_no < 1; block_no++)
  {
    /* Find likelihood ratio for each bit. */

    for (i = 0; i<N; i++)
    { lratio[i] = bsc_data[i]==1 ? (1-ERROR_PROB) / ERROR_PROB
                                  : ERROR_PROB / (1-ERROR_PROB);
    }

    /* Try to decode using the specified method. */

    iters = prprp_decode (H, lratio, dblk, pchk, bitpr);

    /* See if it worked, and how many bits were changed. */

    valid = check(H,dblk,pchk)==0;

    chngd = changed(lratio,dblk,N);

    /* Write decoded block. */

    printf("Decoded data");
    for (j = 0; j < N; j++) {
      if (j%16 == 0) {
        printf ("\n");
      }
      printf ("%d, ", dblk [j]);
    }
    printf ("\n");
  }

  /* Finish up. */

  fprintf(stderr,
  "Decoded %d blocks, %d valid.  Average %d iterations, %f bit changes\n",
   block_no, valid, iters, chngd);

  return 0;
}