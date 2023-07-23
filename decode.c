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
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "open.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "check.h"
#include "dec.h"

#define LDPC_N 192
#define LDPC_M 56

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
  double tot_iter;		/* Double because can be huge for enum */
  double chngd, tot_changed;	/* Double because can be fraction if lratio==1*/
  double error_prob;

  int tot_valid;
  int valid;

  int i, j, k;

  /* Look at initial flag arguments. */

  table = 0;

  /* Look at arguments up to the decoding method specification. */

  error_prob = 0.02;

  /* Look at the specification of the decoding method, which starts at meth and
     continues to the end of the command line (marked by a zero pointer). */

  dec_method = Prprp;
  max_iter = 250;

  /* Read parity check file. */

  H = mod2sparse_read_H();

  M = mod2sparse_rows(H);
  N = mod2sparse_cols(H);

  /* Do the setup for the decoding method. */

  prprp_decode_setup();

  /* Read received blocks, decode, and write decoded blocks. */

  tot_iter = 0;
  tot_valid = 0;
  tot_changed = 0;

  for (block_no = 0; block_no < 1; block_no++)
  {
    /* Find likelihood ratio for each bit. */

    for (i = 0; i<N; i++)
    { lratio[i] = bsc_data[i]==1 ? (1-error_prob) / error_prob
                                  : error_prob / (1-error_prob);
    }

    /* Try to decode using the specified method. */

    iters = prprp_decode (H, lratio, dblk, pchk, bitpr);

    /* See if it worked, and how many bits were changed. */

    valid = check(H,dblk,pchk)==0;

    chngd = changed(lratio,dblk,N);

    tot_iter += iters;
    tot_valid += valid;
    tot_changed += chngd;

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
  "Decoded %d blocks, %d valid.  Average %.1f iterations, %.0f%% bit changes\n",
   block_no, tot_valid, (double)tot_iter/block_no, 
   100.0*(double)tot_changed/(N*block_no));

  exit(0);
}