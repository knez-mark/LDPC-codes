/* CHECK.C - Compute parity checks and other stats on decodings. */

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mod2sparse.h"
#include "check.h"



/* COMPUTE PARITY CHECKS.  Returns the number of parity checks violated by
   dblk.  The results of all the parity checks are stored in pchk. */

int check
( mod2sparse *H,	/* Parity check matrix */
  char *dblk,		/* Guess for codeword */
  char *pchk		/* Place to store parity checks */
)
{
  int M, i, c;

  M = mod2sparse_rows(H);

  mod2sparse_mulvec (H, dblk, pchk);

  c = 0;
  for (i = 0; i<M; i++) 
  { c += pchk[i];
  }

  return c;
}


/* COUNT HOW MANY BITS HAVED CHANGED FROM BIT INDICATED BY LIKELIHOOD.  The
   simple decoding based on likelihood ratio is compared to the given decoding.
   A bit for which the likelihood ratio is exactly one counts as half a 
   change, which explains why the result is a double rather than an int.
 */

double changed
( double *lratio,	/* Likelihood ratios for bits */
  char *dblk,		/* Candidate decoding */
  int N			/* Number of bits */
)
{ 
  double changed;
  int j;

  changed = 0;
  for (j = 0; j<N; j++)
  { changed += lratio[j]==1 ? 0.5 : dblk[j] != (lratio[j]>1); 
  }

  return changed;
}