/* DEC.C - Decoding procedures. */

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


/* NOTE:  See decoding.html for general documentation on the decoding methods */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "alloc.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"
#include "rand.h"
#include "rcode.h"
#include "check.h"
#include "dec.h"


/* GLOBAL VARIABLES.  Declared in dec.h. */

decoding_method dec_method;	/* Decoding method to use */

int table;	/* Trace option, 2 for a table of decoding details */
int block_no;	/* Number of current block, from zero */

int max_iter;	/* Maximum number of iteratons of decoding to do */
char *gen_file;	/* Generator file for Enum_block and Enum_bit */

/* DECODE USING PROBABILITY PROPAGATION.  Tries to find the most probable 
   values for the bits of the codeword, given a parity check matrix (H), and
   likelihood ratios (lratio) for each bit.  If max_iter is positive, up to 
   that many iterations of probability propagation are done, stopping before 
   then if the tentative decoding is a valid codeword.  If max_iter is 
   negative, abs(max_iter) iterations are done, regardless of whether a 
   codeword was found earlier.

   Returns the number of iterations done (as an "unsigned" for consistency
   with enum_decode).  Regardless of whether or not a valid codeword was 
   reached, the bit vector from thresholding the bit-by-bit probabilities is 
   stored in dblk, and the resulting parity checks are stored in pchk (all 
   will be zero if the codeword is valid).  The final probabilities for each 
   bit being a 1 are stored in bprb.

   The setup procedure immediately below outputs headers for the detailed trace
   file, if required.
*/

void prprp_decode_setup (void)
{
  if (table==2)
  { printf(
     "  block  iter  changed  perrs    loglik   Eperrs   Eloglik  entropy\n");
  }
}

unsigned prprp_decode
( mod2sparse *H,	/* Parity check matrix */
  double *lratio,	/* Likelihood ratios for bits */
  char *dblk,		/* Place to store decoding */
  char *pchk,		/* Place to store parity checks */
  double *bprb		/* Place to store bit probabilities */
)
{ 
  int N, n, c;

  N = mod2sparse_cols(H);

  /* Initialize probability and likelihood ratios, and find initial guess. */

  initprp(H,lratio,dblk,bprb);

  /* Do up to abs(max_iter) iterations of probability propagation, stopping
     early if a codeword is found, unless max_iter is negative. */

  for (n = 0; ; n++)
  { 
    c = check(H,dblk,pchk);

    if (table==2)
    { printf("%7d %5d %8.1f %6d %+9.2f %8.1f %+9.2f  %7.1f\n",
       block_no, n, changed(lratio,dblk,N), c, loglikelihood(lratio,dblk,N), 
       expected_parity_errors(H,bprb), expected_loglikelihood(lratio,bprb,N),
       entropy(bprb,N));
    }
   
    if (n==max_iter || n==-max_iter || (max_iter>0 && c==0))
    { break; 
    }

    iterprp(H,lratio,dblk,bprb);
  }

  return n;
}


/* INITIALIZE PROBABILITY PROPAGATION.  Stores initial ratios, probabilities,
   and guess at decoding. */

void initprp
( mod2sparse *H,	/* Parity check matrix */
  double *lratio,	/* Likelihood ratios for bits */
  char *dblk,		/* Place to store decoding */
  double *bprb		/* Place to store bit probabilities, 0 if not wanted */
)
{ 
  mod2entry *e;
  int N;
  int j;

  N = mod2sparse_cols(H);

  for (j = 0; j<N; j++)
  { for (e = mod2sparse_first_in_col(H,j);
         !mod2sparse_at_end(e);
         e = mod2sparse_next_in_col(e))
    { e->pr = lratio[j];
      e->lr = 1;
    }
    if (bprb) bprb[j] = 1 - 1/(1+lratio[j]);
    dblk[j] = lratio[j]>=1;
  }
}


/* DO ONE ITERATION OF PROBABILITY PROPAGATION. */

void iterprp
( mod2sparse *H,	/* Parity check matrix */
  double *lratio,	/* Likelihood ratios for bits */
  char *dblk,		/* Place to store decoding */
  double *bprb		/* Place to store bit probabilities, 0 if not wanted */
)
{
  double pr, dl, t;
  mod2entry *e;
  int N, M;
  int i, j;

  M = mod2sparse_rows(H);
  N = mod2sparse_cols(H);

  /* Recompute likelihood ratios. */

  for (i = 0; i<M; i++)
  { dl = 1;
    for (e = mod2sparse_first_in_row(H,i);
         !mod2sparse_at_end(e);
         e = mod2sparse_next_in_row(e))
    { e->lr = dl;
      dl *= 2/(1+e->pr) - 1;
    }
    dl = 1;
    for (e = mod2sparse_last_in_row(H,i);
         !mod2sparse_at_end(e);
         e = mod2sparse_prev_in_row(e))
    { t = e->lr * dl;
      e->lr = (1-t)/(1+t);
      dl *= 2/(1+e->pr) - 1;
    }
  }

  /* Recompute probability ratios.  Also find the next guess based on the
     individually most likely values. */

  for (j = 0; j<N; j++)
  { pr = lratio[j];
    for (e = mod2sparse_first_in_col(H,j);
         !mod2sparse_at_end(e);
         e = mod2sparse_next_in_col(e))
    { e->pr = pr;
      pr *= e->lr;
    }
    if (isnan(pr))
    { pr = 1;
    }
    if (bprb) bprb[j] = 1 - 1/(1+pr);
    dblk[j] = pr>=1;
    pr = 1;
    for (e = mod2sparse_last_in_col(H,j);
         !mod2sparse_at_end(e);
         e = mod2sparse_prev_in_col(e))
    { e->pr *= pr;
      if (isnan(e->pr)) 
      { e->pr = 1;
      }
      pr *= e->lr;
    }
  }
}
