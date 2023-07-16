/* ENC.C - Encoding procedures. */

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
#include <math.h>

#include "rand.h"
#include "alloc.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"
#include "rcode.h"
#include "enc.h"


/* The procedures in this module obtain the generator matrix to use for
   encoding from the global variables declared in rcode.h */


/* ENCODE A BLOCK USING A SPARSE REPRESENTATION OF THE GENERATOR MATRIX. */

void sparse_encode
( char *sblk,
  char *cblk
)
{
}


/* ENCODE A BLOCK USING DENSE REPRESENTATION OF GENERATOR MATRIX. */

void dense_encode
( char *sblk,
  char *cblk,
  mod2dense *G,
  mod2dense *u,
  mod2dense *v,
  int * cols
)
{
  int j;

  /* Copy source bits to the systematic part of the coded block. */

  for (j = M; j<N; j++) 
  { cblk[cols[j]] = sblk[j-M];
  }

  /* Multiply by Inv(A) X B to produce check bits. */

  for (j = M; j<N; j++)
  { mod2dense_set(u,j-M,0,sblk[j-M]); 
  }
  
  mod2dense_multiply(G,u,v);

  /* Copy check bits to the right places in the coded block. */

  for (j = 0; j<M; j++)
  { cblk[cols[j]] = mod2dense_get(v,j,0);
  }
}


/* ENCODE A BLOCK USING MIXED REPRESENTATION OF GENERATOR MATRIX. */

void mixed_encode
( char *sblk,
  char *cblk,
  mod2dense *u,
  mod2dense *v
)
{
}
