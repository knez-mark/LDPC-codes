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

#include "mod2dense.h"
#include "enc.h"
#include "encode.h"

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

  for (j = LDPC_N - LDPC_K; j<LDPC_N; j++) 
  { cblk[cols[j]] = sblk[j-(LDPC_N - LDPC_K)];
  }

  /* Multiply by Inv(A) X B to produce check bits. */

  for (j = LDPC_N - LDPC_K; j<LDPC_N; j++)
  { mod2dense_set(u,j-(LDPC_N - LDPC_K),0,sblk[j-(LDPC_N - LDPC_K)]); 
  }
  
  mod2dense_multiply(G,u,v);

  /* Copy check bits to the right places in the coded block. */

  for (j = 0; j<(LDPC_N - LDPC_K); j++)
  { cblk[cols[j]] = mod2dense_get(v,j,0);
  }
}