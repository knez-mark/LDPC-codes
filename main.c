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

#include "encode.h"

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
  int i;

  LDPC_encode (sblk, cblk);

  printf ("Original message: ");
  for (i = 0; i < LDPC_K; i++) {
    if (i % 16 == 0) {
      printf ("\n");
    }
    printf ("%d, ", sblk[i]);

  }
  printf("\n");

  printf ("Encoded message: ");
  for (i = 0; i < LDPC_N; i++) {
    if (i % 16 == 0) {
      printf ("\n");
    }
    printf ("%d, ", cblk[i]);

  }
  printf("\n");
}