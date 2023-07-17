/* MOD2DENSE.C - Procedures for handling dense mod2 matrices. */

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


/* NOTE:  See mod2dense.html for documentation on these procedures. */


#include <stdio.h>

#include "mod2dense.h"

/* CLEAR A DENSE MOD2 MATRIX. */

void mod2dense_clear
( mod2dense *r
)
{
  int k, j;

  for (j = 0; j<mod2dense_cols(r); j++)
  { for (k = 0; k<r->n_words; k++)
    { r->col[j][k] = 0;
    }
  }

}

/* GET AN ELEMENT FROM A DENSE MOD2 MATRIX. */

int mod2dense_get  
( mod2dense *m, 	/* Matrix to get element from */
  int row,		/* Row of element (starting with zero) */
  int col		/* Column of element (starting with zero) */
)
{
  if (row<0 || row>=mod2dense_rows(m) || col<0 || col>=mod2dense_cols(m))
  { fprintf(stderr,"mod2dense_get: row or column index out of bounds\n");
  }

  return mod2_getbit (m->col[col][row>>mod2_wordsize_shift], 
                      row&mod2_wordsize_mask);
}


/* SET AN ELEMENT IN A DENSE MOD2 MATRIX. */

void mod2dense_set 
( mod2dense *m, 	/* Matrix to modify element of */
  int row,		/* Row of element (starting with zero) */
  int col,		/* Column of element (starting with zero) */
  int value		/* New value of element (0 or 1) */
)
{ 
  mod2word *w;

  if (row<0 || row>=mod2dense_rows(m) || col<0 || col>=mod2dense_cols(m))
  { fprintf(stderr,"mod2dense_set: row or column index out of bounds\n");
  }

  w = &m->col[col][row>>mod2_wordsize_shift];

  *w = value ? mod2_setbit1(*w,row&mod2_wordsize_mask) 
             : mod2_setbit0(*w,row&mod2_wordsize_mask);
}

/* MULTIPLY TWO DENSE MOD2 MATRICES. 

   The algorithm used runs faster if the second matrix (right operand of the
   multiply) is sparse, but it is also appropriate for dense matrices.  This
   procedure could be speeded up a bit by replacing the call of mod2dense_get
   with in-line code that avoids division, but this doesn't seem worthwhile
   at the moment. 
*/

void mod2dense_multiply 
( mod2dense *m1, 	/* Left operand of multiply */
  mod2dense *m2,	/* Right operand of multiply */
  mod2dense *r		/* Place to store result of multiply */
)
{
  int i, j, k;

  if (mod2dense_cols(m1)!=mod2dense_rows(m2) 
   || mod2dense_rows(m1)!=mod2dense_rows(r) 
   || mod2dense_cols(m2)!=mod2dense_cols(r))
  { fprintf(stderr,
     "mod2dense_multiply: Matrices have incompatible dimensions\n");
  }

  if (r==m1 || r==m2)
  { fprintf(stderr,
      "mod2dense_multiply: Result matrix is the same as one of the operands\n");
  }

  mod2dense_clear(r);

  for (j = 0; j<mod2dense_cols(r); j++)
  { for (i = 0; i<mod2dense_rows(m2); i++)
    { if (mod2dense_get(m2,i,j))
      { for (k = 0; k<r->n_words; k++)
        { r->col[j][k] ^= m1->col[i][k];
        }
      }
    }
  }
}