/* MOD2SPARSE.C - Procedures for handling sparse mod2 matrices. */

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


/* NOTE:  See mod2sparse.html for documentation on these procedures. */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mod2sparse.h"


/* ALLOCATE AN ENTRY WITHIN A MATRIX.  This local procedure is used to
   allocate a new entry, representing a non-zero element, within a given
   matrix.  Entries in this matrix that were previously allocated and
   then freed are re-used.  If there are no such entries, a new block
   of entries is allocated. */

static int num_blocks = 0;
static mod2block b_arr [90];


static mod2entry *alloc_entry
( mod2sparse *m
)
{ 
  mod2block *b;
  mod2entry *e;
  int k;

  if (m->next_free==0)
  { 
    b = b_arr + num_blocks;
    num_blocks++;

    b->next = m->blocks;
    m->blocks = b;

    for (k = 0; k<Mod2sparse_block; k++)
    { b->entry[k].left = m->next_free;
      m->next_free = &b->entry[k];
    }
  }

  e = m->next_free;
  m->next_free = e->left;

  e->pr = 0;
  e->lr = 0;

  return e;
}

/* READ A SPARSE MOD2 MATRIX STORED IN MACHINE-READABLE FORM FROM A FILE. */

mod2sparse *mod2sparse_read_H
()
{
  int n_rows, n_cols;
  int v, row, col;

  n_rows = 56;
  n_cols = 192;

  static mod2sparse m [1];
  static mod2entry m_rows [56];
  static mod2entry m_cols [192];

  mod2entry *e;
  int i, j;

  m->n_rows = n_rows;
  m->n_cols = n_cols;

  m->rows = m_rows;
  m->cols = m_cols;

  m->blocks = 0;
  m->next_free = 0;

  for (i = 0; i<n_rows; i++)
  { e = &m->rows[i];
    e->left = e->right = e->up = e->down = e;
    e->row = e->col = -1;
  }

  for (j = 0; j<n_cols; j++)
  { e = &m->cols[j];
    e->left = e->right = e->up = e->down = e;
    e->row = e->col = -1;
  }

  row = -1;

  int test_arr [] = {-1, 1, 2, 3, 4, 7, 10, 36, 45, 71, 83, 115, 137, 152, 156, 
                    177, 183, -2, 2, 8, 10, 13, 21, 23, 58, 64, 79, 101, 121, 130, 
                    148, 167, 169, 187, -3, 3, 4, 7, 19, 33, 37, 62, 65, 74, 89, 
                    151, 153, 160, 163, 171, 190, -4, 4, 5, 8, 46, 80, 82, 86, 
                    87, 106, 107, 113, 132, 167, 170, 184, 187, -5, 5, 6, 9, 19, 
                    28, 82, 84, 85, 88, 89, 90, 91, 92, 93, 95, 96, -6, 6, 8, 56, 
                    60, 64, 71, 72, 75, 81, 93, 122, 123, 127, 142, 166, 185, -7, 
                    7, 9, 11, 12, 22, 24, 26, 35, 75, 93, 128, 145, 150, 153, 160, 
                    175, -8, 8, 11, 21, 24, 34, 46, 51, 97, 113, 114, 118, 129, 133, 
                    145, 147, 161, -9, 9, 10, 13, 14, 20, 54, 59, 70, 76, 84, 110, 
                    111, 112, 124, 170, 186, -10, 10, 11, 12, 15, 23, 36, 49, 73, 129, 
                    142, 153, 155, 162, 172, 173, 180, -11, 11, 12, 15, 16, 22, 35, 49, 
                    50, 57, 82, 117, 159, 165, 173, 187, 188, -12, 12, 16, 18, 30, 33, 
                    35, 39, 50, 146, 148, 153, 154, 155, 158, 159, 160, -13, 13, 14, 17, 
                    23, 44, 53, 54, 55, 56, 57, 59, 60, 61, 62, 63, 64, -14, 14, 15, 16, 
                    20, 29, 33, 65, 67, 78, 122, 136, 144, 160, 162, 165, 175, -15, 15, 
                    22, 25, 27, 29, 31, 38, 94, 116, 117, 121, 122, 124, 126, 127, 128, 
                    -16, 16, 20, 40, 41, 46, 71, 77, 83, 104, 123, 126, 156, 158, 168, 
                    170, 177, -17, 17, 18, 21, 32, 55, 58, 68, 80, 88, 92, 101, 126, 
                    147, 158, 168, 184, -18, 18, 19, 22, 26, 56, 64, 77, 78, 98, 109, 
                    126, 131, 146, 158, 160, 169, -19, 19, 44, 57, 66, 68, 88, 89, 108, 
                    116, 119, 127, 135, 147, 155, 161, 183, -20, 20, 157, 178, 179, 180, 
                    181, 182, 184, 185, 186, 187, 188, 189, 190, 191, 192, -21, 21, 25, 
                    31, 51, 62, 69, 76, 85, 89, 116, 118, 124, 134, 154, 163, 179, -22, 
                    22, 25, 27, 28, 35, 36, 70, 71, 76, 81, 83, 116, 134, 152, 164, 168, 
                    -23, 23, 24, 26, 36, 40, 43, 66, 79, 104, 105, 137, 138, 140, 144, 
                    157, 183, -24, 24, 26, 43, 50, 59, 64, 74, 89, 114, 119, 124, 161, 
                    172, 173, 186, 188, -25, 25, 30, 34, 53, 83, 91, 95, 102, 103, 130, 
                    135, 136, 148, 150, 174, 189, -26, 26, 42, 47, 58, 65, 81, 109, 119, 
                    120, 123, 149, 150, 151, 166, 170, 172, -27, 27, 28, 29, 48, 52, 96, 
                    108, 119, 121, 129, 134, 149, 152, 179, 181, 191, -28, 28, 33, 47, 59, 
                    81, 84, 86, 102, 122, 156, 175, 182, 184, 187, 188, 189, -29, 29, 48, 
                    49, 51, 74, 76, 80, 87, 96, 107, 109, 112, 122, 129, 154, 180, -30, 30, 
                    45, 51, 62, 84, 97, 99, 103, 109, 138, 139, 153, 154, 157, 166, 191, -31, 
                    31, 39, 43, 62, 63, 69, 78, 81, 84, 99, 102, 106, 144, 152, 163, 192, -32, 
                    32, 37, 87, 130, 131, 132, 134, 135, 136, 138, 139, 140, 141, 142, 143, 144, 
                    -33, 33, 53, 59, 70, 86, 92, 132, 134, 157, 162, 176, 178, 182, 184, 185, 189, 
                    -34, 34, 38, 56, 72, 95, 97, 103, 111, 124, 133, 137, 139, 151, 165, 190, 191, 
                    -35, 35, 38, 48, 65, 70, 90, 92, 101, 111, 115, 143, 145, 159, 164, 168, 176, 
                    -36, 36, 37, 47, 52, 54, 60, 73, 98, 100, 140, 146, 149, 172, 173, 179, 186, 
                    -37, 37, 40, 41, 69, 90, 101, 110, 112, 117, 133, 135, 142, 146, 185, 190, 192, 
                    -38, 38, 45, 73, 77, 87, 100, 118, 120, 132, 167, 169, 172, 177, 178, 180, 182, 
                    -39, 39, 40, 43, 44, 66, 67, 88, 93, 128, 132, 157, 165, 171, 174, 183, 185, -40,
                    40, 41, 49, 52, 162, 163, 164, 165, 167, 168, 169, 171, 173, 174, 175, 176, -41, 
                    41, 42, 44, 47, 67, 69, 85, 90, 97, 107, 131, 141, 146, 161, 162, 164, -42, 42, 
                    55, 82, 87, 91, 100, 105, 110, 113, 121, 141, 152, 167, 170, 177, 191, -43, 43, 
                    98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 110, 111, 112, 125, -44, 44, 
                    58, 61, 63, 72, 73, 85, 91, 93, 100, 126, 131, 137, 169, 174, 177, -45, 45, 53, 
                    54, 57, 68, 80, 94, 99, 103, 104, 112, 114, 117, 136, 166, 186, -46, 46, 50, 52, 
                    58, 60, 61, 71, 72, 106, 123, 125, 127, 140, 149, 150, 161, -47, 47, 54, 55, 66, 
                    88, 96, 104, 105, 108, 128, 150, 155, 159, 163, 180, 182, -48, 48, 66, 67, 68, 69, 
                    70, 72, 73, 74, 75, 76, 77, 78, 79, 80, 86, -49, 49, 56, 57, 68, 75, 77, 79, 94, 95, 
                    118, 121, 130, 140, 145, 149, 189, -50, 50, 55, 82, 85, 90, 94, 105, 113, 114, 118, 
                    125, 131, 136, 143, 181, 188, -51, 51, 61, 78, 98, 106, 113, 114, 115, 125, 128, 133, 
                    139, 145, 155, 156, 176, -52, 1, 4, 5, 6, 14, 15, 17, 25, 28, 30, 61, 65, 67, 91, 110, 
                    120, -53, 2, 12, 17, 32, 34, 38, 52, 53, 92, 98, 102, 123, 143, 156, 158, 178, -54, 3, 
                    4, 13, 16, 30, 34, 63, 95, 120, 133, 135, 137, 138, 151, 181, 190, -55, 5, 6, 7, 9, 14, 
                    17, 27, 74, 75, 79, 107, 119, 130, 139, 148, 178, -56, 6, 7, 9, 10, 11, 13, 19, 21, 39, 
                    60, 120, 125, 127, 129, 141, 151, 0};

  int k = 0;
  
  for (;;)
  {
    v = test_arr [k++];

    if (v==0)
    { return m;
    }
    else if (v<0) 
    { row = -v-1;
    }
    else 
    { col = v-1;
      mod2sparse_insert(m,row,col);
    }
  }
}

/* INSERT AN ENTRY WITH GIVEN ROW AND COLUMN. */

mod2entry *mod2sparse_insert
( mod2sparse *m,
  int row,
  int col
)
{
  mod2entry *re, *ce, *ne;

  /* Find old entry and return it, or allocate new entry and insert into row. */

  re = mod2sparse_last_in_row(m,row);

  if (!mod2sparse_at_end(re) && mod2sparse_col(re)==col) 
  { return re;
  }

  if (mod2sparse_at_end(re) || mod2sparse_col(re)<col) 
  { re = re->right;
  }
  else
  {
    re = mod2sparse_first_in_row(m,row);

    for (;;)
    { 
      if (!mod2sparse_at_end(re) && mod2sparse_col(re)==col) 
      { return re;
      }

      if (mod2sparse_at_end(re) || mod2sparse_col(re)>col)
      { break;
      } 

      re = mod2sparse_next_in_row(re);
    }
  }
  ne = alloc_entry(m);

  ne->row = row;
  ne->col = col;

  ne->left = re->left;
  ne->right = re;
  ne->left->right = ne;
  ne->right->left = ne;

  /* Insert new entry into column.  If we find an existing entry here,
     the matrix must be garbled, since we didn't find it in the row. */

  ce = mod2sparse_last_in_col(m,col);

  if (!mod2sparse_at_end(ce) && mod2sparse_row(ce)==row) 
  { fprintf(stderr,"mod2sparse_insert: Garbled matrix\n");
    exit(1);
  }

  if (mod2sparse_at_end(ce) || mod2sparse_row(ce)<row) 
  { ce = ce->down;
  }
  else
  {
    ce = mod2sparse_first_in_col(m,col);

    for (;;)
    { 
      if (!mod2sparse_at_end(ce) && mod2sparse_row(ce)==row) 
      { fprintf(stderr,"mod2sparse_insert: Garbled matrix\n");
        exit(1);
      }

      if (mod2sparse_at_end(ce) || mod2sparse_row(ce)>row)
      { break;
      } 

      ce = mod2sparse_next_in_col(ce);
    }
  }
    
  ne->up = ce->up;
  ne->down = ce;
  ne->up->down = ne;
  ne->down->up = ne;

  /* Return the new entry. */

  return ne;
}

/* MULTIPLY VECTOR BY SPARSE MATRIX. */

void mod2sparse_mulvec
( mod2sparse *m,	/* The sparse matrix, with M rows and N columns */
  char *u,		/* The input vector, N long */
  char *v		/* Place to store the result, M long */
)
{
  mod2entry *e;
  int M, N;
  int i, j;

  M = mod2sparse_rows(m);
  N = mod2sparse_cols(m);

  for (i = 0; i<M; i++) v[i] = 0;

  for (j = 0; j<N; j++)
  { if (u[j])
    { for (e = mod2sparse_first_in_col(m,j);
           !mod2sparse_at_end(e);
           e = mod2sparse_next_in_col(e))
      { v[mod2sparse_row(e)] ^= 1;
      }
    }
  }
}