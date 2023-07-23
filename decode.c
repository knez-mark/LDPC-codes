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

#include "rand.h"
#include "alloc.h"
#include "blockio.h"
#include "open.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "channel.h"
#include "check.h"
#include "dec.h"

#define LDPC_N 192
#define LDPC_M 56

void usage(void);

/* MAIN PROGRAM. */

int main
( int argc,
  char **argv
)
{
  int M, N;
  mod2sparse *H;

  char *pchk_file, *rfile, *dfile, *pfile;
  char **meth;
  FILE *rf, *df, *pf;

  char *dblk, *pchk;
  double *lratio;
  double *bitpr;

  //Four bits altered
  int bsc_data [] = {1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 
                     0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 
                     1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 
                     0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 
                     0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 
                     0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 
                     0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 
                     0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 
                     0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 
                     0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0};

  unsigned iters;		/* Unsigned because can be huge for enum */
  double tot_iter;		/* Double because can be huge for enum */
  double chngd, tot_changed;	/* Double because can be fraction if lratio==1*/

  int tot_valid;
  int valid;

  int i, j, k;

  /* Look at initial flag arguments. */

  table = 0;
  blockio_flush = 0;

  while (argc>1)
  {
    if (strcmp(argv[1],"-t")==0)
    { if (table!=0) usage();
      table = 1;
    }
    else if (strcmp(argv[1],"-T")==0)
    { if (table!=0) usage();
      table = 2;
    }
    else if (strcmp(argv[1],"-f")==0)
    { if (blockio_flush!=0) usage();
      blockio_flush = 1;
    }
    else 
    { break;
    }

    argc -= 1;
    argv += 1;
  }

  /* Look at arguments up to the decoding method specification. */

  if (!(pchk_file = argv[1])
   || !(rfile = argv[2])
   || !(dfile = argv[3]))
  { usage();
  }

  if (argv[4]==0 || argv[5]==0) usage();

  k = channel_parse(argv+4,argc-4);
  if (k<=0)
  { pfile = argv[4];
    k = channel_parse(argv+5,argc-5);
    if (k<=0) usage();
    meth = argv+5+k;
  }
  else
  { pfile = 0;
    meth = argv+4+k;
  }

  /* Look at the specification of the decoding method, which starts at meth and
     continues to the end of the command line (marked by a zero pointer). */

  dec_method = Prprp;
  max_iter = 250;

  /* Check that we aren't overusing standard input or output. */

  if ((strcmp(pchk_file,"-")==0) 
    + (strcmp(rfile,"-")==0) > 1)
  { fprintf(stderr,"Can't read more than one stream from standard input\n");
    exit(1);
  }
  if ((table>0) 
    + (strcmp(dfile,"-")==0) 
    + (pfile!=0 && strcmp(pfile,"-")==0) > 1)
  { fprintf(stderr,"Can't send more than one stream to standard output\n");
    exit(1);
  }

  /* Read parity check file. */

  H = mod2sparse_read_H();

  M = mod2sparse_rows(H);
  N = mod2sparse_cols(H);

  if (N<=M)
  { fprintf(stderr,
     "Number of bits (%d) should be greater than number of checks (%d)\n",N,M);
    exit(1);
  }

  /* Create file for decoded data. */

  df = open_file_std(dfile,"w");
  if (df==NULL)
  { fprintf(stderr,"Can't create file for decoded data: %s\n",dfile);
    exit(1);
  }

  /* Allocate other space. */

  dblk   = chk_alloc (N, sizeof *dblk);
  lratio = chk_alloc (N, sizeof *lratio);
  pchk   = chk_alloc (M, sizeof *pchk);
  bitpr  = chk_alloc (N, sizeof *bitpr);

  /* Print header for summary table. */

  if (table==1)
  { printf("  block iterations valid  changed\n");
  }

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

    /* Print summary table entry. */

    if (table==1)
    { printf ("%7d %10f    %d  %8.1f\n",
        block_no, (double)iters, valid, (double)chngd);
        /* iters is printed as a double to avoid problems if it's >= 2^31 */
      fflush(stdout);
    }

    /* Write decoded block. */

    blockio_write(df,dblk,N);

    /* Write bit probabilities, if asked to. */

    if (pfile)
    { for (j = 0; j<N; j++)
      { fprintf(pf," %.5f",bitpr[j]);
      }
      fprintf(pf,"\n");
    }

    /* Check for errors when writing. */

    if (ferror(df) || pfile && ferror(pf))
    { break;
    }
  }

  /* Finish up. */

done: 
  fprintf(stderr,
  "Decoded %d blocks, %d valid.  Average %.1f iterations, %.0f%% bit changes\n",
   block_no, tot_valid, (double)tot_iter/block_no, 
   100.0*(double)tot_changed/(N*block_no));

  if (ferror(df) || fclose(df)!=0)
  { fprintf(stderr,"Error writing decoded blocks to %s\n",dfile);
    exit(1);
  }

  if (pfile)
  { if (ferror(pf) || fclose(pf)!=0)
    { fprintf(stderr,"Error writing bit probabilities to %s\n",dfile);
      exit(1);
    }
  }

  exit(0);
}


/* PRINT USAGE MESSAGE AND EXIT. */


void usage(void)
{ fprintf(stderr,"Usage:\n");
  fprintf(stderr,
"  decode [ -f ] [ -t | -T ] pchk-file received-file decoded-file [ bp-file ] channel method\n");
  channel_usage();
  fprintf(stderr,
"Method:  enum-block gen-file | enum-bit gen-file | prprp [-]max-iterations\n");
  exit(1);
}
