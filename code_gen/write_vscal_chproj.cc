
/*
 *
 *  Copyright UKQCD Collaboration, November 2000.
 *  Written by Peter Boyle.
 *  This software is provided for NON-COMMERCIAL use only,
 *  and may not be redistributed without permission.
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "processor.h"
extern struct processor* PROC;

#include "registers.h"


void qdp_vscal_chp( char *);
void qdp_vscal_chm( char *);


char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/

int main ( int argc, char *argv[])
{
  int arg;
  char *c;
  int proj_plus=0;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"n:P:pm")) != EOF){
    switch(arg){
    case 'n': 
      {
	if (strlen(optarg)<20) strcpy(name,optarg); 
      }
      break;
    case 'P': 
      {
	if (strlen(optarg)<20) strcpy(procname,optarg); 
      }
      break;
    case 'p': 
      {
	proj_plus=1; 
      } 
      break;
    default:
      {
	fprintf(stderr,"Usage: %s -nroutine_name -P proc [-p]\t",argv[0]); 
	fprintf(stderr, "\nenable -p option for chiral proj Plus\n");
	fprintf(stderr, "Otherwise it'll be chiral proj minus\n");
	exit (1);
      }
      break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  if( proj_plus == 1 ) { 
    fprintf(stderr, "Duing chiral proj plus\n");
    qdp_vscal_chp(name);
  }
  else { 
    fprintf(stderr, "Duing chiral proj minus\n");
    qdp_vscal_chm(name);
  }

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

/* qdp_vscal_chp(y, a, x, n_4vec)
 *
 * Operation: y = a * (1/2)(1 + g5) x
 *
 *  g5 is in chiral basis diag(1,1,-1, -1)
 * (1/2)(1+g5) = diag(1,1,0,0)
 *
 * A is a Scalar (Float *)
 * x is first element of a vector (Float *)
 * y is first element of a vector (Float *)
 * n_4vec is the number of 4 vectors
 *
 * for(i=0; i < n_4vec; i++) {
 *   for(spin=0; spin < 2: Nspin++) { 
 *     for(col=0; col < Ncol; col++) {
 *        y[i][spin][col] = a * x[i][spin][col];
 *     }
 *   }
 *   // Skip 2 spin components 
 * }
 */
void qdp_vscal_chp( char *name)
{
  int dum = defargcount(4);

  /*Integer register usage*/
  /* Pointers to vec1, vec2, vec3, a and  b*/
  /* and the counter */
  alreg(xptr,Iregs);
  alreg(yptr,Iregs);
  alreg(counter,Iregs);
  alreg(outptr,Iregs);
  alreg(Aptr, Iregs); /* Pointer to A */
  /*Floating register usage*/

  alreg(A, Fregs); /* A itself */
  reg_array_2d(x,Fregs,3,2);
  reg_array_2d(y,Fregs,3,2);
  alreg(zero, Fregs);
  alreg(tmp, Iregs);

  /* I will want to zero out the result (at outreptr and outimptr) */
  /* I need the word_size for that ... */
  int Isize = PROC->I_size;
  def_off(word_size,Byte, Isize);
  
  /* Various offsets. Start with 0 */
  def_dp_off(ZERO,0);

  /* A vector is of length 6 */
  def_dp_off(VEC_ATOM,6);

  /* This I am guessing is some pattern to describe the vectors
     and that it has to match the 2d register allocation above */
  offset_2d_dp(VEC_IMM,3,2);

  /* Declare memory streams for prefetching */
  struct stream *PreVecX;
  struct stream *PreVecY;

  int brchno,retno; /*Branch target handles*/
  int co,rei;



  /* Boilerplate */
  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /* Map the arguments -- in order */
  getarg(yptr); 
  getarg(Aptr);
  getarg(xptr);
  getarg(counter);

  /* get A into A from Aptr */
  queue_fload(A, ZERO, Aptr);


  /* Create a zero */
  /* Ugly hack to get floating zero */
  queue_iload_imm(tmp, ZERO);
  queue_istore(tmp, ZERO, yptr);
  if( PROC->I_size < PROC->FP_size ) { 
    // Write it to the first element of output -- will overwrite later
    queue_istore(tmp, word_size, yptr);
  }
  /*
   * Insert a label to prevent reordering
   */
  make_inst(DIRECTIVE,Target,get_target_label());
  

  /* Create the streams in and out */
  /* A 4spinor contians 4 atoms of length 3x2 */
  PreVecX= create_stream(VEC_ATOM,xptr, 4*counter,STREAM_IN ,LINEAR);
  PreVecY= create_stream(VEC_ATOM,yptr, 4*counter,STREAM_OUT,LINEAR);
  
  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */

  /* Load the zero value back from yptr */
  queue_fload(zero, ZERO, yptr);


  brchno = start_loop(counter);

  pragma(DCBT_SPACE,5);
  /* Load spin component 0 */
  for(co=0;co < 3; co++) { 
    for (rei=0;rei<2;rei++){
      queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
    }
  }

  iterate_stream(PreVecX); // Next col vec atom in X 4 spinor
  queue_prefetch(PreVecX); // Queue X prefetch
  

  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,yptr,ZERO);
  
  /* Multiply spin component 0 */
  for(co=0; co < 3; co++) { 
    for(rei=0;rei<2;rei++){
      queue_fmul(y[co][rei],A,x[co][rei]);
    }
  }


  for(co=0;co < 3; co++) {
    for(rei=0; rei < 2; rei++) {
      queue_fstore(y[co][rei],VEC_IMM[co][rei],outptr);
    }
  }

  iterate_stream(PreVecY); // Next col vec atom in Y 4 spinor

  /* Load spin component 1 */
  for(co=0;co < 3; co++) { 
    for (rei=0;rei<2;rei++){
      queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
    }
  }

  iterate_stream(PreVecX); // Next color component

  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,yptr,ZERO);

  /* Multiply spin component 0 */
  for(co=0; co < 3; co++) { 
    for(rei=0;rei<2;rei++){
      queue_fmul(y[co][rei],A,x[co][rei]);
    }
  }

  for(co=0;co < 3; co++) {
    for(rei=0; rei < 2; rei++) {
      queue_fstore(y[co][rei],VEC_IMM[co][rei],outptr);
    }
  }
  

  iterate_stream(PreVecY); // Next color component
  iterate_stream(PreVecX); // Next color component X
  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);

  /* Spin component 2 */
  /* Blank it with zeroes */
  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,yptr,ZERO);

  for(co=0; co < 3; co++) { 
    for(rei=0; rei < 2; rei++) { 
      queue_fstore(zero, VEC_IMM[co][rei], outptr);
    }
  }

  iterate_stream(PreVecY); // Next color component Y

  /* Spin component 3 */
  /* Blank it with zeroes */
  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,yptr,ZERO);

  for(co=0; co < 3; co++) { 
    for(rei=0; rei < 2; rei++) { 
      queue_fstore(zero, VEC_IMM[co][rei], outptr);
    }
  }

  iterate_stream(PreVecY);

  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}

/* qdp_vscal_chm(y, a, x, n_4vec)
 *
 * Operation: y = a * (1/2)(1 - g5) x
 * 
 *  g5 is in chiral basis diag(1,1,-1, -1)
 * (1/2)(1-g5) = diag(0,0,1,1)
 *
 * A is a Scalar (Float *)
 * x is first element of a vector (Float *)
 * y is first element of a vector (Float *)
 * n_4vec is the number of 4 vectors
 *
 * for(i=0; i < n_4vec; i++) {
 *   // Skip top 2 spin components
 *   for(spin=2; spin < Nspin: spin++) { 
 *     for(col=0; col < Ncol; col++) {
 *        y[i][spin][col] = a * x[i][spin][col];
 *     }
 *   }
 * }
 */
void qdp_vscal_chm( char *name)
{
  int dum = defargcount(4);

  /*Integer register usage*/
  /* Pointers to vec1, vec2, vec3, a and  b*/
  /* and the counter */
  alreg(xptr,Iregs);
  alreg(yptr,Iregs);
  alreg(counter,Iregs);
  alreg(outptr,Iregs);
  alreg(Aptr, Iregs); /* Pointer to A */
  /*Floating register usage*/

  alreg(A, Fregs); /* A itself */
  reg_array_2d(x,Fregs,3,2);
  reg_array_2d(y,Fregs,3,2);
  alreg(zero, Fregs);
  alreg(tmp, Iregs);

  /* I will want to zero out the result (at outreptr and outimptr) */
  /* I need the word_size for that ... */
  int Isize = PROC->I_size;
  def_off(word_size , Byte,Isize);
  
  /* Various offsets. Start with 0 */
  def_dp_off(ZERO,0);

  /* A vector is of length 6 */
  def_dp_off(VEC_ATOM,6);

  /* This I am guessing is some pattern to describe the vectors
     and that it has to match the 2d register allocation above */
  offset_2d_dp(VEC_IMM,3,2);

  /* Declare memory streams for prefetching */
  struct stream *PreVecX;
  struct stream *PreVecY;

  int brchno,retno; /*Branch target handles*/
  int co,rei;



  /* Boilerplate */
  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /* Map the arguments -- in order */
  getarg(yptr); 
  getarg(Aptr);
  getarg(xptr);
  getarg(counter);

  /* get A into A from Aptr */
  queue_fload(A, ZERO, Aptr);


  /* Create a zero */
  /* Ugly hack to get floating zero */
  queue_iload_imm(tmp, ZERO);
  queue_istore(tmp, ZERO, yptr);
  if( PROC->I_size < PROC->FP_size ) { 
    // Write it to the first element of output -- will overwrite later
    queue_istore(tmp, word_size, yptr);
  }
  /*
   * Insert a label to prevent reordering
   */
  make_inst(DIRECTIVE,Target,get_target_label());
  

  /* Create the streams in and out */
  /* A 4spinor contians 4 atoms of length 3x2 */
  PreVecX= create_stream(VEC_ATOM,xptr, 4*counter,STREAM_IN ,LINEAR);
  PreVecY= create_stream(VEC_ATOM,yptr, 4*counter,STREAM_OUT,LINEAR);
  
  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */
  /* Load the zero value back from yptr */
  queue_fload(zero, ZERO, yptr);

  brchno = start_loop(counter);

  pragma(DCBT_SPACE,5);

  /* Spin component 0 */
  /* Blank it with zeroes */
  /* This sets the outptr to the vecptr - with an add immediate */
  iterate_stream(PreVecX);
  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);

  // Skip 2 components in x and queue a prefetch
  queue_iadd_imm(outptr,yptr,ZERO);

  for(co=0; co < 3; co++) { 
    for(rei=0; rei < 2; rei++) { 
      queue_fstore(zero, VEC_IMM[co][rei], outptr);
    }
  }

  iterate_stream(PreVecY); // Next color component Y
  /* Spin component 1 */
  /* Blank it with zeroes */
  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,yptr,ZERO);

  for(co=0; co < 3; co++) { 
    for(rei=0; rei < 2; rei++) { 
      queue_fstore(zero, VEC_IMM[co][rei], outptr);
    }
  }


  iterate_stream(PreVecY);


  /* Load spin component 2 */
  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }

  co = 2;
  for (rei=0; rei<2; rei++) { 
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
  iterate_stream(PreVecX); // Next col vec atom in X 4 spinor
  queue_prefetch(PreVecX); // Queue X prefetch

  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,yptr,ZERO);
  
  /* Multiply spin component 2 */
  co = 0;
  for(rei=0;rei<2;rei++){
    queue_fmul(y[co][rei],A,x[co][rei]);
    queue_fstore(y[co][rei],VEC_IMM[co][rei],outptr);
  }

  co = 1;
  for(rei=0;rei<2;rei++){
    queue_fmul(y[co][rei],A,x[co][rei]);
    queue_fstore(y[co][rei],VEC_IMM[co][rei],outptr);
  }

  co = 2;
  for(rei=0;rei<2;rei++){
    queue_fmul(y[co][rei],A,x[co][rei]);
    queue_fstore(y[co][rei],VEC_IMM[co][rei],outptr);
  }



  iterate_stream(PreVecY); // Next col vec atom in Y 4 spinor



  /* Load spin component 3 */
  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }

  co = 2;
  for (rei=0; rei<2; rei++) { 
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
  iterate_stream(PreVecX); // Next color component
  queue_prefetch(PreVecX);

  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,yptr,ZERO);
  
  /* Multiply spin component 1 */
  co = 0;
  for(rei=0;rei<2;rei++){
    queue_fmul(y[co][rei],A,x[co][rei]);
    queue_fstore(y[co][rei],VEC_IMM[co][rei],outptr);
  }

  co = 1;
  for(rei=0;rei<2;rei++){
    queue_fmul(y[co][rei],A,x[co][rei]);
    queue_fstore(y[co][rei],VEC_IMM[co][rei],outptr);
  }

  co = 2;
  for(rei=0;rei<2;rei++){
    queue_fmul(y[co][rei],A,x[co][rei]);
    queue_fstore(y[co][rei],VEC_IMM[co][rei],outptr);
  }


  iterate_stream(PreVecY); // Next color component


  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








