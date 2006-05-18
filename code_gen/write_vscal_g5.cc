
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


void qdp_vscal_g5( char *);

char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"n:P:")) != EOF){
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
    default:
      {
	fprintf(stderr,"Usage: %s -nroutine_name -P proc \t\n",argv[0]); 
	exit (1);
      }
      break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  qdp_vscal_g5(name);


  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

/* qdp_vscal_g5(y, a, x, n_4vec)
 *
 * Operation: y = a g5 x
 *
 *  g5 is in chiral basis diag(1,1,-1, -1)
 *
 * A is a Scalar (Float *)
 * x is first element of a vector (Float *)
 * y is first element of a vector (Float *)
 * n_4vec is the number of 4 vectors
 *
 * for(i=0; i < n_4vec; i++) {
 *   for(spin=0; spin < 2: spin++) { 
 *     for(col=0; col < Ncol; col++) {
 *        y[i][spin][col] = a * x[i][spin][col];
 *     }
 *   }
 *   for(spin=2; spin < 4: spin++) { 
 *     for(col=0; col < Ncol; col++) {
 *        y[i][spin][col] = -a * x[i][spin][col];
 *     }
 *   }
 * }
 */
void qdp_vscal_g5( char *name)
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
  reg_array_2d(my,Fregs, 3, 2); /* For negated y-s */
  
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
  queue_prefetch(PreVecX);
  brchno = start_loop(counter);

  pragma(DCBT_SPACE,5);

  /* ******************************* SPIN 0 *************************/
  /* Load spin component */
  for(co=0;co < 3; co++) { 
    for (rei=0;rei<2;rei++){
      queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
    }
  }

  iterate_stream(PreVecX); // Next col vec atom in X 4 spinor
  queue_prefetch(PreVecX); // Queue X prefetch
  

  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,yptr,ZERO);
  
  /* Multiply spin component */
  for(co=0; co < 3; co++) { 
    for(rei=0;rei<2;rei++){
      queue_fmul(y[co][rei],A,x[co][rei]);
    }
  }

  /* Store spin component */
  for(co=0;co < 3; co++) {
    for(rei=0; rei < 2; rei++) {
      queue_fstore(y[co][rei],VEC_IMM[co][rei],outptr);
    }
  }

  iterate_stream(PreVecY); // Next col vec atom in Y 4 spinor
  /************************ END OF SPIN COMPONENT ****************/

  /************************ SPIN COMPONENT 1 *********************/

  /* Load spin component */
  for(co=0;co < 3; co++) { 
    for (rei=0;rei<2;rei++){
      queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
    }
  }

  iterate_stream(PreVecX); // Next color component
  queue_prefetch(PreVecX);

  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,yptr,ZERO);

  /* Multiply spin component */
  for(co=0; co < 3; co++) { 
    for(rei=0;rei<2;rei++){
      queue_fmul(y[co][rei],A,x[co][rei]);
    }
  }

  /* Store spin component  */
  for(co=0;co < 3; co++) {
    for(rei=0; rei < 2; rei++) {
      queue_fstore(y[co][rei],VEC_IMM[co][rei],outptr);
    }
  }
  
  iterate_stream(PreVecY); // Next color component
  /************************ END OF SPIN COMPONENT ***************/

  /************************ SPIN COMPONENT 2 *********************/

  /* Load spin component */
  for(co=0;co < 3; co++) { 
    for (rei=0;rei<2;rei++){
      queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
    }
  }

  iterate_stream(PreVecX); // Next color component
  queue_prefetch(PreVecX);

  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,yptr,ZERO);

  /* Multiply spin component */
  for(co=0; co < 3; co++) { 
    for(rei=0;rei<2;rei++){
      queue_fmul(y[co][rei],A,x[co][rei]);
    }
  }

  /* Store spin component  */
  for(co=0;co < 3; co++) {
    for(rei=0; rei < 2; rei++) {
      queue_fneg(my[co][rei], y[co][rei]);
      queue_fstore(my[co][rei],VEC_IMM[co][rei],outptr);
    }
  }
  
  iterate_stream(PreVecY); // Next color component
  /************************ END OF SPIN COMPONENT ***************/

  /************************ SPIN COMPONENT 3 *********************/

  /* Load spin component */
  for(co=0;co < 3; co++) { 
    for (rei=0;rei<2;rei++){
      queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
    }
  }

  iterate_stream(PreVecX); // Next color component
  queue_prefetch(PreVecX);

  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,yptr,ZERO);

  /* Multiply spin component */
  for(co=0; co < 3; co++) { 
    for(rei=0;rei<2;rei++){
      queue_fmul(y[co][rei],A,x[co][rei]);
    }
  }

  /* Store spin component  */
  for(co=0;co < 3; co++) {
    for(rei=0; rei < 2; rei++) {
      queue_fneg(my[co][rei], y[co][rei]);
      queue_fstore(my[co][rei],VEC_IMM[co][rei],outptr);
    }
  }
  
  iterate_stream(PreVecY); // Next color component
  /************************ END OF SPIN COMPONENT ***************/

  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}









