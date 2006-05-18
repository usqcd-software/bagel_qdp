
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


void qdp_lsum2( char *);

char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"n:P:m")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<20) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<20) strcpy(procname,optarg); break;
    default: fprintf(stderr,"Usage: %s -nroutine_name -P proc \t",argv[0]); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  qdp_lsum2(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * lsum2(out, vec1, n2vec)
   * 
   * fpoint vec1[nvec][Ncol][rei]
   * asmint n3vec
   *
   * for(i=0;i<nvec;i++)
   *   out[i] =a* vec1[i];
   *
   */

void qdp_lsum2( char *name)
{
  int dum = defargcount(3);

  /*Integer register usage*/
  /* Pointers to vec1, vec2, vec3, a and  b*/
  /* and the counter */
  alreg(sumptr,Iregs);  /* Pointer to the sum */
  alreg(vec1ptr,Iregs); /* Pointer to input vector */
  alreg(counter,Iregs);
  alreg(zero, Iregs); /* Used to zero sum */
  alreg(one, Iregs);

  /*Floating register usage*/
  reg_array_2d(vec1,Fregs,3,2); /* An array to hold an atom of input */
  reg_array_2d(sums,Fregs,3,2);

  int Isize = PROC->I_size;
  def_off(word_size,Byte,Isize);
  
  /* Various offsets. Start with 0 */
  def_dp_off(ZERO,0);

  /* A vector is of length 6 */
  def_dp_off(VEC_ATOM,6);

  /* This I am guessing is some pattern to describe the vectors
     and that it has to match the 2d register allocation above */
  offset_2d_dp(VEC_IMM,3,2);

  /* Declare memory streams for prefetching */
  struct stream *PreVec1;

  int brchno,retno; /*Branch target handles*/
  int co,rei;

  /* Boilerplate */
  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /* Map the arguments -- in order */
  getarg(sumptr); 
  getarg(vec1ptr);
  getarg(counter);


  /* Ugly hack to get floating zero */
  queue_iload_imm(zero, ZERO);
  queue_istore(zero, ZERO, sumptr);
  if( PROC->I_size < PROC->FP_size ) { 
    queue_istore(zero, word_size, sumptr);
  }


  /*
   * Insert a label to prevent reordering
   */
  make_inst(DIRECTIVE,Target,get_target_label());

  
  

  /* Create the streams in and out */
  PreVec1= create_stream(VEC_ATOM,vec1ptr ,counter,STREAM_IN ,LINEAR);

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  for(co=0; co < 3; co++) { 
    for(rei=0; rei < 2; rei++) { 
      queue_fload(sums[co][rei],ZERO,sumptr);
    }
  }     

  /*
   * Start software pipeline
   */
 
  for(co=0; co < 3; co++) { 
    for (rei=0;rei<2;rei++){
      queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    }
  }

  iterate_stream(PreVec1);
  queue_prefetch(PreVec1);

  /* Disgusting hack to reduce the counter by 1 */
  def_off(ONE,Byte,1);
  make_inst(LOADPIPE,LOAD_IMM,one,ONE);
  make_inst(IALUPIPE, ISUB, counter, counter, one);

  brchno = start_loop(counter);

  pragma(DCBT_SPACE,5);

  /* First color component */
  for(co=0; co < 3; co++) { 
    for(rei=0; rei < 2; rei++) { 
      queue_fmadd(sums[co][rei],vec1[co][rei],vec1[co][rei],sums[co][rei]);
    }
  }

  for(co=0; co < 3; co++) {
    for (rei=0;rei<2;rei++){
      queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    }
  }
  iterate_stream(PreVec1);
  queue_prefetch(PreVec1);

  stop_loop(brchno,counter);

  /* First color component */
  for(co=0; co < 3; co++) { 
    for(rei=0; rei < 2; rei++) { 
      queue_fmadd(sums[co][rei],vec1[co][rei],vec1[co][rei],sums[co][rei]);
    }
  }
  
  for(co=0; co < 3; co++) { 
    queue_fadd(sums[co][0], sums[co][0], sums[co][1]);
  }
  queue_fadd(sums[0][0], sums[0][0], sums[1][0]);
  queue_fadd(sums[0][0], sums[0][0], sums[2][0]);

  queue_fstore(sums[0][0], ZERO, sumptr);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








