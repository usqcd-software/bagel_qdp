
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

#include "registers.h"


void qdp_vscal3( char *);

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
  qdp_vscal3(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * vscal(out,  a, vec1 , n3vec )
   * fpoint a
   * fpoint vec1[nvec][Ncol][rei]
   * asmint n3vec
   *
   * for(i=0;i<nvec;i++)
   *   out[i] =a* vec1[i];
   *
   */

void qdp_vscal3( char *name)
{
  int dum = defargcount(4);

  /*Integer register usage*/
  /* Pointers to vec1, vec2, vec3, a and  b*/
  /* and the counter */
  alreg(vec1ptr,Iregs);
  alreg(vec3ptr,Iregs);
  alreg(counter,Iregs);
  alreg(outptr,Iregs);
  alreg(Aptr, Iregs); /* Pointer to A */
  /*Floating register usage*/

  alreg(A, Fregs); /* A itself */
  reg_array_2d(vec1,Fregs,3,2);
  reg_array_2d(vec3,Fregs,3,2);

  /* Various offsets. Start with 0 */
  def_dp_off(ZERO,0);

  /* A vector is of length 6 */
  def_dp_off(VEC_ATOM,6);

  /* This I am guessing is some pattern to describe the vectors
     and that it has to match the 2d register allocation above */
  offset_2d_dp(VEC_IMM,3,2);

  /* Declare memory streams for prefetching */
  struct stream *PreVec1;
  struct stream *PreVec3;

  int brchno,retno; /*Branch target handles*/
  int co,rei;

  /* Boilerplate */
  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /* Map the arguments -- in order */
  getarg(vec3ptr); 
  getarg(Aptr);
  getarg(vec1ptr);
  getarg(counter);

  /* get A into A from Aptr */
  queue_fload(A, ZERO, Aptr);

  /* Create the streams in and out */
  PreVec1= create_stream(VEC_ATOM,vec1ptr ,counter,STREAM_IN ,LINEAR);
  PreVec3= create_stream(VEC_ATOM,vec3ptr,counter,STREAM_OUT,LINEAR);
  
  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */

  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
  }

  brchno = start_loop(counter);

  pragma(DCBT_SPACE,5);

  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,vec3ptr,ZERO);
  
  /* First color component */
  co = 0;
  for(rei=0;rei<2;rei++){
    queue_fmul(vec3[co][rei],A,vec1[co][rei]);
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr);
  }

  /* Second color component */
  co = 1;
  for(rei=0;rei<2;rei++){
    queue_fmul(vec3[co][rei],A,vec1[co][rei]);
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr);
  }

  co = 2;
  for(rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    queue_fmul(vec3[co][rei],A,vec1[co][rei]);
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr);
  }


  iterate_stream(PreVec1);
  iterate_stream(PreVec3);

  queue_prefetch(PreVec1);


  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);

  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
  }


  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








