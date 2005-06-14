
/*
 *
 *  Copyright UKQCD Collaboration, November 2000.
 *  Written by Peter Boyle.
 *  This software is provided for NON-COMMERCIAL use only,
 *  and may not be redistributed without permission.
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "processor.h"
}

#include "registers.h"


void qdp_vaxpby( char *);

int overwrite1 = 0;
int overwrite2 = 0;
int do_vaxmy = 0;

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
    case 'm': do_vaxmy = 1; break;
    default: fprintf(stderr,"Usage: %s -nroutine_name -P proc -m\t",argv[0]); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  qdp_vaxpby(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * vaxpby (out,  a, vec1, b, vec2 , nvec )
   * fpoint a
   * fpoint vec1[nvec][Ncol][rei]
   * fpoint b
   * fpoint vec2[nvec][Ncol][rei] 
   * asmint nvec
   *
   * for(i=0;i<nvec;i++)
   *   vec1[i] = a * vec1[i] + b*vec2 [i];
   *
   */

void qdp_vaxpby( char *name)
{
  int dum = defargcount(6);

  /*Integer register usage*/
  /* Pointers to vec1, vec2, vec3, a and  b*/
  /* and the counter */
  alreg(vec1ptr,Iregs);
  alreg(vec2ptr,Iregs);
  alreg(vec3ptr,Iregs);
  alreg(counter,Iregs);
  alreg(Aptr,Iregs);
  alreg(Bptr,Iregs);
  alreg(outptr,Iregs);

  /*Floating register usage*/
  /* Alltogether 3*2*3+2 = 20 registers */
  alreg(A,Fregs);
  alreg(B,Fregs);
  alreg(tmp,Fregs);

  reg_array_2d(vec1,Fregs,3,2);
  reg_array_2d(vec2,Fregs,3,2);
  reg_array_2d(vec3,Fregs,3,2);

  /* Various offsets. Start with 0 */
  def_off(ZERO,0);

  /* A vector is of length 6 */
  def_off(VEC_ATOM,6);

  /* This I am guessing is some pattern to describe the vectors
     and that it has to match the 2d register allocation above */
  offset_2d(VEC_IMM,3,2);

  /* Declare memory streams for prefetching */
  struct stream *PreVec1;
  struct stream *PreVec2;
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
  getarg(Bptr);
  getarg(vec2ptr);           
  getarg(counter);

  /* Actually these are reals not complex-s*/
  queue_fload(A,ZERO,Aptr); /*Find the complex scale factor*/
  queue_fload(B,ZERO,Bptr); /*Find the other scle factor */

  /* Create the streams in and out */
  PreVec1= create_stream(VEC_ATOM,vec1ptr ,counter,STREAM_IN ,LINEAR);
  PreVec2= create_stream(VEC_ATOM,vec2ptr ,counter,STREAM_IN ,LINEAR);
  PreVec3 = create_stream(VEC_ATOM,vec3ptr,counter,STREAM_OUT,LINEAR);
  
  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */

  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);

  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);

  }

  brchno = start_loop(counter);

  pragma(DCBT_SPACE,5);

  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,vec3ptr,ZERO);
  
  /* First color component */
  co = 0;
  for(rei=0;rei<2;rei++){
    queue_fmul(tmp,B,vec2[co][rei]);
    if ( !do_vaxmy ) { 
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],tmp);
    } else {
      make_inst(FMACPIPE,FMSUB,vec3[co][rei],A,vec1[co][rei],tmp);
    }
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr);
  }

  /* Second color component */
  co = 1;
  for(rei=0;rei<2;rei++){
    queue_fmul(tmp,B,vec2[co][rei]);
    if ( !do_vaxmy ) { 
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],tmp);
    } else {
      make_inst(FMACPIPE,FMSUB,vec3[co][rei],A,vec1[co][rei],tmp);
    }
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr);
  }

  co = 2;
  for(rei=0;rei<2;rei++){
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);

    queue_fmul(tmp,B,vec2[co][rei]);
    if ( !do_vaxmy ) { 
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],tmp);
    } else {
      make_inst(FMACPIPE,FMSUB,vec3[co][rei],A,vec1[co][rei],tmp);
    }
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr);
  }


  iterate_stream(PreVec1);
  iterate_stream(PreVec2);
  iterate_stream(PreVec3);

  queue_prefetch(PreVec2);
  queue_prefetch(PreVec1);


  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);

  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
  }


  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








