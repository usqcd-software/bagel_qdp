
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


void qdp_vadd3( char *);

int do_minus = 0;

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
    case 'm': do_minus = 1; break;
    default: fprintf(stderr,"Usage: %s -nroutine_name -P proc -m\t",argv[0]); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  qdp_vadd3(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * vadd (out,  vec1, vec2 , n3vec )
   * fpoint out [nvec][Ncol][rei]
   * fpoint vec1[nvec][Ncol][rei]
   * fpoint vec2[nvec][Ncol][rei] 
   * asmint n3vec
   *
   * for(i=0;i<nvec;i++)
   *   out[i] = vec1[i] +/- vec2 [i];
   *
   */

void qdp_vadd3( char *name)
{
  int dum = defargcount(4);

  /*Integer register usage*/
  /* Pointers to vec1, vec2, vec3, a and  b*/
  /* and the counter */
  alreg(vec1ptr,Iregs);
  alreg(vec2ptr,Iregs);
  alreg(vec3ptr,Iregs);
  alreg(counter,Iregs);
  alreg(outptr,Iregs);
  alreg(tmp, Iregs);
  alreg(one, Iregs);

  /*Floating register usage*/
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
  getarg(vec1ptr);
  getarg(vec2ptr);           
  getarg(counter);

  /* Create the streams in and out */
  PreVec1= create_stream(VEC_ATOM,vec1ptr ,counter,STREAM_IN ,LINEAR);
  PreVec2= create_stream(VEC_ATOM,vec2ptr ,counter,STREAM_IN ,LINEAR);
  PreVec3= create_stream(VEC_ATOM,vec3ptr,counter,STREAM_OUT,LINEAR);
  
  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */
  pragma(DCBT_SPACE,6);

  queue_fload(vec1[0][0],VEC_IMM[0][0],vec1ptr);
  queue_fload(vec1[0][1],VEC_IMM[0][1],vec1ptr);
  queue_fload(vec1[1][0],VEC_IMM[1][0],vec1ptr);
  queue_fload(vec1[1][1],VEC_IMM[1][1],vec1ptr);
  queue_fload(vec1[2][0],VEC_IMM[2][0],vec1ptr);
  queue_fload(vec1[2][1],VEC_IMM[2][1],vec1ptr);
  iterate_stream(PreVec1);
  queue_prefetch(PreVec1);

  queue_fload(vec2[0][0],VEC_IMM[0][0],vec2ptr);
  queue_fload(vec2[0][1],VEC_IMM[0][1],vec2ptr);
  queue_fload(vec2[1][0],VEC_IMM[1][0],vec2ptr);
  queue_fload(vec2[1][1],VEC_IMM[1][1],vec2ptr);
  queue_fload(vec2[2][0],VEC_IMM[2][0],vec2ptr);
  queue_fload(vec2[2][1],VEC_IMM[2][1],vec2ptr);
  iterate_stream(PreVec2);
  queue_prefetch(PreVec2);


  /* counter - 1 loop */
  /* last bit is done at end */
  /* Peter doesn't provide a queue_isub so I do it this way */
  /* Disgusting hack to reduce the counter by 1 */
  int ONE=def_byte_offset(1,"ONE");
  make_inst(LOADPIPE,LOAD_IMM,one,ONE);
  make_inst(IALUPIPE, ISUB, counter, counter, one);

  brchno = start_loop(counter);



  /* This sets the outptr to the vecptr - with an add immediate */
  queue_iadd_imm(outptr,vec3ptr,ZERO);

  if ( !do_minus ) { 
    queue_fadd(vec3[0][0],vec1[0][0],vec2[0][0]);
  } else {
    queue_fsub(vec3[0][0],vec1[0][0],vec2[0][0]);
  }
  if ( !do_minus ) { 
    queue_fadd(vec3[0][1],vec1[0][1],vec2[0][1]);
  } else {
    queue_fsub(vec3[0][1],vec1[0][1],vec2[0][1]);
  }


  if ( !do_minus ) { 
    queue_fadd(vec3[1][0],vec1[1][0],vec2[1][0]);
  } else {
    queue_fsub(vec3[1][0],vec1[1][0],vec2[1][0]);
  }
  if ( !do_minus ) { 
    queue_fadd(vec3[1][1],vec1[1][1],vec2[1][1]);
  } else {
    queue_fsub(vec3[1][1],vec1[1][1],vec2[1][1]);
  }


  if ( !do_minus ) { 
    queue_fadd(vec3[2][0],vec1[2][0],vec2[2][0]);
  } else {
    queue_fsub(vec3[2][0],vec1[2][0],vec2[2][0]);
  }

  if ( !do_minus ) { 
    queue_fadd(vec3[2][1],vec1[2][1],vec2[2][1]);
  } else {
    queue_fsub(vec3[2][1],vec1[2][1],vec2[2][1]);
  }
  

  queue_fstore(vec3[0][0],VEC_IMM[0][0],outptr);
  queue_fstore(vec3[0][1],VEC_IMM[0][1],outptr);
  queue_fstore(vec3[1][0],VEC_IMM[1][0],outptr);
  queue_fstore(vec3[1][1],VEC_IMM[1][1],outptr);
  queue_fstore(vec3[2][0],VEC_IMM[2][0],outptr);
  queue_fstore(vec3[2][1],VEC_IMM[2][1],outptr);

  iterate_stream(PreVec3);
  
  queue_fload(vec1[0][0],VEC_IMM[0][0],vec1ptr);
  queue_fload(vec1[0][1],VEC_IMM[0][1],vec1ptr);
  queue_fload(vec1[1][0],VEC_IMM[1][0],vec1ptr);
  queue_fload(vec1[1][1],VEC_IMM[1][1],vec1ptr);
  queue_fload(vec1[2][0],VEC_IMM[2][0],vec1ptr);
  queue_fload(vec1[2][1],VEC_IMM[2][1],vec1ptr);
  iterate_stream(PreVec1);
  queue_prefetch(PreVec1);

  queue_fload(vec2[0][0],VEC_IMM[0][0],vec2ptr);
  queue_fload(vec2[0][1],VEC_IMM[0][1],vec2ptr);
  queue_fload(vec2[1][0],VEC_IMM[1][0],vec2ptr);
  queue_fload(vec2[1][1],VEC_IMM[1][1],vec2ptr);
  queue_fload(vec2[2][0],VEC_IMM[2][0],vec2ptr);
  queue_fload(vec2[2][1],VEC_IMM[2][1],vec2ptr);
  iterate_stream(PreVec2);
  queue_prefetch(PreVec2);

  stop_loop(brchno,counter);

  queue_iadd_imm(outptr,vec3ptr,ZERO);
  
  if ( !do_minus ) { 
    queue_fadd(vec3[0][0],vec1[0][0],vec2[0][0]);
  } else {
    queue_fsub(vec3[0][0],vec1[0][0],vec2[0][0]);
  }
  if ( !do_minus ) { 
    queue_fadd(vec3[0][1],vec1[0][1],vec2[0][1]);
  } else {
    queue_fsub(vec3[0][1],vec1[0][1],vec2[0][1]);
  }


  if ( !do_minus ) { 
    queue_fadd(vec3[1][0],vec1[1][0],vec2[1][0]);
  } else {
    queue_fsub(vec3[1][0],vec1[1][0],vec2[1][0]);
  }
  if ( !do_minus ) { 
    queue_fadd(vec3[1][1],vec1[1][1],vec2[1][1]);
  } else {
    queue_fsub(vec3[1][1],vec1[1][1],vec2[1][1]);
  }

  if ( !do_minus ) { 
    queue_fadd(vec3[2][0],vec1[2][0],vec2[2][0]);
  } else {
    queue_fsub(vec3[2][0],vec1[2][0],vec2[2][0]);
  }
  if ( !do_minus ) { 
    queue_fadd(vec3[2][1],vec1[2][1],vec2[2][1]);
  } else {
    queue_fsub(vec3[2][1],vec1[2][1],vec2[2][1]);
  }


  queue_fstore(vec3[0][0],VEC_IMM[0][0],outptr);
  queue_fstore(vec3[0][1],VEC_IMM[0][1],outptr);

  queue_fstore(vec3[1][0],VEC_IMM[1][0],outptr);
  queue_fstore(vec3[1][1],VEC_IMM[1][1],outptr);

  queue_fstore(vec3[2][0],VEC_IMM[2][0],outptr);
  queue_fstore(vec3[2][1],VEC_IMM[2][1],outptr);
 
 

  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








