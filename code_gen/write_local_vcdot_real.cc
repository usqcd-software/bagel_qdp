 
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
  extern struct processor* PROC;
}

#include "registers.h"


void qdp_lcdotr( char *);

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
  qdp_lcdotr(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * lcdot(out,  vec1, vec2,  n3vec)
   * fpoint *out,
   *
   * fpoint vec1[nvec][Ncol][rei]
   * fpoint vec2[nvec][Ncol][rei]
   * asmint n3vec
   *
   * out.re=0;
   * for(i=0;i<nvec;i++) {
   *   out.re += vec1[i].re*vec2[i].re + vec1[i].im*vec2[i].im;
   * }
   */

void qdp_lcdotr( char *name)
{
  int dum = defargcount(4);

  /*Integer register usage*/
  /* Pointers to vec1, vec2, vec3, a and  b*/
  /* and the counter */
  alreg(outptr,Iregs);  /* Pointer to the real part of output */
  alreg(vec1ptr,Iregs);   /* Pointer to input vector1 */
  alreg(vec2ptr,Iregs);   /* Pointer to input vector2 */
  alreg(counter,Iregs);   /* Pointer to counter */
  alreg(tmp, Iregs); /* Used to zero things */
  alreg(one, Iregs);

  /*Floating register usage */
  reg_array_2d(vec1,Fregs,3,2); /* An array to hold an atom of vec 1 */
  reg_array_2d(vec2,Fregs,3,2); /* An array to hold an atom of vec 2 */
  /*  alreg(out,Fregs);
      alreg(out1,Fregs);
  */
  reg_array_2d(out,Fregs,3,2);
  alreg(zero,Fregs);

  /* I will want to zero out the result (at outreptr and outimptr) */
  /* I need the word_size for that ... */
  int Isize = PROC->I_size;
  int word_size = def_byte_offset(Isize,"word_size");
  
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

  int brchno,retno; /*Branch target handles*/
  int co,rei;

  /* Boilerplate */
  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /* Map the arguments -- in order */
  getarg(outptr); 
  getarg(vec1ptr);
  getarg(vec2ptr);
  getarg(counter);

  /* Ugly hack to get floating zero */
  queue_iload_imm(tmp, ZERO);
  queue_istore(tmp, ZERO, outptr);
  if( PROC->I_size < PROC->FP_size ) { 
    queue_istore(tmp, word_size, outptr);
  }
  /* Load zero value back into outptr */
  for(co=0; co < 3; co++) { 
    for(rei=0; rei < 2 ; rei++) { 
      queue_fload(out[co][rei], ZERO, outptr);
    }
  }
  
  /*
   * Insert a label to prevent reordering
   */
  make_inst(DIRECTIVE,Target,get_target_label());  


  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  // pragma(DCBT_SPACE,0);
  /* Create the streams in and out */
  PreVec1= create_stream(VEC_ATOM,vec1ptr ,counter,STREAM_IN, LINEAR);
  PreVec2= create_stream(VEC_ATOM,vec2ptr ,counter,STREAM_IN, LINEAR);

  /*
   * Start software pipeline
   */
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

  /* 
     result += v1_0r * v2_0r;
  */
  queue_fmadd(out[0][0], vec1[0][0], vec2[0][0],out[0][0]);

  /* result = result + v1_0i*v2_0i; */

  queue_fmadd(out[0][1], vec1[0][1], vec2[0][1], out[0][1]);

  /* result = result + v1_1r*v2_1r; */
  queue_fmadd(out[1][0], vec1[1][0], vec2[1][0], out[1][0]);

  /* result = result + v1_1i*v2_1i; */
  queue_fmadd(out[1][1], vec1[1][1], vec2[1][1], out[1][1]);

  /* result = result + v1_2r*v2_2r; */
  queue_fmadd(out[2][0], vec1[2][0], vec2[2][0], out[2][0]);

  /* 
      result = result + v1_2i*v2_2i;
  */
  queue_fmadd(out[2][1], vec1[2][1], vec2[2][1], out[2][1]);

  

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

  
  queue_fmadd(out[0][0], vec1[0][0], vec2[0][0] ,out[0][0]);
  queue_fmadd(out[0][1], vec1[0][1], vec2[0][1] ,out[0][1]);
  queue_fmadd(out[1][0], vec1[1][0], vec2[1][0] ,out[1][0]);
  queue_fmadd(out[1][1], vec1[1][1], vec2[1][1] ,out[1][1]);
  queue_fmadd(out[2][0], vec1[2][0], vec2[2][0] ,out[2][0]);
  queue_fmadd(out[2][1], vec1[2][1], vec2[2][1] ,out[2][1]);

  queue_fadd(out[0][0], out[0][0], out[0][1]);
  queue_fadd(out[1][0], out[1][0], out[1][1]);
  queue_fadd(out[2][0], out[2][0], out[2][1]);
  queue_fadd(out[0][0], out[0][0], out[1][0]);
  queue_fadd(out[0][0], out[0][0], out[2][0]);

  make_inst(DIRECTIVE,Target,retno);  
  queue_fstore(out[0][0], ZERO, outptr);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








