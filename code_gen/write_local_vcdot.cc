
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


void qdp_lcdot( char *);

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
  qdp_lcdot(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * lcdot(out_re, out_im,  vec1, vec2,  n3vec)
   * fpoint *out_re,
   * fpoint *out_im;
   *
   * fpoint vec1[nvec][Ncol][rei]
   * fpoint vec2[nvec][Ncol][rei]
   * asmint n3vec
   *
   * out.re=0;
   * out.im=0;
   * for(i=0;i<nvec;i++) {
   *   out.re += vec1[i].re*vec2[i].re + vec1[i].im*vec2[i].im;
   *   out_im += vec1[i].re*vec2[i].im - vec1[i].im*vec2[i].re;
   * }
   */

void qdp_lcdot( char *name)
{
  int dum = defargcount(5);

  /*Integer register usage*/
  /* Pointers to vec1, vec2, vec3, a and  b*/
  /* and the counter */
  alreg(outreptr,Iregs);  /* Pointer to the real part of output */
  alreg(outimptr,Iregs);  /* Pointer to imag part of output */
  alreg(vec1ptr,Iregs);   /* Pointer to input vector1 */
  alreg(vec2ptr,Iregs);   /* Pointer to input vector2 */
  alreg(counter,Iregs);   /* Pointer to counter */
  alreg(tmp, Iregs); /* Used to zero things */

  /*Floating register usage - altogether 26 fp registers.. */
  alreg(outre,Fregs);     /* The real part of the cdot */
  alreg(outim,Fregs);     /* The imag part of the cdot */
  reg_array_2d(vec1,Fregs,3,2); /* An array to hold an atom of vec 1 */
  reg_array_2d(vec2,Fregs,3,2); /* An array to hold an atom of vec 2 */
  reg_array_2d(outs,Fregs,3,2); /* An array to hold partial accumulators */
                             
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
  getarg(outreptr); 
  getarg(outimptr);
  getarg(vec1ptr);
  getarg(vec2ptr);
  getarg(counter);


  /* Ugly hack to get floating zero */
  queue_iload_imm(tmp, ZERO);
  queue_istore(tmp, ZERO, outreptr);
  queue_istore(tmp, ZERO, outimptr);
  if( PROC->I_size < PROC->FP_size ) { 
    queue_istore(tmp, word_size, outreptr);
    queue_istore(tmp, word_size, outimptr);
  }

  /*
   * Insert a label to prevent reordering
   */
  make_inst(DIRECTIVE,Target,get_target_label());  

  /* Create the streams in and out */
  PreVec1= create_stream(VEC_ATOM,vec1ptr ,counter,STREAM_IN, LINEAR);
  PreVec2= create_stream(VEC_ATOM,vec2ptr ,counter,STREAM_IN, LINEAR);

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /* Load the "zero answer" into outres[0][re] */
  queue_fload(outs[0][0],ZERO,outreptr);
  /* Copy it into outres[0][im] and outims[0][re,im] */
  queue_fmov(outs[0][1], outs[0][0]);
  
  /* Now do the other co-s */
  for(co=1; co < 3; co++) { 
    for(rei=0; rei < 2; rei++) { 
      queue_fmov(outs[co][rei],outs[0][0]);
    }
  }     

  /* The accumulators should now hold zeros */

  /*
   * Start software pipeline
   */
 
  co=0;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
  }

  co=1;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
  }
  

  brchno = start_loop(counter);

  pragma(DCBT_SPACE,5);

  /* First color component */
  co=0;
  // out_re += v1.r*v2.r
  queue_fmadd(outs[co][0], vec1[co][0],vec2[co][0], outs[co][0]);
  // out_im += v1.r*v2.i
  queue_fmadd(outs[co][1], vec1[co][0],vec2[co][1], outs[co][1]);
  // out_re += v1.i*v2.i
  queue_fmadd(outs[co][0], vec1[co][1],vec2[co][1], outs[co][0]);
  // out_im -= v1.i*v2.re
  queue_fnmsub(outs[co][1], vec1[co][1],vec2[co][0], outs[co][1]);

  co=1;
  // out_re += v1.r*v2.r
  queue_fmadd(outs[co][0], vec1[co][0],vec2[co][0], outs[co][0]);
  // out_im += v1.r*v2.i
  queue_fmadd(outs[co][1], vec1[co][0],vec2[co][1], outs[co][1]);
  // out_re += v1.i*v2.i
  queue_fmadd(outs[co][0], vec1[co][1],vec2[co][1], outs[co][0]);
  // out_im -= v1.i*v2.re
  queue_fnmsub(outs[co][1], vec1[co][1],vec2[co][0], outs[co][1]);

  co=2;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
  }
  
  // out_re += v1.r*v2.r
  queue_fmadd(outs[co][0], vec1[co][0],vec2[co][0], outs[co][0]);
  // out_im += v1.r*v2.i
  queue_fmadd(outs[co][1], vec1[co][0],vec2[co][1], outs[co][1]);
  // out_re += v1.i*v2.i
  queue_fmadd(outs[co][0], vec1[co][1],vec2[co][1], outs[co][0]);
  // out_im -= v1.i*v2.re
  queue_fnmsub(outs[co][1], vec1[co][1],vec2[co][0], outs[co][1]);
  
  iterate_stream(PreVec1);
  iterate_stream(PreVec2);

  queue_prefetch(PreVec1);
  queue_prefetch(PreVec2);

  co=0;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
  }

  co=1;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
  }
 
  stop_loop(brchno,counter);

  queue_fmov(outre, outs[0][0]);
  queue_fadd(outre, outre,outs[1][0]);
  queue_fadd(outre, outre,outs[2][0]);

  queue_fmov(outim, outs[0][1]);
  queue_fadd(outim, outim, outs[1][1]);
  queue_fadd(outim, outim, outs[2][1]);


  queue_fstore(outre, ZERO, outreptr);
  queue_fstore(outim, ZERO, outimptr);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








