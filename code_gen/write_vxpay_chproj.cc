
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


void qdp_vxpay_chp( char *);
void qdp_vxpay_chm( char *);


char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/

static int do_minus=0;

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
    case 'm':
      {
	do_minus=1;
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
    qdp_vxpay_chp(name);
  }
  else { 
    fprintf(stderr, "Duing chiral proj minus\n");
    qdp_vxpay_chm(name); 
  }

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

/* qdp_vxpay_chp(y, a, x, y, n_4vec)
 *
 * Operation: z = x +/- (1/2)(1 + g5)a * y
 *
 *  g5 is in chiral basis diag(1,1,-1, -1)
 * (1/2)(1+g5) = diag(1,1,0,0)
 *
 * z is first element of output vector (Float *)
 * x is first element of unprojected input vector (Float *)
 * y is first element of projected input  vector (Float *)
 * n_4vec is the number of 4 vectors
 *
 * for(i=0; i < n_4vec; i++) {
 *   for(spin=0; spin < 2: spin++) { 
 *     for(col=0; col < Ncol; col++) {
 *        z[i][spin][col] = x[i][spin][col] +/- a * y[i][spin][col];
 *     }
 *   }
 *   for(spin=2; spin < Nspin; spin++) {
 *     for(col=0; col < Ncol; col++) { 
 *       z[i][spin][col] = x[i][spin][col]
 *     }
 *   }
 * }
 */

void qdp_vxpay_chp( char *name)
{
  int dum = defargcount(5);

  /*Integer register usage*/
  /* Pointers to vec1, vec2, vec3, a and  b*/
  /* and the counter */
  alreg(xptr,Iregs);
  alreg(yptr,Iregs);
  alreg(zptr,Iregs);
  alreg(aptr,Iregs);
  alreg(counter,Iregs);
  alreg(outptr,Iregs);

  /*Floating register usage*/
  alreg(a, Fregs);
  reg_array_2d(x,Fregs,3,2);
  reg_array_2d(y,Fregs,3,2);
  reg_array_2d(z,Fregs,3,2);

  /* Various offsets. Start with 0 */
  def_off(ZERO,0);

  /* A vector is of length 6 */
  def_off(VEC_ATOM,6);

  /* This I am guessing is some pattern to describe the vectors
     and that it has to match the 2d register allocation above */
  offset_2d(VEC_IMM,3,2);

  /* Declare memory streams for prefetching */
  struct stream *PreVecX;
  struct stream *PreVecY;
  struct stream *PreVecZ;

  int brchno,retno; /*Branch target handles*/
  int co,rei;



  /* Boilerplate */
  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /* Map the arguments -- in order */
  getarg(zptr);
  getarg(aptr);
  getarg(xptr); 
  getarg(yptr);
  getarg(counter);


  /* Create the streams in and out */
  /* A 4spinor contians 4 atoms of length 3x2 */
  PreVecX= create_stream(VEC_ATOM,xptr, 4*counter,STREAM_IN, LINEAR);
  PreVecY= create_stream(VEC_ATOM,yptr, 4*counter,STREAM_IN, LINEAR);
  PreVecZ= create_stream(VEC_ATOM,zptr, 4*counter,STREAM_OUT, LINEAR);

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  queue_fload(a, ZERO, aptr);
  /*
   * Start software pipeline
   */
  /******************* START LOOP **********************/
  brchno = start_loop(counter);

  pragma(DCBT_SPACE,5);

  /* ************* PRELOAD SPIN COMPONENT 0 *******************/
  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }

  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  /********** SPIN COMPONENT 0 ************************/
  /* Point outptr to base of current zptr atom for storing */
  queue_iadd_imm(outptr,zptr,ZERO);

  co = 0;
  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  co = 1;
  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  co = 2;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
 
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);

  iterate_stream(PreVecY);
  queue_prefetch(PreVecY);

  iterate_stream(PreVecZ);
  
  /*************** PRELOAD SPIN COMPONENT 1 *****************/
  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }

  co=0;
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  co=1;
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  /************* DO SPIN COMPONENT 1 ************************/
  queue_iadd_imm(outptr,zptr,ZERO);

  co = 0;
  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  co = 1;
  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  co = 2;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
 
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);

  iterate_stream(PreVecY);
  iterate_stream(PreVecZ);


  /* Preload X for spin component 2 ***********************/
  co=0;
  for(rei=0; rei<2;rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }
  co=1;
  for(rei=0; rei<2;rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }
  co=2;
  for(rei=0; rei<2; rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }

  /**************** SPIN COMPONENT 2 **********************/
  queue_iadd_imm(outptr,zptr,ZERO);
  co=0;
  for(rei=0; rei<2;rei++) {
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);
  }
  co=1;
  for(rei=0; rei<2; rei++) { 
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);
  }
  co=2;
  for(rei=0; rei<2; rei++) { 
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);
  }

  /* Iterate only once as we use Zptr for the outptr thing */
  iterate_stream(PreVecZ);

  /* Iterate once and queue a prefetch as we need in next loop */
  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);

  /* Iterate twice */
  iterate_stream(PreVecY); 
  iterate_stream(PreVecY);

  /* Preload X for spin component 3 ***********************/
  co=0;
  for(rei=0; rei<2;rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }
  co=1;
  for(rei=0; rei<2;rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }
  co=2;
  for(rei=0; rei<2; rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }

  /**************** SPIN COMPONENT 3 **********************/
  queue_iadd_imm(outptr,zptr,ZERO);
  co=0;
  for(rei=0; rei<2;rei++) {
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);
  }
  co=1;
  for(rei=0; rei<2; rei++) { 
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);
  }
  co=2;
  for(rei=0; rei<2; rei++) { 
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);
  }
  iterate_stream(PreVecZ);
  queue_prefetch(PreVecY);
  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);
  // Y already iterated and prefetch already stored

  /***************** END LOOP *****************************/
  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}

/* qdp_vxpay_chm(z, a, x, y, n_4vec)
 *
 * Operation: y = x +/-  (1/2)(1 - g5)a* x
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
 *   for(spin=0; spin < 2; spin++) {
 *     for(col=0; col < Ncol; col++) { 
 *       z[i][spin][col] = x[i][spin][col] +/- a*y[i][spin][col];
 *     }
 *   }
 *   for(spin=2; spin < Nspin: spin++) { 
 *     for(col=0; col < Ncol; col++) {
 *        y[i][spin][col] = x[i][spin][col];
 *     }
 *   }
 * }
 */

void qdp_vxpay_chm( char *name)
{
  int dum = defargcount(5);

  /*Integer register usage*/
  /* Pointers to vec1, vec2, vec3, a and  b*/
  /* and the counter */
  alreg(xptr,Iregs);
  alreg(yptr,Iregs);
  alreg(zptr,Iregs);
  alreg(aptr,Iregs);

  alreg(counter,Iregs);
  alreg(outptr,Iregs);

  /*Floating register usage*/
  reg_array_2d(x,Fregs,3,2);
  reg_array_2d(y,Fregs,3,2);
  reg_array_2d(z,Fregs,3,2);
  alreg(a, Fregs);
  /* Various offsets. Start with 0 */
  def_off(ZERO,0);

  /* A vector is of length 6 */
  def_off(VEC_ATOM,6);

  /* This I am guessing is some pattern to describe the vectors
     and that it has to match the 2d register allocation above */
  offset_2d(VEC_IMM,3,2);

  /* Declare memory streams for prefetching */
  struct stream *PreVecX;
  struct stream *PreVecY;
  struct stream *PreVecZ;

  int brchno,retno; /*Branch target handles*/
  int co,rei;



  /* Boilerplate */
  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  /* Map the arguments -- in order */
  getarg(zptr);
  getarg(aptr);
  getarg(xptr); 
  getarg(yptr);
  getarg(counter);


  /* Create the streams in and out */
  /* A 4spinor contians 4 atoms of length 3x2 */
  PreVecX= create_stream(VEC_ATOM,xptr, 4*counter,STREAM_IN, LINEAR);
  PreVecY= create_stream(VEC_ATOM,yptr, 4*counter,STREAM_IN, LINEAR);
  PreVecZ= create_stream(VEC_ATOM,zptr, 4*counter,STREAM_OUT, LINEAR);

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  queue_fload(a, ZERO, aptr);

  /*
   * Start software pipeline
   */
  /******************* START LOOP **********************/
  brchno = start_loop(counter);

  iterate_stream(PreVecY);
  iterate_stream(PreVecY);

  pragma(DCBT_SPACE,6);
  /* Preload X for spin component 0 ***********************/
  co=0;
  for(rei=0; rei<2;rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }
  co=1;
  for(rei=0; rei<2;rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }
  co=2;
  for(rei=0; rei<2; rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }

  /**************** SPIN COMPONENT 0 **********************/
  queue_iadd_imm(outptr,zptr,ZERO);
  co=0;
  for(rei=0; rei<2;rei++) {
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);
  }
  co=1;
  for(rei=0; rei<2; rei++) { 
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);
  }
  co=2;
  for(rei=0; rei<2; rei++) { 
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);
  }

  /* Iterate only once as we use Zptr for the outptr thing */
  iterate_stream(PreVecZ);

  /* Iterate once and queue a prefetch as we need in next loop */
  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);

  /* Preload X for spin component 1 ***********************/
  co=0;
  for(rei=0; rei<2;rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }
  co=1;
  for(rei=0; rei<2;rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }
  co=2;
  for(rei=0; rei<2; rei++) { 
    queue_fload(x[co][rei], VEC_IMM[co][rei],xptr);
  }

  /**************** SPIN COMPONENT 1 **********************/
  queue_iadd_imm(outptr,zptr,ZERO);
  co=0;
  for(rei=0; rei<2;rei++) {
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);    
  }
  co=1;
  for(rei=0; rei<2; rei++) { 
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);
  }
  co=2;
  for(rei=0; rei<2; rei++) { 
    queue_fstore(x[co][rei], VEC_IMM[co][rei], outptr);
  }

  /* queue a prefetch for Y. It is already iterated */
  queue_prefetch(PreVecY);

  /* Iterate only once as we use Zptr for the outptr thing */
  iterate_stream(PreVecZ);


  /* Iterate once and queue a prefetch as we need in next loop */
  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);


  /* ************* PRELOAD SPIN COMPONENT 2 *******************/
  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }

  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  /********** SPIN COMPONENT 2 ************************/
  /* Point outptr to base of current zptr atom for storing */
  queue_iadd_imm(outptr,zptr,ZERO);

  co = 0;
  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  co = 1;
  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  co = 2;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
 
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);

  iterate_stream(PreVecY);
  queue_prefetch(PreVecY);

  iterate_stream(PreVecZ);
  
  /*************** PRELOAD SPIN COMPONENT 3 *****************/
  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }

  co=0;
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  co=1;
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  /************* DO SPIN COMPONENT 3 ************************/
  queue_iadd_imm(outptr,zptr,ZERO);

  co = 0;
  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  co = 1;
  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  co = 2;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
 
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  for( rei=0; rei<2; rei++) {
    if( !do_minus ) { 
      queue_fmadd(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    else {
      queue_fnmsub(z[co][rei], a, y[co][rei], x[co][rei]);
    }
    queue_fstore(z[co][rei], VEC_IMM[co][rei], outptr);
  }

  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);

  iterate_stream(PreVecY);
  iterate_stream(PreVecZ);

  /***************** END LOOP *****************************/
  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}







