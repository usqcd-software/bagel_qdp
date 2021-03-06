
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

#define RE 0
#define IM 1

void qdp_vxpag5iy( char *);

char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/

static int do_minus=0;

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"n:P:m")) != EOF){
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
    case 'm':
      {
	do_minus=1;
      }
      break;
    default:
      {
	fprintf(stderr,"Usage: %s -nroutine_name -P proc [-p]\t",argv[0]); 
	exit (1);
      }
      break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  qdp_vxpag5iy(name);


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
 *        z[i][spin][col][re] = x[i][spin][col][re] 
 *                                         -/+ b * i * y[i][spin][col][im];
 *        z[i][spin][col][im] = x[i][spin][col][im] 
 *                                         +/- b * i * y[i][spin][col][re];
 *     }
 *   }
 *   for(spin=2; spin < 4: spin++) { 
 *     for(col=0; col < Ncol; col++) {
 *        z[i][spin][col][re] = x[i][spin][col][re] 
 *                                         +/- b * y[i][spin][col][im];
 *        z[i][spin][col][im] = x[i][spin][col][im] 
 *                                         -/+ b * y[i][spin][col][re];
 *     }
 *   }
 * }
 */

void qdp_vxpag5iy(char *name)
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
  def_dp_off(ZERO,0);

  /* A vector is of length 6 */
  def_dp_off(VEC_ATOM,6);

  /* This I am guessing is some pattern to describe the vectors
     and that it has to match the 2d register allocation above */
  offset_2d_dp(VEC_IMM,3,2);

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

  if( !do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( !do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);
  

  co = 1;

  if( !do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( !do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);


  co = 2;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
 
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  if( !do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( !do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);

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

  if( !do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( !do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);

  co = 1;

  if( !do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( !do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);

  co = 2;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
 
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  if( !do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( !do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);

  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);

  iterate_stream(PreVecY);
  queue_prefetch(PreVecY);

  iterate_stream(PreVecZ);

  /*************** PRELOAD SPIN COMPONENT 2 *****************/
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

  /************* DO SPIN COMPONENT 2 ************************/
  queue_iadd_imm(outptr,zptr,ZERO);

  co = 0;
  if( do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);

  co = 1;
  if( do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);

  co = 2;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
 
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  if( do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);

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
  if( do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);

  co = 1;

  if( do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);

  co = 2;
  for (rei=0;rei<2;rei++){
    queue_fload(x[co][rei],VEC_IMM[co][rei],xptr);
  }
 
  for (rei=0;rei<2;rei++){
    queue_fload(y[co][rei],VEC_IMM[co][rei],yptr);
  }

  if( do_minus ) { 
    queue_fnmsub(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  else {
    queue_fmadd(z[co][RE], a, y[co][IM], x[co][RE]);
  }
  queue_fstore(z[co][RE], VEC_IMM[co][RE], outptr);
  
  if( do_minus ) { 
    queue_fmadd(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  else {
    queue_fnmsub(z[co][IM], a, y[co][RE], x[co][IM]);
  }
  queue_fstore(z[co][IM], VEC_IMM[co][IM], outptr);


  iterate_stream(PreVecX);
  queue_prefetch(PreVecX);

  iterate_stream(PreVecY);
  queue_prefetch(PreVecY);

  iterate_stream(PreVecZ);



  /***************** END LOOP *****************************/
  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








