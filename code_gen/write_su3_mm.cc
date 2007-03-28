
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


void bagel_su3( char *);

#include "bagel_qdp_options.h"

#ifdef USE_DOUBLE
#warning Using Double precision for Matrix
Datum GaugeType=Double;
#else
#warning Using Single Precision for Matrix
Datum GaugeType=Single;
#endif

char name[80]="";

char procname[80] = "UNSPECIFIED"; /*Default processor type*/

int main ( int argc, char *argv[])
{
  int arg;
  char *c;

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"lrn:P:")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<20) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<20) strcpy(procname,optarg); break;
    default: fprintf(stderr,"Usage: %s -[h] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  bagel_su3(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * bagel_su3( A, B, C ) 
   * Cregs A[nvec][Nrow][Ncol]
   * Cregs B[nvec][Nrow][Ncol]
   * Cregs C[nvec][Nrow][Ncol]
   *
   * asmint nvec
   *
   * Load C[i, 0, :] - 18 regs 
   * Load B[i, 0, :]  - 6 regs 

   * for(i=0;i<nvec-1;i++) {
   *     A[i,0,0]  = B[i, 0, 0] * C[i, 0, 0]
   *     A[i,0,1]  = B[i, 0, 0] * C[i, 0, 1]
   *     A[i,0,2]  = B[i, 0, 0] * C[i, 0, 2]
   *     A[i,0,0] += B[i, 0, 1] * C[i, 1, 0]
   *     A[i,0,1] += B[i, 0, 1] * C[i, 1, 1]
   *     A[i,0,2] += B[i, 0, 1] * C[i, 1, 2]
   *     A[i,0,0] += B[i, 0, 2] * C[i, 2, 0]
   *     A[i,0,1] += B[i, 0, 2] * C[i, 2, 1]
   *     A[i,0,2] += B[i, 0, 2] * C[i, 2, 2]
   *     Steram out A[i,0,:]

   *     Load B[i, 1, :]
   *     A[i,1,0]  = B[i, 1, 0] * C[i, 0, 0]
   *     A[i,1,1]  = B[i, 1, 0] * C[i, 0, 1]
   *     A[i,1,2]  = B[i, 1, 0] * C[i, 0, 2]
   *     A[i,1,0] += B[i, 1, 1] * C[i, 1, 0]
   *     A[i,1,1] += B[i, 1, 1] * C[i, 1, 1]
   *     A[i,1,2] += B[i, 1, 1] * C[i, 1, 2]
   *     A[i,1,0] += B[i, 1, 2] * C[i, 2, 0]
   *     A[i,1,1] += B[i, 1, 2] * C[i, 2, 1]
   *     A[i,1,2] += B[i, 1, 2] * C[i, 2, 2]
   *     Stream out A[i,1,:]

   *     Load B[i, 2, :]
   *     A[i,2,0]  = B[i, 2, 0] * C[i, 0, 0]
   *     A[i,2,1]  = B[i, 2, 0] * C[i, 0, 1]
   *     A[i,2,2]  = B[i, 2, 0] * C[i, 0, 2]
   *     A[i,2,0] += B[i, 2, 1] * C[i, 1, 0]
   *     A[i,2,1] += B[i, 2, 1] * C[i, 1, 1]
   *     A[i,2,2] += B[i, 2, 1] * C[i, 1, 2]
   *     A[i,2,0] += B[i, 2, 2] * C[i, 2, 0]
   *     A[i,2,1] += B[i, 2, 2] * C[i, 2, 1]
   *     A[i,2,2] += B[i, 2, 2] * C[i, 2, 2]
   *     Stream out A[i,2,:]
   *     Load C[i+1, :, :] - 18 regs 
   *     Load B[i+1, 0, :]  - 6 regs 
   *
   *  }
   *
   *  A[nvec-1,0,0]  = B[nvec-1, 0, 0] * C[nvec-1, 0, 0]
   *  A[nvec-1,0,1]  = B[nvec-1, 0, 0] * C[nvec-1, 0, 1]
   *  A[nvec-1,0,2]  = B[nvec-1, 0, 0] * C[nvec-1, 0, 2]
   *  A[nvec-1,0,0] += B[nvec-1, 0, 1] * C[nvec-1, 1, 0]
   *  A[nvec-1,0,1] += B[nvec-1, 0, 1] * C[nvec-1, 1, 1]
   *  A[nvec-1,0,2] += B[nvec-1, 0, 1] * C[nvec-1, 1, 2]
   *  A[nvec-1,0,0] += B[nvec-1, 0, 2] * C[nvec-1, 2, 0]
   *  A[nvec-1,0,1] += B[nvec-1, 0, 2] * C[nvec-1, 2, 1]
   *  A[nvec-1,0,2] += B[nvec-1, 0, 2] * C[nvec-1, 2, 2]
   *  Steram out A[nvec-1,0,:]

   *  Load B[nvec-1, 1, :]
   *  A[nvec-1,1,0]  = B[nvec-1, 1, 0] * C[nvec-1, 0, 0]
   *  A[nvec-1,1,1]  = B[nvec-1, 1, 0] * C[nvec-1, 0, 1]
   *  A[nvec-1,1,2]  = B[nvec-1, 1, 0] * C[nvec-1, 0, 2]
   *  A[nvec-1,1,0] += B[nvec-1, 1, 1] * C[nvec-1, 1, 0]
   *  A[nvec-1,1,1] += B[nvec-1, 1, 1] * C[nvec-1, 1, 1]
   *  A[nvec-1,1,2] += B[nvec-1, 1, 1] * C[nvec-1, 1, 2]
   *  A[nvec-1,1,0] += B[nvec-1, 1, 2] * C[nvec-1, 2, 0]
   *  A[nvec-1,1,1] += B[nvec-1, 1, 2] * C[nvec-1, 2, 1]
   *  A[nvec-1,1,2] += B[nvec-1, 1, 2] * C[nvec-1, 2, 2]
   *  Stream out A[nvec-1,1,:]

   *  Load B[nvec-1, 2, :]
   *  A[nvec-1,2,0]  = B[nvec-1, 2, 0] * C[nvec-1, 0, 0]
   *  A[nvec-1,2,1]  = B[nvec-1, 2, 0] * C[nvec-1, 0, 1]
   *  A[nvec-1,2,2]  = B[nvec-1, 2, 0] * C[nvec-1, 0, 2]
   *  A[nvec-1,2,0] += B[nvec-1, 2, 1] * C[nvec-1, 1, 0]
   *  A[nvec-1,2,1] += B[nvec-1, 2, 1] * C[nvec-1, 1, 1]
   *  A[nvec-1,2,2] += B[nvec-1, 2, 1] * C[nvec-1, 1, 2]
   *  A[nvec-1,2,0] += B[nvec-1, 2, 2] * C[nvec-1, 2, 0]
   *  A[nvec-1,2,1] += B[nvec-1, 2, 2] * C[nvec-1, 2, 1]
   *  A[nvec-1,2,2] += B[nvec-1, 2, 2] * C[nvec-1, 2, 2]
   *  Stream out A[nvec-1,2,:]

   * }
   *
   */

void bagel_su3(char *name)
{
  int dum = defargcount(5);

  /*Integer register usage*/
  alreg(Aptr,Iregs);
  alreg(Bptr,Iregs);
  alreg(Cptr,Iregs);
  alreg(counter,Iregs);

  /*Floating register usage*/
  reg_array_1d(A,Cregs,3);
  reg_array_1d(B,Cregs,3);
  reg_array_2d(C,Cregs,3,3);

  offset_2d(ROW,GaugeType,3,2);
  offset_3d(MATRIX,GaugeType, 3,3,2);

  struct stream *PreA;
  struct stream *PreB;
  struct stream *PreC;

  int brchno,retno; /*Branch target handles*/
  int co,rei;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  getarg(Aptr);
  getarg(Bptr);           /*Get args*/
  getarg(Cptr);           /*Get args*/
  getarg(counter);

  def_off(MATRIX_ATOM, GaugeType, 18);
  def_off(ROW_ATOM, GaugeType, 6);

  PreA = create_stream(ROW_ATOM, Aptr,3*counter,STREAM_OUT,LINEAR);
  PreB= create_stream(ROW_ATOM,  Bptr, 3*counter, STREAM_IN ,LINEAR);
  PreC= create_stream(MATRIX_ATOM, Cptr, counter, STREAM_IN ,LINEAR);

  

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */

  brchno = start_loop(counter);

  // Load C
  for(int row =0; row < 3; row++) { 
    for(int col=0; col < 3; col++) { 
      complex_load(C[row][col], MATRIX[row][col][0], Cptr);
    }
  }
  iterate_stream(PreC);

  // Load B 
  for(int col=0; col < 3; col++) { 
    complex_load(B[col], ROW[col][0], Bptr);
  }
  iterate_stream(PreB);
  /*
   *     A[i,0,0]  = B[i, 0, 0] * C[i, 0, 0]
   *     A[i,0,1]  = B[i, 0, 0] * C[i, 0, 1]
   *     A[i,0,2]  = B[i, 0, 0] * C[i, 0, 2]
   *     A[i,0,0] += B[i, 0, 1] * C[i, 1, 0]
   *     A[i,0,1] += B[i, 0, 1] * C[i, 1, 1]
   *     A[i,0,2] += B[i, 0, 1] * C[i, 1, 2]
   *     A[i,0,0] += B[i, 0, 2] * C[i, 2, 0]
   *     A[i,0,1] += B[i, 0, 2] * C[i, 2, 1]
   *     A[i,0,2] += B[i, 0, 2] * C[i, 2, 2]
   */
  complex_mul(A[0], B[0], C[0][0]);
  complex_mul(A[1], B[0], C[0][1]);
  complex_mul(A[2], B[0], C[0][2]);
  complex_madd(A[0], B[1], C[1][0]);
  complex_madd(A[1], B[1], C[1][1]);
  complex_madd(A[2], B[1], C[1][2]);
  complex_madd(A[0], B[2], C[2][0]);
  complex_madd(A[1], B[2], C[2][1]);
  complex_madd(A[2], B[2], C[2][2]);

  for(int col=0; col < 3; col++) {
    complex_store(A[col], ROW[col][0], Aptr);
  }
  iterate_stream(PreA);

  // Next row of B
  for(int col=0; col < 3; col++) { 
    complex_load(B[col], ROW[col][0], Bptr);
  }
  iterate_stream(PreB);

  /*
   *     A[i,0,0]  = B[i, 0, 0] * C[i, 0, 0]
   *     A[i,0,1]  = B[i, 0, 0] * C[i, 0, 1]
   *     A[i,0,2]  = B[i, 0, 0] * C[i, 0, 2]
   *     A[i,0,0] += B[i, 0, 1] * C[i, 1, 0]
   *     A[i,0,1] += B[i, 0, 1] * C[i, 1, 1]
   *     A[i,0,2] += B[i, 0, 1] * C[i, 1, 2]
   *     A[i,0,0] += B[i, 0, 2] * C[i, 2, 0]
   *     A[i,0,1] += B[i, 0, 2] * C[i, 2, 1]
   *     A[i,0,2] += B[i, 0, 2] * C[i, 2, 2]
   */
  complex_mul(A[0], B[0], C[0][0]);
  complex_mul(A[1], B[0], C[0][1]);
  complex_mul(A[2], B[0], C[0][2]);
  complex_madd(A[0], B[1], C[1][0]);
  complex_madd(A[1], B[1], C[1][1]);
  complex_madd(A[2], B[1], C[1][2]);
  complex_madd(A[0], B[2], C[2][0]);
  complex_madd(A[1], B[2], C[2][1]);
  complex_madd(A[2], B[2], C[2][2]);

  for(int col=0; col < 3; col++) {
    complex_store(A[col], ROW[col][0], Aptr);
  }
  iterate_stream(PreA);

  // Next row of B
  for(int col=0; col < 3; col++) { 
    complex_load(B[col], ROW[col][0], Bptr);
  }
  iterate_stream(PreB);

  /*
   *     A[i,0,0]  = B[i, 0, 0] * C[i, 0, 0]
   *     A[i,0,1]  = B[i, 0, 0] * C[i, 0, 1]
   *     A[i,0,2]  = B[i, 0, 0] * C[i, 0, 2]
q   *     A[i,0,0] += B[i, 0, 1] * C[i, 1, 0]
   *     A[i,0,1] += B[i, 0, 1] * C[i, 1, 1]
   *     A[i,0,2] += B[i, 0, 1] * C[i, 1, 2]
   *     A[i,0,0] += B[i, 0, 2] * C[i, 2, 0]
   *     A[i,0,1] += B[i, 0, 2] * C[i, 2, 1]
   *     A[i,0,2] += B[i, 0, 2] * C[i, 2, 2]
   */
  complex_mul(A[0], B[0], C[0][0]);
  complex_mul(A[1], B[0], C[0][1]);
  complex_mul(A[2], B[0], C[0][2]);
  complex_madd(A[0], B[1], C[1][0]);
  complex_madd(A[1], B[1], C[1][1]);
  complex_madd(A[2], B[1], C[1][2]);
  complex_madd(A[0], B[2], C[2][0]);
  complex_madd(A[1], B[2], C[2][1]);
  complex_madd(A[2], B[2], C[2][2]);

  for(int col=0; col < 3; col++) {
    complex_store(A[col], ROW[col][0], Aptr);
  }
  iterate_stream(PreA);


  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








