
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

extern struct rotating_reg *CMADregs;

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

bool conjB = false;
bool conjC = false;
bool bothConj = false;

int main ( int argc, char *argv[])
{
  int arg;
  char *c;
  name[0]='\0';

  /*Process command line*/
  while ( ( arg = getopt(argc,argv,"lrhHn:P:")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<20) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<20) strcpy(procname,optarg); break;
    case 'h': conjB = true; break;
    case 'H': conjC = true; break;

    default: fprintf(stderr,"Usage: %s -[h] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }

  if( conjB && conjC) bothConj=true;
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

void bagel_su3(char *name)
{
  int dum = defargcount(5);

  /*Integer register usage*/
  /*For arguments */
  alreg(Aptr,Iregs);
  alreg(Bptr,Iregs);
  alreg(Cptr,Iregs);
  alreg(counter,Iregs);
  alreg(OMIPtr, Iregs);
  alreg(Aprev,Iregs);

  

  /* All the matrices */
  alreg(OMI,Cregs);
  
  reg_array_2d(A, Cregs, 3, 3);
  reg_array_2d(B, Cregs, 3, 3);
  reg_array_2d(C, Cregs, 3, 3);

  /* Space between successive matrices */
  def_off(MATRIX_ATOM, GaugeType, 18);

  /* Byte offset for one_minus_i */
  int ioff = def_offset(0, Byte, "ioff");

  /* Offsets in the register arrays */
  offset_3d(MATRIX, GaugeType, 3, 3, 2);  


  struct stream *PreA; /* Prefetching */
  struct stream *PreB;
  struct stream *PreC;

  int brchno,retno; /*Branch target handles*/
  int co,rei; /* Internal counters */
  int i,j,k;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  getarg(Aptr);
  getarg(Bptr);           /*Get args*/
  getarg(Cptr);           /*Get args*/
  getarg(counter);
  getarg(OMIPtr);

  for(int i=0; i < 18; i++) { 
    need_constant(i*2*SizeofDatum(GaugeType));
  }
  
  PreA = create_stream(MATRIX_ATOM, Aptr,counter,STREAM_OUT,LINEAR);
  PreB = create_stream(MATRIX_ATOM, Bptr, counter, STREAM_IN, LINEAR);
  PreC = create_stream(MATRIX_ATOM, Cptr, counter, STREAM_IN, LINEAR);
 
  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */
  if( bothConj ) {
   complex_load(OMI, ioff, OMIPtr, GaugeType);  
  }

  // Make the first Adude the previous one.
  make_inst (IALUPIPE,IOR,Aprev,Aptr,Aptr);

  // Start loop
  brchno = start_loop(counter);

  // Store Last row of A from previous iteration.
  // In first iteration it writes junk but who cares
  for(j = 0;j<3;j++ ){
    complex_store(A[2][j],MATRIX[2][j][0],Aprev,GaugeType);
  }


  /* Load ALL of A & B */
  for(i=0; i < 3; i++) {
    for(j=0; j < 3; j++) { 
      if( conjB == true ) {
	complex_load(B[j][i], MATRIX[i][j][0], Bptr, GaugeType);
      }
      else { 
	complex_load(B[i][j], MATRIX[i][j][0], Bptr, GaugeType);
      }

      if (bothConj ) {
          if( have_hummer() ) {  
	    make_inst( SIMD2HUMMERPIPE, FMUL2, B[j][i], OMI, B[j][i] );
          }
	  else { 
	    queue_fneg(B[j][i]+1, B[j][i]+1);
          }
      }

      if( conjC == true) {
	complex_load(C[j][i], MATRIX[i][j][0], Cptr, GaugeType);
      }
      else { 
	complex_load(C[i][j], MATRIX[i][j][0], Cptr, GaugeType);
      }

      if (bothConj ) { 
	if( have_hummer() ) { 
	   make_inst( SIMD2HUMMERPIPE, FMUL2, C[j][i], OMI, C[j][i]);
        }
	else { 
	  queue_fneg(C[j][i]+1, C[j][i]+1);
        }
      }
    }
  }



  //First row
  i=0;
  
  if ( conjB && (!conjC) ) {
    k=0;
    complex_three_conjmuls(A[i][0], B[i][0], C[0][0],
			   A[i][1], B[i][0], C[0][1],
			   A[i][2], B[i][0], C[0][2]);
    k=1;
    complex_three_conjmadds(A[i][0], B[i][1], C[1][0],
			    A[i][1], B[i][1], C[1][1],
			    A[i][2], B[i][1], C[1][2]);
    k=2;
    complex_three_conjmadds(A[i][0], B[i][2], C[2][0],
			    A[i][1], B[i][2], C[2][1],
			    A[i][2], B[i][2], C[2][2]);
  }
  else if( (!conjB) && (conjC) ) {
    k=0;
    complex_three_conjmuls(A[i][0], C[0][0], B[i][0],
			   A[i][1], C[0][1], B[i][0],
			   A[i][2], C[0][2], B[i][0]);
    k=1;
    complex_three_conjmadds(A[i][0], C[1][0], B[i][1],
			    A[i][1], C[1][1], B[i][1],
			    A[i][2], C[1][2], B[i][1]);
    k=2;
    complex_three_conjmadds(A[i][0], C[2][0], B[i][2],
			    A[i][1], C[2][1], B[i][2],
			    A[i][2], C[2][2], B[i][2]);
  }
  else { 
    k=0;
    complex_three_cmuls(A[i][0], B[i][0], C[0][0],
			A[i][1], B[i][0], C[0][1],
			A[i][2], B[i][0], C[0][2]);
    k=1;
    complex_three_cmadds(A[i][0], B[i][1], C[1][0],
			 A[i][1], B[i][1], C[1][1],
			 A[i][2], B[i][1], C[1][2]);
    k=2;
    complex_three_cmadds(A[i][0], B[i][2], C[2][0],
			 A[i][1], B[i][2], C[2][1],
			 A[i][2], B[i][2], C[2][2]);

 
    
  }

  
  // Store i-th row
  for(j = 0;j<3;j++ ) {
    complex_store(A[i][j],MATRIX[i][j][0],Aptr,GaugeType);
  }
  
  queue_prefetch(PreB);
  queue_prefetch(PreC);

  i=1;
  if ( conjB  && (!conjC) ) {
    k=0;
    complex_three_conjmuls(A[i][0], B[i][0], C[0][0],
			   A[i][1], B[i][0], C[0][1],
			   A[i][2], B[i][0], C[0][2]);
    k=1;
    complex_three_conjmadds(A[i][0], B[i][1], C[1][0],
			    A[i][1], B[i][1], C[1][1],
			    A[i][2], B[i][1], C[1][2]);
    k=2;
    complex_three_conjmadds(A[i][0], B[i][2], C[2][0],
			    A[i][1], B[i][2], C[2][1],
			    A[i][2], B[i][2], C[2][2]);
  }
  else if( (!conjB) && (conjC) ) {
    k=0;
    complex_three_conjmuls(A[i][0], C[0][0], B[i][0],
			   A[i][1], C[0][1], B[i][0],
			   A[i][2], C[0][2], B[i][0]);
    k=1;
    complex_three_conjmadds(A[i][0], C[1][0], B[i][1],
			    A[i][1], C[1][1], B[i][1],
			    A[i][2], C[1][2], B[i][1]);
    k=2;
    complex_three_conjmadds(A[i][0], C[2][0], B[i][2],
			    A[i][1], C[2][1], B[i][2],
			    A[i][2], C[2][2], B[i][2]);
  }
  else { 
    k=0;
    complex_three_cmuls(A[i][0], B[i][0], C[0][0],
			A[i][1], B[i][0], C[0][1],
			A[i][2], B[i][0], C[0][2]);
    k=1;
    complex_three_cmadds(A[i][0], B[i][1], C[1][0],
			 A[i][1], B[i][1], C[1][1],
			 A[i][2], B[i][1], C[1][2]);
    k=2;
    complex_three_cmadds(A[i][0], B[i][2], C[2][0],
			 A[i][1], B[i][2], C[2][1],
			 A[i][2], B[i][2], C[2][2]);

  
  }

  
  // Store i-th row
  for(j = 0;j<3;j++ ){
    complex_store(A[i][j],MATRIX[i][j][0],Aptr,GaugeType);
  }


  i=2;
  if ( conjB  && (!conjC) ) {
    k=0;
    complex_three_conjmuls(A[i][0], B[i][0], C[0][0],
			   A[i][1], B[i][0], C[0][1],
			   A[i][2], B[i][0], C[0][2]);
    k=1;
    complex_three_conjmadds(A[i][0], B[i][1], C[1][0],
			    A[i][1], B[i][1], C[1][1],
			    A[i][2], B[i][1], C[1][2]);
    k=2;
    complex_three_conjmadds(A[i][0], B[i][2], C[2][0],
			    A[i][1], B[i][2], C[2][1],
			    A[i][2], B[i][2], C[2][2]);
  }
  else if( (!conjB) && (conjC) ) {
    k=0;
    complex_three_conjmuls(A[i][0], C[0][0], B[i][0],
			   A[i][1], C[0][1], B[i][0],
			   A[i][2], C[0][2], B[i][0]);
    k=1;
    complex_three_conjmadds(A[i][0], C[1][0], B[i][1],
			    A[i][1], C[1][1], B[i][1],
			    A[i][2], C[1][2], B[i][1]);
    k=2;
    complex_three_conjmadds(A[i][0], C[2][0], B[i][2],
			    A[i][1], C[2][1], B[i][2],
			    A[i][2], C[2][2], B[i][2]);
  }
  else { 
    // This is what to do if either none, or both are conjugated
    // in the case where both are conjugated we need different loading 
    // and storing
    k=0;
    complex_three_cmuls(A[i][0], B[i][0], C[0][0],
			A[i][1], B[i][0], C[0][1],
			A[i][2], B[i][0], C[0][2]);
    k=1;
    complex_three_cmadds(A[i][0], B[i][1], C[1][0],
			 A[i][1], B[i][1], C[1][1],
			 A[i][2], B[i][1], C[1][2]);
    k=2;
    complex_three_cmadds(A[i][0], B[i][2], C[2][0],
			 A[i][1], B[i][2], C[2][1],
			 A[i][2], B[i][2], C[2][2]);


  }

  
  make_inst(IALUPIPE,IOR,Aprev,Aptr,Aptr);
 

  iterate_stream(PreA);
  iterate_stream(PreB);
  iterate_stream(PreC);


  stop_loop(brchno,counter);

  // The last row of the last one, falls through
  for(j = 0;j<3;j++ ){
    complex_store(A[2][j],MATRIX[2][j][0],Aprev,GaugeType);
  }

  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}
 







