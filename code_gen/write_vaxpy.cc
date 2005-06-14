
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


void qcdoc_vaxpy( char *);

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
  while ( ( arg = getopt(argc,argv,"lrn:P:m")) != EOF){
    switch(arg){
    case 'n': if (strlen(optarg)<20) strcpy(name,optarg); break;
    case 'P': if (strlen(optarg)<20) strcpy(procname,optarg); break;
    case 'r': overwrite2 = 1; break;
    case 'l': overwrite1 = 1; break;
    case 'm': do_vaxmy = 1; break;
    default: fprintf(stderr,"Usage: %s -[h] -n routine_name -P proc\t",argv[0]); 
             exit (1); break;
    }
  }

  if ( overwrite1 && overwrite2 ) {
    printf("Cannot overwrite both x and y\n");
    exit(-1);
  }

  /*Control output according to user set up args*/
  set_processor_optarg(procname);
  setup_cmadds();

  /*queue the naive asm code*/
  qcdoc_vaxpy(name);

  /*Filter through an virtual out of order processor */
  schedule_for_proc();
  /*Dump the resulting code*/
  dump_instruction_queue();
  return(0);

}

  /*
   * vaxpy ( a, vec1, vec2 , nvec )
   * fpoint a
   * fpoint vec1[nvec][Ncol][rei] 
   * fpoint vec2[nvec][Ncol][rei] 
   * asmint nvec
   *
   * for(i=0;i<nvec;i++)
   *   vec1[i] = A * vec1[i] + vec2 [i];
   *
   */

void qcdoc_vaxpy( char *name)
{
  int dum = defargcount(5);

  /*Integer register usage*/
  alreg(vec1ptr,Iregs);
  alreg(vec2ptr,Iregs);
  alreg(vec3ptr,Iregs);
  alreg(counter,Iregs);
  alreg(Aptr,Iregs);
  alreg(outptr,Iregs);

  /*Floating register usage*/
  alreg(A,Fregs);
  reg_array_2d(vec1,Fregs,3,2);
  reg_array_2d(vec2,Fregs,3,2);
  reg_array_2d(vec3,Fregs,3,2);

  def_off(ZERO,0);
  def_off(VEC_ATOM,6);
  offset_2d(VEC_IMM,3,2);

  struct stream *PreVec1;
  struct stream *PreVec2;
  struct stream *PreVec3;

  int brchno,retno; /*Branch target handles*/
  int co,rei;

  make_inst(DIRECTIVE,Enter_Routine,name);
  grab_stack(0);
  save_regs();

  if ( (!overwrite1) && (!overwrite2) ) {
    getarg(vec3ptr);
  }
  getarg(Aptr);
  getarg(vec1ptr);           /*Get args*/
  getarg(vec2ptr);           /*Get args*/
  getarg(counter);

  queue_fload(A,ZERO,Aptr); /*Find the complex scale factor*/

  PreVec1= create_stream(VEC_ATOM,vec1ptr ,counter,STREAM_IN ,LINEAR);
  PreVec2= create_stream(VEC_ATOM,vec2ptr ,counter,STREAM_IN ,LINEAR);
  if ( (!overwrite1) && (!overwrite2) ) {
    PreVec3 = create_stream(VEC_ATOM,vec3ptr,counter,STREAM_OUT,LINEAR);
  }

  /*Branch to stack restore if length <1*/
  retno = get_target_label();
  check_iterations(counter,retno); 

  /*
   * Start software pipeline
   */

  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
  }

  brchno = start_loop(counter);

  pragma(DCBT_SPACE,5);


  if ( overwrite1 ) {
    queue_iadd_imm(outptr,vec1ptr,ZERO);
  } else if ( overwrite2 ) {
    queue_iadd_imm(outptr,vec2ptr,ZERO);
  } else {
    queue_iadd_imm(outptr,vec3ptr,ZERO);
  }

  co = 0;
  for(rei=0;rei<2;rei++){
    if ( !do_vaxmy ) { 
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    } else {
      make_inst(FMACPIPE,FMSUB,vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    }
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr);
  }
  co = 1;
  for(rei=0;rei<2;rei++){
    if ( !do_vaxmy ) { 
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    } else {
      make_inst(FMACPIPE,FMSUB,vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    }
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr);
  }
  co = 2;
  for(rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
    if ( !do_vaxmy ) { 
      queue_fmadd(vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    } else {
      make_inst(FMACPIPE,FMSUB,vec3[co][rei],A,vec1[co][rei],vec2[co][rei]);
    }
    queue_fstore(vec3[co][rei],VEC_IMM[co][rei],outptr);
  }


  /*
  if ( (!overwrite1) && (!overwrite2) ) {
    queue_prefetch(PreVec3);
  }
  */

  iterate_stream(PreVec1);
  iterate_stream(PreVec2);

  if ( (!overwrite1) && (!overwrite2) ) {
    iterate_stream(PreVec3);
  }

  queue_prefetch(PreVec1);
  queue_prefetch(PreVec2);

  co = 0;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
  }
  co = 1;
  for (rei=0;rei<2;rei++){
    queue_fload(vec1[co][rei],VEC_IMM[co][rei],vec1ptr);
    queue_fload(vec2[co][rei],VEC_IMM[co][rei],vec2ptr);
  }


  stop_loop(brchno,counter);


  make_inst(DIRECTIVE,Target,retno);

  restore_regs();
  free_stack();
  make_inst(DIRECTIVE,Exit_Routine,name);

  return;
}








