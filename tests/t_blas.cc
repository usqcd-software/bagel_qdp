// $Id: t_blas.cc,v 1.1 2005-06-14 17:13:32 edwards Exp $

#include <iostream>
#include <iomanip>
#include <cstdio>

#include <time.h>

#include "qdp.h"
#include "bagel_qdp.h"

#include "scalarsite_generic/generic_blas_vadd.h"
#include "scalarsite_generic/generic_blas_vscal.h"
#include "scalarsite_generic/generic_blas_vscal_g5.h"
using namespace std;
using namespace QDP;

// time axpy like ops.
void time_five_op(void f1(REAL*, REAL*, REAL*, REAL*, int),
		  std::string name_f1,
		  void f2(REAL*, REAL*, REAL*, REAL*, int),
		  std::string name_f2,
		  REAL *arg1,
		  REAL *arg2,
		  REAL *arg3,
		  REAL *arg4,
		  int  arg5,
		  int site_flops);

void time_six_op(void f1(REAL*, REAL*, REAL*, REAL*, REAL*, int),
		 std::string name_f1,
		 void f2(REAL*, REAL*, REAL*, REAL*, REAL*, int),
		 std::string name_f2,
		 REAL *arg1,
		 REAL *arg2,
		 REAL *arg3,
		 REAL *arg4,
		 REAL *arg5,
		 int  arg6,
		 int site_flops);

void time_four_op(void f1(REAL*, REAL*, REAL*, int),
		  std::string name_f1,
		  void f2(REAL*, REAL*, REAL*, int),
		  std::string name_f2,
		  REAL *arg1,
		  REAL *arg2,
		  REAL *arg3,
		  int  arg4,
		  int site_flops);

// Time 3 op gsum thingies
void time_norm2(LatticeFermion& x);

void time_innerProd(LatticeFermion& x, LatticeFermion& y);

void time_innerProdReal(LatticeFermion& x, LatticeFermion& y);

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);
  
  // Setup the layout
  // const int foo[] = {8,8,16,16};
  const int local_vol[] = { 4,4,4,4 };
  multi1d<int> nrow(Nd);
  nrow=local_vol;

  Layout::setLattSize(nrow);
  Layout::create();
  for(int i=0; i < Nd; i++) {
    nrow[i] = local_vol[i]*Layout::logicalSize()[i];
  }
  Layout::setLattSize(nrow);
  Layout::create();

  Real a=Real(1.5);
  Real b=Real(-4.2);
  LatticeFermion qx;
  LatticeFermion qy;
  LatticeFermion qz;
  LatticeFermion qz2;
  LatticeFermion qtmp;
  Double dnorm;
  Double dnorm2;
  LatticeFermion diff;
  REAL *aptr;
  REAL *bptr;
  REAL *x;
  REAL *y;
  REAL *z;
  const int Ncmpx=2;
  int n_3vec;
  int n_4vec;


  // Test y += a*x
  gaussian(qx);
  gaussian(qy);

  qtmp = a*qx;
  qz = qy + qtmp;  
  aptr = (REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_3vec = (all.end() - all.start() + 1)*Ns;
  qdp_vaxpy3(z, aptr, x, y, n_3vec);
  diff = qz - qy;
  dnorm =norm2(diff);
  QDPIO::cout << "axpy3 diff = " << dnorm  << endl;
  time_five_op( vaxpy3, "generic vaxpy", 
		qdp_vaxpy3, "qdp vaxpy", 
		z, aptr, x, y, n_3vec,
		2*Ncmpx*Nc*Ns);
  
  gaussian(qx);
  gaussian(qy);
  qtmp = a*qx;
  qz = qtmp - qy;  
  aptr = (REAL *)&(a.elem().elem().elem().elem());
  x  = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_3vec = (all.end() - all.start() + 1)*Ns;
  qdp_vaxmy3(z, aptr, x, y, n_3vec);
  diff = qz - qy;
  dnorm=norm2(diff);
  QDPIO::cout << "axmy3 diff = " << dnorm << endl;

  gaussian(qx);
  gaussian(qy);

  qz = a*qx + b*qy;
  aptr = (REAL *)&(a.elem().elem().elem().elem());
  bptr = (REAL *)&(b.elem().elem().elem().elem());
  x  = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_3vec = (all.end() - all.start() + 1)*Ns;
  qdp_vaxpby3(z,aptr,x,bptr,y,n_3vec);
  diff = qz2 - qz;
  dnorm=norm2(diff);
  QDPIO::cout << "axpby3 diff = " << dnorm << endl;
  time_six_op( vaxpby3, "generic_vaxpby", 
	       qdp_vaxpby3, "qdp_vaxpby", 
	       z, aptr, x, bptr, y, n_3vec,
	       3*Ncmpx*Nc*Ns);
  
  gaussian(qx);
  gaussian(qy);

  qz = a*qx - b*qy;
  aptr = (REAL *)&(a.elem().elem().elem().elem());
  bptr = (REAL *)&(b.elem().elem().elem().elem());
  x  = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_3vec = (all.end() - all.start() + 1)*Ns;
  qdp_vaxmby3(z,aptr,x,bptr,y,n_3vec);
  diff = qz2 - qz;
  dnorm=norm2(diff);
  QDPIO::cout << "axmby3 diff = " << dnorm << endl;


  gaussian(qx);
  gaussian(qy);
  qz = qx + qy;
  x  = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_3vec = (all.end() - all.start() + 1)*Ns;
  qdp_vadd3(z,x,y,n_3vec);
  diff = qz2 - qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vadd diff = " << dnorm << endl;
  time_four_op( vadd, "generic_vadd", 
		qdp_vadd3, "qdp_vadd", 
	       z, x, y, n_3vec,
	       1*Ncmpx*Nc*Ns);
  
  gaussian(qx);
  gaussian(qy);
  qz = qx - qy;
  x  = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_3vec = (all.end() - all.start() + 1)*Ns;
  qdp_vsub3(z,x,y,n_3vec);
  diff = qz2 - qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vsub diff = " << dnorm << endl;

  gaussian(qx);
  qz = a*qx;
  aptr = (REAL *)&(a.elem().elem().elem().elem());  
  x  = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_3vec = (all.end() - all.start() + 1)*Ns;
  qdp_vscal3(z,aptr,x,n_3vec);
  Internal::globalSum(*aptr);
  diff = qz2 - qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vscal diff = " << dnorm << endl;
  time_four_op( vscal, "generic_vscal", 
		qdp_vscal3, "qdp_vscal", 
	       z, aptr, y, n_3vec,
	       1*Ncmpx*Nc*Ns);

  gaussian(qx);
  Real n2=norm2(qx);
  
  x  =( REAL *)&(qx.elem(0).elem(0).elem(0).real());
  bptr = (REAL *)&(b.elem().elem().elem().elem());  
  n_3vec = (all.end() - all.start() + 1)*Ns;
  qdp_lsum2(bptr, x, n_3vec);
  Internal::globalSum(*bptr);
  QDPIO::cout << " norm2() = " << n2 << " b=" << b << endl;
  Real diff_r = n2-b;
  QDPIO::cout << "norm diff = " << diff_r << endl;
  time_norm2(qx);

  gaussian(qx);
  gaussian(qy);

  Complex cinner= innerProduct(qx,qy);
  Complex c2;
  aptr = (REAL*)&(c2.elem().elem().elem().real());
  bptr = (REAL*)&(c2.elem().elem().elem().imag());
  x  =( REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y  =( REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_3vec = (all.end() - all.start() + 1)*Ns;
  qdp_lcdot(aptr,bptr,x,y, n_3vec);
  Internal::globalSumArray(aptr,2);
  Complex cdiff = c2 -cinner;
  QDPIO::cout << " innerProd diff = " << cdiff << endl;
  time_innerProd(qx,qy);

  gaussian(qx);
  gaussian(qy);
  Double cinnerr = innerProductReal(qx,qy);
  Real c2r;
  
  aptr = (REAL*)&(c2r.elem().elem().elem().elem());
  x  =( REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y  =( REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_3vec = ((all.end() - all.start() + 1)*Ns);
  *aptr=0;
  qdp_lcdotr(aptr,x,y, n_3vec);
  Internal::globalSum(*aptr);
  Real rdiff=cinnerr -c2r;
  QDPIO::cout << "innerProductReal diff = " << rdiff << endl;
  time_innerProdReal(qx,qy);




  gaussian(qx);
  a=Real(2.5);
  qz = a*chiralProjectPlus(qx);
  aptr = (REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vscal_chp(y,aptr,x, n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vscal ProjPlus diff = " << dnorm << endl;

  gaussian(qx);
  a=Real(2.5);
  qz = a*chiralProjectMinus(qx);
  aptr = (REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vscal_chm(y,aptr,x, n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vscal ProjMinus diff = " << dnorm << endl;


  time_four_op(scal_g5ProjPlus, "scal_g5ProjPlus",
	       qdp_vscal_chp, "qdp_vscal_chp",
	       y,aptr,x,n_4vec,
	       4*Nc);
  time_four_op(scal_g5ProjMinus, "scal_g5ProjMinus",
	       qdp_vscal_chp, "qdp_vscal_chm",
	       y,aptr,x,n_4vec,
	       4*Nc);

  gaussian(qx);
  gaussian(qy);
  
  qz = qx + chiralProjectPlus(qy);
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vadd_chp(z,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vadd ProjPlus diff = " << dnorm << endl;

  qz = qx + chiralProjectMinus(qy);
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vadd_chm(z,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vadd ProjMinus diff = " << dnorm << endl;

  gaussian(qx);
  gaussian(qy);
  
  qz = qx - chiralProjectPlus(qy);
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vsub_chp(z,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vsub ProjPlus diff = " << dnorm << endl;

  qz = qx - chiralProjectMinus(qy);
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vsub_chm(z,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vsub ProjMinus diff = " << dnorm << endl;


  time_four_op(add_g5ProjPlus, "add_g5ProjPlus",
	       qdp_vadd_chp, "qdp_vadd_chp",
	       z,x,y,n_4vec,
	       4*Nc);

  time_four_op(add_g5ProjMinus, "add_g5ProjMinus",
	       qdp_vadd_chm, "qdp_vadd_chm",
	       z,x,y,n_4vec,
	       4*Nc);

  time_four_op(sub_g5ProjPlus, "sub_g5ProjPlus",
	       qdp_vsub_chp, "qdp_vsub_chp",
	       z,x,y,n_4vec,
	       4*Nc);

  time_four_op(sub_g5ProjMinus, "sub_g5ProjMinus",
	       qdp_vsub_chm, "qdp_vsub_chm",
	       z,x,y,n_4vec,
	       4*Nc);


  a=Real(2.4);
  qz = a*qx + chiralProjectPlus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vaxpy_chp(z,aptr,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vaxpy ProjPlus diff = " << dnorm << endl;

  qz = a*qx + chiralProjectMinus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vaxpy_chm(z,aptr, x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vaxpy ProjMinus diff = " << dnorm << endl;

  gaussian(qx);
  gaussian(qy);
  
  qz = a*qx - chiralProjectPlus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vaxmy_chp(z,aptr,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vaxmy ProjPlus diff = " << dnorm << endl;

  qz = a*qx - chiralProjectMinus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vaxmy_chm(z,aptr,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vaxmy ProjMinus diff = " << dnorm << endl;


  time_five_op(axpyz_g5ProjPlus, "axpyz_g5ProjPlus",
	       qdp_vaxpy_chp, "qdp_vaxpy_chp",
	       z,aptr,x,y,n_4vec,
	       6*Nc);

  time_five_op(axpyz_g5ProjMinus, "axpyz_g5ProjMinus",
	       qdp_vaxpy_chm, "qdp_vaxpy_chm",
	       z,aptr,x,y,n_4vec,
	       6*Nc);

  time_five_op(axmyz_g5ProjPlus, "axmyz_g5ProjPlus",
	       qdp_vaxmy_chp, "qdp_vaxmy_chp",
	       z,aptr,x,y,n_4vec,
	       6*Nc+Nc);

  time_five_op(axmyz_g5ProjMinus, "axmyz_g5ProjMinus",
	       qdp_vaxmy_chm, "qdp_vaxmy_chm",
	       z,aptr,x,y,n_4vec,
	       6*Nc);


  gaussian(qx);
  gaussian(qy);

  a=Real(-1.5);
  qz = qx + a*chiralProjectPlus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());

  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vxpay_chp(z,aptr,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vxpay ProjPlus diff = " << dnorm << endl;

  qz = qx + a* chiralProjectMinus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vxpay_chm(z,aptr, x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vxpay ProjMinus diff = " << dnorm << endl;

  gaussian(qx);
  gaussian(qy);
  
  qz = qx - a*chiralProjectPlus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vxmay_chp(z,aptr,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vxmay ProjPlus diff = " << dnorm << endl;

  qz = qx - a*chiralProjectMinus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vxmay_chm(z,aptr,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vxmay ProjMinus diff = " << dnorm << endl;


  time_five_op(xpayz_g5ProjPlus, "xpayz_g5ProjPlus",
	       qdp_vxpay_chp, "qdp_vxpay_chp",
	       z,aptr,x,y,n_4vec,
	       4*Nc);

  time_five_op(xpayz_g5ProjMinus, "xpayz_g5ProjMinus",
	       qdp_vxpay_chm, "qdp_vxpay_chm",
	       z,aptr,x,y,n_4vec,
	       4*Nc);

  time_five_op(xmayz_g5ProjPlus, "xmyaz_g5ProjPlus",
	       qdp_vxmay_chp, "qdp_vxmay_chp",
	       z,aptr,x,y,n_4vec,
	       4*Nc);

  time_five_op(xmayz_g5ProjMinus, "xmayz_g5ProjMinus",
	       qdp_vxmay_chm, "qdp_vxmay_chm",
	       z,aptr,x,y,n_4vec,
	       4*Nc);

  gaussian(qx);
  gaussian(qy);

  a=Real(-1.5);
  b=Real(.7894);

  qz = a*qx + b*chiralProjectPlus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  bptr=(REAL *)&(b.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vaxpby_chp(z,aptr,x,bptr,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vaxpby ProjPlus diff = " << dnorm << endl;

  qz = a*qx + b*chiralProjectMinus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  bptr=(REAL *)&(b.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vaxpby_chm(z,aptr, x, bptr, y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vaxpby ProjMinus diff = " << dnorm << endl;

  gaussian(qx);
  gaussian(qy);
  
  qz = a*qx - b*chiralProjectPlus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  bptr=(REAL *)&(b.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vaxmby_chp(z,aptr,x,bptr,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vaxmby ProjPlus diff = " << dnorm << endl;

  qz = a*qx - b*chiralProjectMinus(qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  bptr=(REAL *)&(b.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vaxmby_chm(z,aptr,x,bptr,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vaxmby ProjMinus diff = " << dnorm << endl;


  time_six_op(axpbyz_g5ProjPlus, "axpbyz_g5ProjPlus",
	       qdp_vaxpby_chp, "qdp_vaxpby_chp",
	       z,aptr,x,bptr,y,n_4vec,
	       14*Nc);

  time_six_op(axpbyz_g5ProjMinus, "axpbyz_g5ProjMinus",
	       qdp_vaxpby_chm, "qdp_vaxpby_chm",
	       z,aptr,x,bptr,y,n_4vec,
	       14*Nc);

  time_six_op(axmbyz_g5ProjPlus, "axmbyz_g5ProjPlus",
	       qdp_vaxmby_chp, "qdp_vaxmby_chp",
	       z,aptr,x,bptr,y,n_4vec,
	       14*Nc);

  time_six_op(axmbyz_g5ProjMinus, "axmbyz_g5ProjMinus",
	       qdp_vaxmby_chm, "qdp_vaxmby_chm",
	       z,aptr,x,bptr,y,n_4vec,
	       14*Nc);



  // Gamma 5 variants
  gaussian(qx);
  a=Real(0.345);
  qz = Gamma(15)*qx;
  qz*=a;
  aptr = (REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec=(all.end() - all.start() + 1);
  qdp_vscal_g5(y,aptr,x,n_4vec);
  diff = qz2-qz;
  dnorm = norm2(diff);
  QDPIO::cout << "vscal_g5 diff = " << dnorm << endl;
  time_four_op(scal_g5, "scal_g5",
	       qdp_vscal_g5, "qdp_vscal_g5",
	       y,aptr,x,n_4vec,
	       2*Nc*Ns);

  qz = a*qx + b*(GammaConst<Ns,Ns*Ns-1>()*qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  bptr=(REAL *)&(b.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vaxpbg5y(z,aptr,x,bptr,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vaxpbg5y diff = " << dnorm << endl;

  qz = a*qx - b*(GammaConst<Ns,Ns*Ns-1>()*qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  bptr=(REAL *)&(b.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vaxmbg5y(z,aptr,x,bptr,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vaxmbg5y diff = " << dnorm << endl;

  time_six_op(axpbyz_g5, "axpbyz_g5",
	       qdp_vaxpbg5y, "qdp_vaxpbg5y",
	       z,aptr,x,bptr,y,n_4vec,
	       3*2*Nc*Ns);

  QDPIO::cout << "No generic for vaxmbg5y for timing" << endl;

  qz = qx + a*(GammaConst<Ns,Ns*Ns-1>()*qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vxpag5y(z,aptr,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vxpag5y diff = " << dnorm << endl;

  qz = qx - a*(GammaConst<Ns,Ns*Ns-1>()*qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  bptr=(REAL *)&(b.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vxmag5y(z,aptr,x,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vxmag5y diff = " << dnorm << endl;
 
  QDPIO::cout << "No generic xpayz_g5 for timing " << endl;
  
  time_five_op(xmayz_g5, "xpayz_g5",
	       qdp_vxpag5y, "qdp_vxpag5y",
	       z,aptr,x,y,n_4vec,
	       4*Nc*Ns);

  qz = a*(GammaConst<Ns,Ns*Ns-1>()*qx) + b*qy;
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  bptr=(REAL *)&(b.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vag5xpby(z,aptr,x,bptr,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vag5xpby diff = " << dnorm << endl;

  qz = a*(GammaConst<Ns,Ns*Ns-1>()*qx) - b*qy;
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  bptr=(REAL *)&(b.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vag5xmby(z,aptr,x,bptr,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vag5xmby diff = " << dnorm << endl;

  QDPIO::cout << "No generic vag5xmby exists for timing. Using generic axpbyz_g5 instead " << endl;
  time_six_op(axpbyz_g5, "axpbyz_g5",
	       qdp_vag5xmby, "qdp_vag5xmby",
	       z,aptr,x,bptr,y,n_4vec,
	       3*2*Nc*Ns);

  QDPIO::cout << "No generic for vaxmbg5y for timing" << endl;

  qz = GammaConst<Ns,Ns*Ns-1>()*(a*qx + b*qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  bptr=(REAL *)&(b.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vg5axpby(z,aptr,x,bptr,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vg5axpby diff = " << dnorm << endl;

  qz = GammaConst<Ns,Ns*Ns-1>()*(a*qx - b*qy);
  aptr=(REAL *)&(a.elem().elem().elem().elem());
  bptr=(REAL *)&(b.elem().elem().elem().elem());
  x = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  y = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  z = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  n_4vec = (all.end() - all.start() + 1);
  qdp_vg5axmby(z,aptr,x,bptr,y,n_4vec);
  diff=qz2-qz;
  dnorm=norm2(diff);
  QDPIO::cout << "vg5axmby diff = " << dnorm << endl;

  QDPIO::cout << "No generic for vg5axpby for timing" << endl;

  time_six_op(g5_axmbyz, "g5_axmbyz",
	       qdp_vg5axmby, "qdp_vg5axmby",
	       z,aptr,x,bptr,y,n_4vec,
	       3*2*Nc*Ns);


 
  // Time to bolt
  QDP_finalize();

  exit(0);
}

// Time five operand functions:
//  eg: axpy3: z=ax + y (z, a, x, y, count)
//      innerProduct(x,y): (out_re, out_im, x, y, count)
void time_five_op(void f1(REAL*, REAL*, REAL*, REAL*, int),
		  std::string name_f1,
		  void f2(REAL*, REAL*, REAL*, REAL*, int),
		  std::string name_f2,
		  REAL *arg1,
		  REAL *arg2,
		  REAL *arg3,
		  REAL *arg4,
		  int  arg5,
		  int site_flops)
{
  double flopcount1;
  double flopcount2;

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    while( time < 1.0 ) { 
      iter *=2;
      //      QDPIO::cout << "Calling " << name_f1 << " with " << iter << " iterations " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	f1(arg1, arg2, arg3, arg4, arg5);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing " << name_f1 << " with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
      f1(arg1,arg2,arg3,arg4,arg5);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();
    
    double Nflops = (double)(site_flops*Layout::sitesOnNode()*iter);
    flopcount1 = Nflops/time;
    QDPIO::cout << "Function " << name_f1 << " time: " << time << " (us), " << flopcount1 << "Mflops per node " << endl;
  }

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    while( time < 1.0 ) { 
      iter *=2;
      // QDPIO::cout << "Calling " << name_f2 << " with " << iter << " iterations " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	f2(arg1, arg2, arg3, arg4, arg5);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing " << name_f2 << " with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
      f2(arg1,arg2,arg3,arg4,arg5);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();
    
    double Nflops = (double)(site_flops*Layout::sitesOnNode()*iter);
    flopcount2 = Nflops/time;
    QDPIO::cout << "Function " << name_f2 << " time: " << time << " (us), " << flopcount2 << "Mflops per node " << endl;
  }

  QDPIO::cout << "Fn: " << name_f1 << " flopcount: "  << flopcount1 << " Mflop/s" << endl;
  QDPIO::cout << "Fn: " << name_f2 << " flopcount: "  << flopcount2 << " Mflop/s" << endl;
}

// Time 6 operand functions:
// eg: z=a*x + b*y (z, a, x, b, y, count)
void time_six_op(void f1(REAL*, REAL*, REAL*, REAL*, REAL*, int),
		std::string name_f1,
		void f2(REAL*, REAL*, REAL*, REAL*, REAL*, int),
		std::string name_f2,
		REAL *arg1,
		REAL *arg2,
		REAL *arg3,
		REAL *arg4,
		REAL *arg5,
		int  arg6,
		int site_flops)
{
  double flopcount1;
  double flopcount2;

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    while( time < 1.0 ) { 
      iter *=2;
      // QDPIO::cout << "Calling " << name_f1 << " with " << iter << " iterations " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	f1(arg1, arg2, arg3, arg4, arg5,arg6);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing " << name_f1 << " with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
      f1(arg1,arg2,arg3,arg4,arg5,arg6);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();
    
    double Nflops = (double)(site_flops*Layout::sitesOnNode()*iter);
    flopcount1 = Nflops/time;
    QDPIO::cout << "Function " << name_f1 << " time: " << time << " (us), " << flopcount1 << "Mflops per node " << endl;
  }

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    while( time < 1.0 ) { 
      iter *=2;
      //      QDPIO::cout << "Calling " << name_f2 << " with " << iter << " iterations " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	f2(arg1, arg2, arg3, arg4, arg5, arg6);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing " << name_f2 << " with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
      f2(arg1,arg2,arg3,arg4,arg5,arg6);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();
    
    double Nflops = (double)(site_flops*Layout::sitesOnNode()*iter);
    flopcount2 = Nflops/time;
    QDPIO::cout << "Function " << name_f2 << " time: " << time << " (us), " << flopcount2 << "Mflops per node " << endl;
  }

  QDPIO::cout << "Fn: " << name_f1 << " flopcount: "  << flopcount1 << " Mflop/s" << endl;
  QDPIO::cout << "Fn: " << name_f2 << " flopcount: "  << flopcount2 << " Mflop/s" << endl;
}

// Time four operand thingies:
// eg: z = x +/- y   (z, x, y, count)
//     z = a*x       (z, a, x, count)
void time_four_op(void f1(REAL*, REAL*, REAL*, int),
		  std::string name_f1,
		  void f2(REAL*, REAL*, REAL*, int),
		  std::string name_f2,
		  REAL *arg1,
		  REAL *arg2,
		  REAL *arg3,
		  int  arg4,
		  int site_flops)
{
  double flopcount1;
  double flopcount2;

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    while( time < 1.0 ) { 
      iter *=2;
      // QDPIO::cout << "Calling " << name_f1 << " with " << iter << " iterations " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	f1(arg1, arg2, arg3, arg4);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing " << name_f1 << " with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
      f1(arg1,arg2,arg3,arg4);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();
    
    double Nflops = (double)(site_flops*Layout::sitesOnNode()*iter);
    flopcount1 = Nflops/time;
    QDPIO::cout << "Function " << name_f1 << " time: " << time << " (us), " << flopcount1 << "Mflops per node " << endl;
  }

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    while( time < 1.0 ) { 
      iter *=2;
      //      QDPIO::cout << "Calling " << name_f2 << " with " << iter << " iterations " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	f2(arg1, arg2, arg3, arg4);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing " << name_f2 << " with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
      f2(arg1,arg2,arg3,arg4);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();
    
    double Nflops = (double)(site_flops*Layout::sitesOnNode()*iter);
    flopcount2 = Nflops/time;
    QDPIO::cout << "Function " << name_f2 << " time: " << time << " (us), " << flopcount2 << "Mflops per node " << endl;
  }

  QDPIO::cout << "Fn: " << name_f1 << " flopcount: "  << flopcount1 << " Mflop/s" << endl;
  QDPIO::cout << "Fn: " << name_f2 << " flopcount: "  << flopcount2 << " Mflop/s" << endl;
}

// Time 3 op gsum thingies
void time_norm2(LatticeFermion& x)
{
  double flopcount1;
  double flopcount2;

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    Double foo;
    while( time < 1.0 ) { 
      iter *=2;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	foo=norm2(x);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing norm2 with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
      foo=norm2(x);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();

    // 3NcNs * sites for the norm2
    // + Nsites-1 additions for the sum
    // (3NcNs + 1)*sites-1
    double Nflops = (double)( ((3*Nc*Ns+1)*Layout::sitesOnNode()-1)*iter);
    flopcount1 = Nflops/time;
    QDPIO::cout << "norm2 time: " << time << " (us), " << flopcount1 << "Mflops per node " << endl;
  }

  Real the_sum;
  REAL* sumptr=(REAL *)&(the_sum.elem().elem().elem().elem());
  REAL* vptr  =(REAL *)&(x.elem(0).elem(0).elem(0).real());
  int n_count = (all.end()-all.start()+1)*Ns;

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    while( time < 1.0 ) { 
      iter *=2;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qdp_lsum2(sumptr, vptr, n_count);
	Internal::globalSum(*sumptr);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing qdp_lsum2 with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
	qdp_lsum2(sumptr, vptr, n_count);
	Internal::globalSum(*sumptr);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();

    double Nflops = (double)( ((3*Nc*Ns+1)*Layout::sitesOnNode()-1)*iter);    
    flopcount2 = Nflops/time;
    QDPIO::cout << "qdp_lsum time: " << time << " (us), " << flopcount2 << "Mflops per node " << endl;
  }

  QDPIO::cout << "Fn: norm2() flopcount: "  << flopcount1 << " Mflop/s" << endl;
  QDPIO::cout << "Fn: qdp_lsum2()+gum flopcount: "  << flopcount2 << " Mflop/s" << endl;
}


void time_innerProd(LatticeFermion& x, LatticeFermion& y)
{
  double flopcount1;
  double flopcount2;

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    Complex foo;
    while( time < 1.0 ) { 
      iter *=2;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	foo=innerProduct(x,y);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing innerProd with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
      foo=innerProduct(x,y);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();

    // 3NcNs * sites for the norm2
    // + Nsites-1 additions for the sum
    // (3NcNs + 1)*sites-1
    double Nflops = (double)(2*((3*Nc*Ns+1)*Layout::sitesOnNode()-1)*iter);
    flopcount1 = Nflops/time;
    QDPIO::cout << "innerProd time: " << time << " (us), " << flopcount1 << "Mflops per node " << endl;
  }

  Complex the_sum;
  REAL* sumrptr=(REAL *)&(the_sum.elem().elem().elem().real());
  REAL* sumiptr=(REAL *)&(the_sum.elem().elem().elem().imag());
  REAL* v1ptr  =(REAL *)&(x.elem(0).elem(0).elem(0).real());
  REAL* v2ptr  =(REAL *)&(y.elem(0).elem(0).elem(0).real());
  int n_count = (all.end()-all.start()+1)*Ns;

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    while( time < 1.0 ) { 
      iter *=2;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qdp_lcdot(sumrptr, sumiptr, v1ptr, v2ptr, n_count);
	Internal::globalSumArray(sumrptr,2);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing qdp_lcdot+gsum with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
	qdp_lcdot(sumrptr, sumiptr, v1ptr, v2ptr, n_count);
	Internal::globalSumArray(sumrptr,2);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();

    double Nflops = (double)( 2*((3*Nc*Ns+1)*Layout::sitesOnNode()-1)*iter);    
    flopcount2 = Nflops/time;
    QDPIO::cout << "qdp_lcdot+gsum time: " << time << " (us), " << flopcount2 << "Mflops per node " << endl;
  }

  QDPIO::cout << "Fn: innerProd() flopcount: "  << flopcount1 << " Mflop/s" << endl;
  QDPIO::cout << "Fn: qdp_lcdot()+gum flopcount: "  << flopcount2 << " Mflop/s" << endl;
}


void time_innerProdReal(LatticeFermion& x, LatticeFermion& y)
{
  double flopcount1;
  double flopcount2;

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    Real foo;
    while( time < 1.0 ) { 
      iter *=2;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	foo=innerProductReal(x,y);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing innerProdReal with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
      foo=innerProductReal(x,y);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();

    // 3NcNs * sites for the norm2
    // + Nsites-1 additions for the sum
    // (3NcNs + 1)*sites-1
    double Nflops = (double)(((3*Nc*Ns+1)*Layout::sitesOnNode()-1)*iter);
    flopcount1 = Nflops/time;
    QDPIO::cout << "innerProdReal time: " << time << " (us), " << flopcount1 << "Mflops per node " << endl;
  }

  Real the_sum;
  REAL* sumptr=(REAL *)&(the_sum.elem().elem().elem().elem());
  REAL* v1ptr  =(REAL *)&(x.elem(0).elem(0).elem(0).real());
  REAL* v2ptr  =(REAL *)&(y.elem(0).elem(0).elem(0).real());
  int n_count = (all.end()-all.start()+1)*Ns;

  {
    StopWatch swatch;
    double time=0;
    int iter =1;
    while( time < 1.0 ) { 
      iter *=2;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qdp_lcdotr(sumptr, v1ptr, v2ptr, n_count);
	Internal::globalSum(the_sum);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      Internal::globalSum(time);
      time /= (double)Layout::numNodes();
    }
    QDPIO::cout << "Timing qdp_lcdotr+gsum with " << iter << " iterations "<< endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) { 
	qdp_lcdotr(sumptr, v1ptr, v2ptr, n_count);
	Internal::globalSum(the_sum);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();

    double Nflops = (double)(((3*Nc*Ns+1)*Layout::sitesOnNode()-1)*iter);    
    flopcount2 = Nflops/time;
    QDPIO::cout << "qdp_lcdotr+gsum time: " << time << " (us), " << flopcount2 << "Mflops per node " << endl;
  }

  QDPIO::cout << "Fn: innerProdReal() flopcount: "  << flopcount1 << " Mflop/s" << endl;
  QDPIO::cout << "Fn: qdp_lcdotr()+gum flopcount: "  << flopcount2 << " Mflop/s" << endl;
}



