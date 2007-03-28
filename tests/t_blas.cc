// $Id: t_blas.cc,v 1.2 2007-03-28 16:22:50 bjoo Exp $

#include <iostream>
#include <iomanip>
#include <cstdio>

#include <time.h>

#include "qdp.h"
#include <bagel_qdp.h>

using namespace std;
using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);
  
  // Setup the layout
  // const int foo[] = {8,8,16,16};
  const int local_vol[] = { 2,2,2,2 };
  multi1d<int> nrow(Nd);
  nrow=local_vol;

  Layout::setLattSize(nrow);
  Layout::create();

  LatticeColorMatrix a,b,c;
  LatticeColorMatrix a1;
  LatticeColorMatrix diff;

  gaussian(b);
  gaussian(c);

  QDPIO::cout << "Doing Old" << endl;
  a=b*c;
  QDPIO::cout << "Doing New" << endl;

  qdp_su3_mm((REAL*)&(a1.elem(0).elem().elem(0,0).real()),
	    (REAL*)&(b.elem(0).elem().elem(0,0).real()),
	    (REAL*)&(c.elem(0).elem().elem(0,0).real()),
	    all.end()-all.start()+1);

  
  diff = a1 - a;
  QDPIO::cout << "Diff : " << sqrt(norm2(diff)) << endl;

  // Time it: 
  StopWatch swatch;
  swatch.reset();
  double seconds;
  int iter = 1;
  do {
    swatch.start();
    for(int i=0; i < iter; i++) {
      qdp_su3_mm((REAL*)&(a1.elem(0).elem().elem(0,0).real()),
		 (REAL*)&(b.elem(0).elem().elem(0,0).real()),
		 (REAL*)&(c.elem(0).elem().elem(0,0).real()),
		 all.end()-all.start()+1);
    }
    swatch.stop();
    seconds = swatch.getTimeInSeconds();
    swatch.reset();
    if (seconds < 1 ) iter *=2;
  }
  while( seconds < 1 );


  swatch.reset();
  swatch.start();
  for(int i=0; i < iter; i++) {
    qdp_su3_mm((REAL*)&(a1.elem(0).elem().elem(0,0).real()),
	       (REAL*)&(b.elem(0).elem().elem(0,0).real()),
	       (REAL*)&(c.elem(0).elem().elem(0,0).real()),
	       all.end()-all.start()+1);
    
  }
  swatch.stop(); 
  
  QDPIO::cout << "NEW Iterations :" << iter << " Time = " << swatch.getTimeInSeconds() << " s" << endl;

  swatch.reset();
  swatch.start();
  for(int i=0; i < iter; i++) {
    a=b*c;
  }
  swatch.stop(); 
  
  QDPIO::cout << "OLD Iterations :" << iter << " Time = " << swatch.getTimeInSeconds() << " s" << endl;

  // Time to bolt
  QDP_finalize();

  exit(0);
}




