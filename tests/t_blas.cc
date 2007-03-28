// $Id: t_blas.cc,v 1.3 2007-03-28 17:19:52 bjoo Exp $

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
  Internal::broadcast(iter);


  swatch.reset();
  swatch.start();
  for(int i=0; i < iter; i++) {
    qdp_su3_mm((REAL*)&(a1.elem(0).elem().elem(0,0).real()),
	       (REAL*)&(b.elem(0).elem().elem(0,0).real()),
	       (REAL*)&(c.elem(0).elem().elem(0,0).real()),
	       all.end()-all.start()+1);
    
  }
  swatch.stop(); 
  seconds = swatch.getTimeInSeconds();
  sum(seconds);
  seconds /= Layout::numNodes();

  double flops = iter*180*Layout::sitesOnNode()/1.0e6;

  QDPIO::cout << "New Way: " << seconds << " (s)  " << flops << " MFLOP,  " << flops/seconds << " MFlop/s per process" << endl;

  swatch.reset();
  swatch.start();
  for(int i=0; i < iter; i++) {
    a=b*c;
  }
  swatch.stop(); 
  seconds = swatch.getTimeInSeconds();
  sum(seconds);
  seconds /= Layout::numNodes();

  QDPIO::cout << "Old Way: " << seconds << " (s)  " << flops << " MFLOP,  " << flops/seconds << " MFlop/s per process" << endl;

  // Time to bolt
  QDP_finalize();

  exit(0);
}




