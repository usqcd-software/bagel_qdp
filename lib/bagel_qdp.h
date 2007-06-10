#ifndef BAGEL_QDP_H
#define BAGEL_QDP_H

#ifdef __cplusplus




extern "C" { 
#include "bagel_qdp_options.h"

  void qdp_vaxpy3(Float* Out, 
		  Float* scalep, 
		  Float* InScale,
		  Float* Out,
		  int n_3vec);

  void qdp_vaxmy3(Float* Out, 
		  Float* scalep, 
		  Float* InScale,
		  Float* Out,
		  int n_3vec);

  void qdp_vaxpby3(Float* Out, 
		   Float* a, 
		   Float* x,
		   Float* b,
		   Float* y,
		   int n_3vec);

  void qdp_vaxmby3(Float* Out, 
		   Float* a,
		   Float* x,
		   Float* b,
		   Float* y,
		   int n_3vec);

  void qdp_vadd3(Float* Out, 
		 Float* x,
		 Float* y,
		 int n_3vec);

  void qdp_vsub3(Float* Out, 
		 Float* x,
		 Float* y,
		 int n_3vec);

  void qdp_vscal3(Float* Out,
		  Float* a,
		  Float* x,
		  int n_3vec);

  void qdp_lsum2(Float* Out,
		 Float* x,
		 int n_3vec);

  void qdp_lcdot(Float* Out_r,
		 Float* Out_i,
		 Float* x,
		 Float* y,
		 int n_3vec);

  void qdp_lcdotr(Float* Out,
		  Float* x,
		  Float* y,
		  int n_3vec);

  void qdp_vscal_chp(Float* Out,
		     Float* a,
		     Float* x,
		     int n_4vec);
  
  void qdp_vscal_chm(Float* Out,
		     Float* a,
		     Float* x,
		     int n_4vec);

  void qdp_vadd_chp(Float* Out, 
		    Float* x,
		    Float* y,
		    int n_4vec);

  void qdp_vadd_chm(Float* Out, 
		    Float* x,
		    Float* y,
		    int n_4vec);

  void qdp_vsub_chp(Float* Out, 
		    Float* x,
		    Float* y,
		    int n_4vec);

  void qdp_vsub_chm(Float* Out, 
		    Float* x,
		    Float* y,
		    int n_4vec);

  void qdp_vaxpy_chp(Float* Out, 
		     Float* a,
		     Float* x,
		     Float* y,
		     int n_4vec);

  void qdp_vaxpy_chm(Float* Out,
		     Float* a, 
		     Float* x,
		     Float* y,
		     int n_4vec);

  void qdp_vaxmy_chp(Float* Out,
		     Float* a, 
		     Float* x,
		     Float* y,
		     int n_4vec);
  
  void qdp_vaxmy_chm(Float* Out,
		     Float* a, 
		     Float* x,
		     Float* y,
		     int n_4vec);

  void qdp_vxpay_chp(Float* Out, 
		     Float* a,
		     Float* x,
		     Float* y,
		     int n_4vec);

  void qdp_vxpay_chm(Float* Out,
		     Float* a, 
		     Float* x,
		     Float* y,
		     int n_4vec);

  void qdp_vxmay_chp(Float* Out,
		     Float* a, 
		     Float* x,
		     Float* y,
		     int n_4vec);
  
  void qdp_vxmay_chm(Float* Out,
		     Float* a, 
		     Float* x,
		     Float* y,
		     int n_4vec);

  void qdp_vaxpby_chp(Float* Out, 
		      Float* a,
		      Float* x,
		      Float* b,
		      Float* y,
		      int n_4vec);

  void qdp_vaxpby_chm(Float* Out,
		      Float* a, 
		      Float* x,
		      Float* b,
		      Float* y,
		      int n_4vec);

  void qdp_vaxmby_chp(Float* Out,
		      Float* a, 
		      Float* x,
		      Float* b,
		      Float* y,
		      int n_4vec);
  
  void qdp_vaxmby_chm(Float* Out,
		      Float* a, 
		      Float* x,
		      Float* b,
		      Float* y,
		      int n_4vec);

  void qdp_vscal_g5(Float* Out,
		    Float* a,
		    Float* x,
		    int n_4vec);

  void qdp_vaxpbg5y(Float* Out,
		    Float* a, 
		    Float* x,
		    Float* b,
		    Float* y,
		    int n_4vec);

  void qdp_vaxmbg5y(Float* Out,
		    Float* a, 
		    Float* x,
		    Float* b,
		    Float* y,
		    int n_4vec);

  void qdp_vag5xpby(Float* Out,
		    Float* a, 
		    Float* x,
		    Float* b,
		    Float* y,
		    int n_4vec);

  void qdp_vag5xmby(Float* Out,
		    Float* a, 
		    Float* x,
		    Float* b,
		    Float* y,
		    int n_4vec);

  void qdp_vg5axpby(Float* Out,
		    Float* a, 
		    Float* x,
		    Float* b,
		    Float* y,
		    int n_4vec);

  void qdp_vg5axmby(Float* Out,
		    Float* a, 
		    Float* x,
		    Float* b,
		    Float* y,
		    int n_4vec);

 void qdp_vxpag5y(Float* Out, 
		  Float* a,
		  Float* x,
		  Float* y,
		  int n_4vec);

 void qdp_vxmag5y(Float* Out, 
		  Float* a,
		  Float* x,
		  Float* y,
		  int n_4vec);

  void qdp_vaxpbg5iy(Float* Out,
		    Float* a, 
		    Float* x,
		    Float* b,
		    Float* y,
		    int n_4vec);

  void qdp_vaxmbg5iy(Float* Out,
		    Float* a, 
		    Float* x,
		    Float* b,
		    Float* y,
		    int n_4vec);

 void qdp_vxpag5iy(Float* Out, 
		   Float* a,
		   Float* x,
		   Float* y,
		   int n_4vec);

 void qdp_vxmag5iy(Float* Out, 
		   Float* a,
		   Float* x,
		   Float* y,
		   int n_4vec);

  void qdp_su3_mm(Float* A,
		  Float* B,
		  Float* C,
		  unsigned long n_mat,
		  unsigned long one_plus_minus_i_ptr);

  void qdp_su3_ma(Float* A,
		  Float* B,
		  Float* C,
		  unsigned long n_mat,
		  unsigned long one_plus_minus_i_ptr);

  void qdp_su3_am(Float* A,
		  Float* B,
		  Float* C,
		  unsigned long n_mat,
		  unsigned long one_plus_minus_i_ptr);

  void qdp_su3_aa(Float* A,
		  Float* B,
		  Float* C,
		  unsigned long n_mat,
		  unsigned long one_plus_minus_i_ptr);
};


#endif
#endif
