#ifndef BAGEL_QDP_H
#define BAGEL_QDP_H

#ifdef __cplusplus




extern "C" { 
#include "bagel_qdp_options.h"

  void qdp_vaxpy3(BAGELQDPFloat* Out, 
		  BAGELQDPFloat* scalep, 
		  BAGELQDPFloat* InScale,
		  BAGELQDPFloat* Out,
		  int n_3vec);

  void qdp_vaxmy3(BAGELQDPFloat* Out, 
		  BAGELQDPFloat* scalep, 
		  BAGELQDPFloat* InScale,
		  BAGELQDPFloat* Out,
		  int n_3vec);

  void qdp_vaxpby3(BAGELQDPFloat* Out, 
		   BAGELQDPFloat* a, 
		   BAGELQDPFloat* x,
		   BAGELQDPFloat* b,
		   BAGELQDPFloat* y,
		   int n_3vec);

  void qdp_vaxmby3(BAGELQDPFloat* Out, 
		   BAGELQDPFloat* a,
		   BAGELQDPFloat* x,
		   BAGELQDPFloat* b,
		   BAGELQDPFloat* y,
		   int n_3vec);

  void qdp_vadd3(BAGELQDPFloat* Out, 
		 BAGELQDPFloat* x,
		 BAGELQDPFloat* y,
		 int n_3vec);

  void qdp_vsub3(BAGELQDPFloat* Out, 
		 BAGELQDPFloat* x,
		 BAGELQDPFloat* y,
		 int n_3vec);

  void qdp_vscal3(BAGELQDPFloat* Out,
		  BAGELQDPFloat* a,
		  BAGELQDPFloat* x,
		  int n_3vec);

  void qdp_lsum2(BAGELQDPFloat* Out,
		 BAGELQDPFloat* x,
		 int n_3vec);

  void qdp_lcdot(BAGELQDPFloat* Out_r,
		 BAGELQDPFloat* Out_i,
		 BAGELQDPFloat* x,
		 BAGELQDPFloat* y,
		 int n_3vec);

  void qdp_lcdotr(BAGELQDPFloat* Out,
		  BAGELQDPFloat* x,
		  BAGELQDPFloat* y,
		  int n_3vec);

  void qdp_vscal_chp(BAGELQDPFloat* Out,
		     BAGELQDPFloat* a,
		     BAGELQDPFloat* x,
		     int n_4vec);
  
  void qdp_vscal_chm(BAGELQDPFloat* Out,
		     BAGELQDPFloat* a,
		     BAGELQDPFloat* x,
		     int n_4vec);

  void qdp_vadd_chp(BAGELQDPFloat* Out, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* y,
		    int n_4vec);

  void qdp_vadd_chm(BAGELQDPFloat* Out, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* y,
		    int n_4vec);

  void qdp_vsub_chp(BAGELQDPFloat* Out, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* y,
		    int n_4vec);

  void qdp_vsub_chm(BAGELQDPFloat* Out, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* y,
		    int n_4vec);

  void qdp_vaxpy_chp(BAGELQDPFloat* Out, 
		     BAGELQDPFloat* a,
		     BAGELQDPFloat* x,
		     BAGELQDPFloat* y,
		     int n_4vec);

  void qdp_vaxpy_chm(BAGELQDPFloat* Out,
		     BAGELQDPFloat* a, 
		     BAGELQDPFloat* x,
		     BAGELQDPFloat* y,
		     int n_4vec);

  void qdp_vaxmy_chp(BAGELQDPFloat* Out,
		     BAGELQDPFloat* a, 
		     BAGELQDPFloat* x,
		     BAGELQDPFloat* y,
		     int n_4vec);
  
  void qdp_vaxmy_chm(BAGELQDPFloat* Out,
		     BAGELQDPFloat* a, 
		     BAGELQDPFloat* x,
		     BAGELQDPFloat* y,
		     int n_4vec);

  void qdp_vxpay_chp(BAGELQDPFloat* Out, 
		     BAGELQDPFloat* a,
		     BAGELQDPFloat* x,
		     BAGELQDPFloat* y,
		     int n_4vec);

  void qdp_vxpay_chm(BAGELQDPFloat* Out,
		     BAGELQDPFloat* a, 
		     BAGELQDPFloat* x,
		     BAGELQDPFloat* y,
		     int n_4vec);

  void qdp_vxmay_chp(BAGELQDPFloat* Out,
		     BAGELQDPFloat* a, 
		     BAGELQDPFloat* x,
		     BAGELQDPFloat* y,
		     int n_4vec);
  
  void qdp_vxmay_chm(BAGELQDPFloat* Out,
		     BAGELQDPFloat* a, 
		     BAGELQDPFloat* x,
		     BAGELQDPFloat* y,
		     int n_4vec);

  void qdp_vaxpby_chp(BAGELQDPFloat* Out, 
		      BAGELQDPFloat* a,
		      BAGELQDPFloat* x,
		      BAGELQDPFloat* b,
		      BAGELQDPFloat* y,
		      int n_4vec);

  void qdp_vaxpby_chm(BAGELQDPFloat* Out,
		      BAGELQDPFloat* a, 
		      BAGELQDPFloat* x,
		      BAGELQDPFloat* b,
		      BAGELQDPFloat* y,
		      int n_4vec);

  void qdp_vaxmby_chp(BAGELQDPFloat* Out,
		      BAGELQDPFloat* a, 
		      BAGELQDPFloat* x,
		      BAGELQDPFloat* b,
		      BAGELQDPFloat* y,
		      int n_4vec);
  
  void qdp_vaxmby_chm(BAGELQDPFloat* Out,
		      BAGELQDPFloat* a, 
		      BAGELQDPFloat* x,
		      BAGELQDPFloat* b,
		      BAGELQDPFloat* y,
		      int n_4vec);

  void qdp_vscal_g5(BAGELQDPFloat* Out,
		    BAGELQDPFloat* a,
		    BAGELQDPFloat* x,
		    int n_4vec);

  void qdp_vaxpbg5y(BAGELQDPFloat* Out,
		    BAGELQDPFloat* a, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* b,
		    BAGELQDPFloat* y,
		    int n_4vec);

  void qdp_vaxmbg5y(BAGELQDPFloat* Out,
		    BAGELQDPFloat* a, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* b,
		    BAGELQDPFloat* y,
		    int n_4vec);

  void qdp_vag5xpby(BAGELQDPFloat* Out,
		    BAGELQDPFloat* a, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* b,
		    BAGELQDPFloat* y,
		    int n_4vec);

  void qdp_vag5xmby(BAGELQDPFloat* Out,
		    BAGELQDPFloat* a, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* b,
		    BAGELQDPFloat* y,
		    int n_4vec);

  void qdp_vg5axpby(BAGELQDPFloat* Out,
		    BAGELQDPFloat* a, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* b,
		    BAGELQDPFloat* y,
		    int n_4vec);

  void qdp_vg5axmby(BAGELQDPFloat* Out,
		    BAGELQDPFloat* a, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* b,
		    BAGELQDPFloat* y,
		    int n_4vec);

 void qdp_vxpag5y(BAGELQDPFloat* Out, 
		  BAGELQDPFloat* a,
		  BAGELQDPFloat* x,
		  BAGELQDPFloat* y,
		  int n_4vec);

 void qdp_vxmag5y(BAGELQDPFloat* Out, 
		  BAGELQDPFloat* a,
		  BAGELQDPFloat* x,
		  BAGELQDPFloat* y,
		  int n_4vec);

  void qdp_vaxpbg5iy(BAGELQDPFloat* Out,
		    BAGELQDPFloat* a, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* b,
		    BAGELQDPFloat* y,
		    int n_4vec);

  void qdp_vaxmbg5iy(BAGELQDPFloat* Out,
		    BAGELQDPFloat* a, 
		    BAGELQDPFloat* x,
		    BAGELQDPFloat* b,
		    BAGELQDPFloat* y,
		    int n_4vec);

 void qdp_vxpag5iy(BAGELQDPFloat* Out, 
		   BAGELQDPFloat* a,
		   BAGELQDPFloat* x,
		   BAGELQDPFloat* y,
		   int n_4vec);

 void qdp_vxmag5iy(BAGELQDPFloat* Out, 
		   BAGELQDPFloat* a,
		   BAGELQDPFloat* x,
		   BAGELQDPFloat* y,
		   int n_4vec);

  // A = B*C  - A, B, C are vectors of SU(3) matrices
  //          - one_minus_i_ptr is ignored
  //
  void qdp_su3_mm(BAGELQDPFloat* A,
		  BAGELQDPFloat* B,
		  BAGELQDPFloat* C,
		  unsigned long n_mat,
		  unsigned long one_minus_i_ptr);

  // A = B*(C)^\dagger - A, B, C are vectors of SU(3) matrices
  //                   - one_minus_i_ptr is ignored
  //
  void qdp_su3_ma(BAGELQDPFloat* A,
		  BAGELQDPFloat* B,
		  BAGELQDPFloat* C,
		  unsigned long n_mat,
		  unsigned long one_minus_i_ptr);

  // A = (B)^\dagger * C - A, B, C are vectors of SU(3) matrices
  //                   - one_minus_i_ptr is ignored
  //
  void qdp_su3_am(BAGELQDPFloat* A,
		  BAGELQDPFloat* B,
		  BAGELQDPFloat* C,
		  unsigned long n_mat,
		  unsigned long one_minus_i_ptr);

  // A = (B)^\dagger * (C)^\dagger
  //
  //                 - A, B, C are vectors of SU(3) matrices
  //                 - one_minus_i_ptr MUST point to 
  //                   an array of 2 doubles containing (1,i)
  //                
  //                   ie: one_minus_i_ptr must point to
  //                      double one_minus_i[2] QDP_ALIGN16;
  //                      one_minus_i[0] = (double)1;
  //                      one_minus_i[1] = (double)(-1);
  //
  //                   This is used to conjugate B and C
  //                   There may be a better way, but I am not sure
  //                   how to do it without modifying BAGEL
  //
  void qdp_su3_aa(BAGELQDPFloat* A,
		  BAGELQDPFloat* B,
		  BAGELQDPFloat* C,
		  unsigned long n_mat,
		  unsigned long one_plus_minus_i_ptr);
  

  // A = D + a B*C  - A, B, C are vectors of SU(3) matrices
  //             - a is a scalar float (eg -1)
  //             - one_minus_i_ptr is ignored
  //
  void qdp_su3_mm_peq(BAGELQDPFloat* A,
		      BAGELQDPFloat* a,
		      BAGELQDPFloat* B,
		      BAGELQDPFloat* C,
			unsigned long n_mat,
			unsigned long one_plus_minus_i_ptr);

  // A = D + a B*(C)^\dagger - A, B, C are vectors of SU(3) matrices
  //                      - a is a scalar float: eg -1
  //                      - one_minus_i_ptr is ignored
  //
  void qdp_su3_ma_peq(BAGELQDPFloat* A,
		      BAGELQDPFloat* a,
		      BAGELQDPFloat* B,
		      BAGELQDPFloat* C,
		      unsigned long n_mat,
		      unsigned long one_plus_minus_i_ptr);

  // A = D + a (B)^\dagger * C - A, B, C are vectors of SU(3) matrices
  //                        - a is a scalar float (eg -1)
  //                         - one_minus_i_ptr is ignored
  //
  void qdp_su3_am_peq(BAGELQDPFloat* A,
		      BAGELQDPFloat* a,
		      BAGELQDPFloat* B,
		      BAGELQDPFloat* C,
		      unsigned long n_mat,
		      unsigned long one_plus_minus_i_ptr);

  // A = D + a (B)^\dagger * (C)^\dagger
  //
  //                 - A, B, C are vectors of SU(3) matrices
  //                 - a is a scalar float (eg -1)
  //                 - one_minus_i_ptr MUST point to 
  //                   an array of 2 doubles containing (1,i)
  //                
  //                   ie: one_minus_i_ptr must point to
  //                      double one_minus_i[2] QDP_ALIGN16;
  //                      one_minus_i[0] = (double)1;
  //                      one_minus_i[1] = (double)(-1);
  //
  //                   This is used to conjugate B and C
  //                   There may be a better way, but I am not sure
  //                   how to do it without modifying BAGEL
  //
  void qdp_su3_aa_peq(BAGELQDPFloat* A,
		      BAGELQDPFloat* a,
		      BAGELQDPFloat* B,
		      BAGELQDPFloat* C,
		      unsigned long n_mat,
		      unsigned long one_plus_minus_i_ptr);
};


#endif
#endif
