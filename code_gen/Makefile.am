# Substitutions
buildCXX = @bagel_CXX@
buildCXXFLAGS = @bagel_CXXFLAGS@
buildLDFLAGS = @bagel_LDFLAGS@
buildLIBS = @bagel_LIBS@
buildLINKFLAGS = $(buildLDFLAGS) $(buildLIBS)

# Override compile and links for build system compilation
SUFFIXES: .C .x

noinst_PROGRAMS = \
	write_vaxpy.x \
	write_vaxpby.x \
	write_vadd.x \
	write_vscal.x \
	write_local_sumsq.x \
	write_local_vcdot.x \
	write_local_vcdot_real.x \
	write_vscal_chproj.x \
	write_vadd_chproj.x \
	write_vaxpy_chproj.x \
	write_vxpay_chproj.x \
	write_vaxpby_chproj.x \
	write_vscal_g5.x \
	write_vaxpbg5y.x \
	write_vag5xpby.x \
	write_vxpag5y.x \
	write_vg5axpby.x

write_vaxpy_x_SOURCES = write_vaxpy.cc
write_vaxpy.x: write_vaxpy.cc

write_vaxpby_x_SOURCES = write_vaxpby.cc
write_vaxpby.x: write_vaxpby.cc

write_vadd_x_SOURCES = write_vadd.cc
write_vadd.x: write_vadd.cc

write_vscal_x_SOURCES = write_vscal.cc
write_vscal.x: write_vscal.cc

write_local_sumsq_x_SOURCES = write_local_sumsq.cc
write_local_sumsq.x: write_local_sumsq.cc

write_local_vcdot_x_SOURCES = write_local_vcdot.cc
write_local_vcdot.x: write_local_vcdot.cc

write_local_vcdot_real_x_SOURCES = write_local_vcdot_real.cc
write_local_vcdot_real.x: write_local_vcdot_real.cc

write_vscal_chproj_x_SOURCES = write_vscal_chproj.cc
write_vscal_chproj.x: write_vscal_chproj.cc

write_vadd_chproj_x_SOURCES = write_vadd_chproj.cc
write_vadd_chproj.x: write_vadd_chproj.cc

write_vaxpy_chproj_x_SOURCES = write_vaxpy_chproj.cc
write_vaxpy_chproj.x: write_vaxpy_chproj.cc

write_vxpay_chproj_x_SOURCES = write_vxpay_chproj.cc
write_vxpay_chproj.x: write_vxpay_chproj.cc

write_vaxpby_chproj_x_SOURCES = write_vaxpby_chproj.cc
write_vaxpby_chproj.x: write_vaxpby_chproj.cc

write_vscal_g5_x_SOURCES = write_vscal_g5.cc
write_vscal_g5.x: write_vscal_g5.cc

write_vaxpbg5y_x_SOURCES = write_vaxpbg5y.cc
write_vaxpbg5y.x: write_vaxpbg5y.cc

write_vxpag5y_x_SOURCES = write_vxpag5y.cc
write_vxpag5y.x: write_vxpag5y.cc

write_vag5xpby_x_SOURCES = write_vag5xpby.cc
write_vag5xpby.x: write_vag5xpby.cc

write_vg5axpby_x_SOURCES = write_vg5axpby.cc
write_vg5axpby.x: write_vg5axpby.cc

%.x: %.cc
	$(buildCXX) $(buildCXXFLAGS) -o $@ $< $(buildLDFLAGS) $(buildLIBS)
	chmod u+x $@

