# Substitutions
buildCXX = @bagel_CXX@
buildCXXFLAGS = @bagel_CXXFLAGS@
buildLDFLAGS = @bagel_LDFLAGS@
buildLIBS = @bagel_LIBS@
buildLINKFLAGS = $(buildLDFLAGS) $(buildLIBS)

# Override compile and links for build system compilation
SUFFIXES: .C $(EXEEXT)

noinst_PROGRAMS = \
	write_vaxpy$(EXEEXT) \
	write_vaxpby$(EXEEXT) \
	write_vadd$(EXEEXT) \
	write_vscal$(EXEEXT) \
	write_local_sumsq$(EXEEXT) \
	write_local_vcdot$(EXEEXT) \
	write_local_vcdot_real$(EXEEXT) \
	write_vscal_chproj$(EXEEXT) \
	write_vadd_chproj$(EXEEXT) \
	write_vaxpy_chproj$(EXEEXT) \
	write_vxpay_chproj$(EXEEXT) \
	write_vaxpby_chproj$(EXEEXT) \
	write_vscal_g5$(EXEEXT) \
	write_vaxpbg5y$(EXEEXT) \
	write_vag5xpby$(EXEEXT) \
	write_vxpag5y$(EXEEXT) \
	write_vg5axpby$(EXEEXT)

write_vaxpy_SOURCES = write_vaxpy.cc
write_vaxpy$(EXEEXT): write_vaxpy.cc

write_vaxpby_SOURCES = write_vaxpby.cc
write_vaxpby$(EXEEXT): write_vaxpby.cc

write_vadd_SOURCES = write_vadd.cc
write_vadd$(EXEEXT): write_vadd.cc

write_vscal_SOURCES = write_vscal.cc
write_vscal$(EXEEXT): write_vscal.cc

write_local_sumsq_SOURCES = write_local_sumsq.cc
write_local_sumsq$(EXEEXT): write_local_sumsq.cc

write_local_vcdot_SOURCES = write_local_vcdot.cc
write_local_vcdot$(EXEEXT): write_local_vcdot.cc

write_local_vcdot_real_SOURCES = write_local_vcdot_real.cc
write_local_vcdot_real$(EXEEXT): write_local_vcdot_real.cc

write_vscal_chproj_SOURCES = write_vscal_chproj.cc
write_vscal_chproj$(EXEEXT): write_vscal_chproj.cc

write_vadd_chproj_SOURCES = write_vadd_chproj.cc
write_vadd_chproj$(EXEEXT): write_vadd_chproj.cc

write_vaxpy_chproj_SOURCES = write_vaxpy_chproj.cc
write_vaxpy_chproj$(EXEEXT): write_vaxpy_chproj.cc

write_vxpay_chproj_SOURCES = write_vxpay_chproj.cc
write_vxpay_chproj$(EXEEXT): write_vxpay_chproj.cc

write_vaxpby_chproj_SOURCES = write_vaxpby_chproj.cc
write_vaxpby_chproj$(EXEEXT): write_vaxpby_chproj.cc

write_vscal_g5_SOURCES = write_vscal_g5.cc
write_vscal_g5$(EXEEXT): write_vscal_g5.cc

write_vaxpbg5y_SOURCES = write_vaxpbg5y.cc
write_vaxpbg5y$(EXEEXT): write_vaxpbg5y.cc

write_vxpag5y_SOURCES = write_vxpag5y.cc
write_vxpag5y$(EXEEXT): write_vxpag5y.cc

write_vag5xpby_SOURCES = write_vag5xpby.cc
write_vag5xpby$(EXEEXT): write_vag5xpby.cc

write_vg5axpby_SOURCES = write_vg5axpby.cc
write_vg5axpby$(EXEEXT): write_vg5axpby.cc

%$(EXEEXT): %.cc
	$(buildCXX) $(buildCXXFLAGS) -o $@ $< $(buildLDFLAGS) $(buildLIBS)
	chmod u+x $@

