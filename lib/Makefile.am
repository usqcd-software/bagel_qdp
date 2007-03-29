
nobase_include_HEADERS = bagel_qdp.h \
	bagel_qdp_options.h

nodist_include_HEADERS = \
	bagel_qdp_options_internal.h

lib_LIBRARIES = libbagel_qdp.a

libbagel_qdp_a_SOURCES = 

arch_sources =  \
	bagel_vaxpy3.S \
	bagel_vaxmy3.S \
	bagel_vaxpby3.S \
	bagel_vaxmby3.S \
	bagel_vadd3.S \
	bagel_vsub3.S \
	bagel_vscal3.S \
	bagel_local_sumsq.S \
	bagel_local_vcdot.S \
	bagel_local_vcdot_real.S \
	bagel_vscal_ch_plus.S \
	bagel_vscal_ch_minus.S \
	bagel_vadd_ch_plus.S \
	bagel_vadd_ch_minus.S \
	bagel_vsub_ch_plus.S \
	bagel_vsub_ch_minus.S \
	bagel_vaxpy_ch_plus.S \
	bagel_vaxpy_ch_minus.S \
	bagel_vaxmy_ch_plus.S \
	bagel_vaxmy_ch_minus.S \
	bagel_vxpay_ch_plus.S \
	bagel_vxpay_ch_minus.S \
	bagel_vxmay_ch_plus.S \
	bagel_vxmay_ch_minus.S \
	bagel_vaxpby_ch_plus.S \
	bagel_vaxpby_ch_minus.S \
	bagel_vaxmby_ch_plus.S \
	bagel_vaxmby_ch_minus.S \
	bagel_vscal_g5.S \
	bagel_vaxpbg5y.S \
	bagel_vaxmbg5y.S \
	bagel_vag5xpby.S \
	bagel_vag5xmby.S \
	bagel_vxpag5y.S \
	bagel_vxmag5y.S \
	bagel_vg5axpby.S \
	bagel_vg5axmby.S \
	bagel_vaxpbg5iy.S \
	bagel_vaxmbg5iy.S \
	bagel_vxpag5iy.S \
	bagel_vxmag5iy.S

# Noarch architectures
noarch_sources = \
	noarch_bagel_vaxpy3.cc \
	noarch_bagel_vaxmy3.cc \
	noarch_bagel_vaxpby3.cc \
	noarch_bagel_vaxmby3.cc \
	noarch_bagel_vadd3.cc \
	noarch_bagel_vsub3.cc \
	noarch_bagel_vscal3.cc \
	noarch_bagel_local_sumsq.cc \
	noarch_bagel_local_vcdot.cc \
	noarch_bagel_local_vcdot_real.cc  \
	noarch_bagel_vscal_ch_plus.cc \
	noarch_bagel_vscal_ch_minus.cc \
	noarch_bagel_vadd_ch_plus.cc \
	noarch_bagel_vadd_ch_minus.cc \
	noarch_bagel_vsub_ch_plus.cc \
	noarch_bagel_vsub_ch_minus.cc \
	noarch_bagel_vaxpy_ch_plus.cc \
	noarch_bagel_vaxpy_ch_minus.cc \
	noarch_bagel_vaxmy_ch_plus.cc \
	noarch_bagel_vaxmy_ch_minus.cc \
	noarch_bagel_vxpay_ch_plus.cc \
	noarch_bagel_vxpay_ch_minus.cc \
	noarch_bagel_vxmay_ch_plus.cc \
	noarch_bagel_vxmay_ch_minus.cc \
	noarch_bagel_vaxpby_ch_plus.cc \
	noarch_bagel_vaxpby_ch_minus.cc \
	noarch_bagel_vaxmby_ch_plus.cc \
	noarch_bagel_vaxmby_ch_minus.cc \
	noarch_bagel_vscal_g5.cc \
	noarch_bagel_vaxpbg5y.cc \
	noarch_bagel_vaxmbg5y.cc \
	noarch_bagel_vag5xpby.cc \
	noarch_bagel_vag5xmby.cc \
	noarch_bagel_vxpag5y.cc \
	noarch_bagel_vxmag5y.cc \
	noarch_bagel_vg5axpby.cc \
	noarch_bagel_vg5axmby.cc \
	noarch_bagel_vaxpbg5iy.cc \
	noarch_bagel_vaxmbg5iy.cc \
	noarch_bagel_vxpag5iy.cc \
	noarch_bagel_vxmag5iy.cc

# Conditional definition of extra sources
if BUILD_ARCH_SOURCES
extra_sources = $(arch_sources)
endif

# The special noarch architecture 
if BUILD_NOARCH_SOURCES
extra_sources = $(noarch_sources)
endif

# Mark the extra sources for cleaning
CLEANFILES = $(extra_sources)

# Add the extra sources to the library, but not to the dist
nodist_libbagel_qdp_a_SOURCES = $(extra_sources)
libbagel_qdp_a_DEPENDENCIES= $(extra_sources) \
	$(bagelWriterDir)/write_vaxpy$(EXEEXT) \
        $(bagelWriterDir)/write_vaxpby$(EXEEXT) \
	$(bagelWriterDir)/write_vadd$(EXEEXT) \
	$(bagelWriterDir)/write_vscal$(EXEEXT) \
	$(bagelWriterDir)/write_local_sumsq$(EXEEXT) \
	$(bagelWriterDir)/write_local_vcdot$(EXEEXT) \
	$(bagelWriterDir)/write_local_vcdot_real$(EXEEXT) \
	$(bagelWriterDir)/write_vscal_chproj$(EXEEXT) \
	$(bagelWriterDir)/write_vadd_chproj$(EXEEXT) \
	$(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT) \
	$(bagelWriterDir)/write_vxpay_chproj$(EXEEXT) \
	$(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT) \
	$(bagelWriterDir)/write_vscal_g5$(EXEEXT) \
	$(bagelWriterDir)/write_vaxpbg5y$(EXEEXT) \
	$(bagelWriterDir)/write_vag5xpby$(EXEEXT) \
	$(bagelWriterDir)/write_vxpag5y$(EXEEXT) \
	$(bagelWriterDir)/write_vg5axpby$(EXEEXT) \
	$(bagelWriterDir)/write_vaxpbg5iy$(EXEEXT) \
	$(bagelWriterDir)/write_vxpag5iy$(EXEEXT)

# Funky rules for making the .S and .C sources
bagelTargetName = @bagel_MUNGED_TARGET@
bagelWriterDir  = @top_builddir@/code_gen

bagel_vaxpy3.S: $(bagelWriterDir)/write_vaxpy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy$(EXEEXT) -nqdp_vaxpy3 -P$(bagelTargetName) > $@

bagel_vaxmy3.S: $(bagelWriterDir)/write_vaxpy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy$(EXEEXT) -m -nqdp_vaxmy3 -P$(bagelTargetName) > $@

bagel_vaxpby3.S: $(bagelWriterDir)/write_vaxpby$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby$(EXEEXT) -nqdp_vaxpby3 -P$(bagelTargetName) > $@

bagel_vaxmby3.S: $(bagelWriterDir)/write_vaxpby$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby$(EXEEXT) -m -nqdp_vaxmby3 -P$(bagelTargetName) > $@

bagel_vadd3.S: $(bagelWriterDir)/write_vadd$(EXEEXT)
	$(bagelWriterDir)/write_vadd$(EXEEXT) -nqdp_vadd3 -P$(bagelTargetName) > $@

bagel_vsub3.S: $(bagelWriterDir)/write_vadd$(EXEEXT)
	$(bagelWriterDir)/write_vadd$(EXEEXT) -m -nqdp_vsub3 -P$(bagelTargetName) > $@


bagel_vscal3.S: $(bagelWriterDir)/write_vscal$(EXEEXT)
	$(bagelWriterDir)/write_vscal$(EXEEXT) -nqdp_vscal3 -P$(bagelTargetName) > $@

bagel_local_sumsq.S: $(bagelWriterDir)/write_local_sumsq$(EXEEXT)
	$(bagelWriterDir)/write_local_sumsq$(EXEEXT) -nqdp_lsum2 -P$(bagelTargetName) > $@

bagel_local_vcdot.S: $(bagelWriterDir)/write_local_vcdot$(EXEEXT)
	$(bagelWriterDir)/write_local_vcdot$(EXEEXT) -nqdp_lcdot -P$(bagelTargetName) > $@

bagel_local_vcdot_real.S: $(bagelWriterDir)/write_local_vcdot_real$(EXEEXT)
	$(bagelWriterDir)/write_local_vcdot_real$(EXEEXT) -nqdp_lcdotr -P$(bagelTargetName) > $@

bagel_vscal_ch_plus.S: $(bagelWriterDir)/write_vscal_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vscal_chproj$(EXEEXT) -nqdp_vscal_chp -P$(bagelTargetName) -p > $@

bagel_vscal_ch_minus.S: $(bagelWriterDir)/write_vscal_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vscal_chproj$(EXEEXT) -nqdp_vscal_chm -P$(bagelTargetName) > $@


bagel_vadd_ch_plus.S: $(bagelWriterDir)/write_vadd_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vadd_chproj$(EXEEXT)  -nqdp_vadd_chp -P$(bagelTargetName) -p > $@


bagel_vadd_ch_minus.S: $(bagelWriterDir)/write_vadd_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vadd_chproj$(EXEEXT)  -nqdp_vadd_chm -P$(bagelTargetName) > $@

bagel_vsub_ch_plus.S: $(bagelWriterDir)/write_vadd_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vadd_chproj$(EXEEXT)  -m -nqdp_vsub_chp -P$(bagelTargetName) -p > $@

bagel_vsub_ch_minus.S: $(bagelWriterDir)/write_vadd_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vadd_chproj$(EXEEXT) -m -nqdp_vsub_chm -P$(bagelTargetName)  > $@

bagel_vaxpy_ch_plus.S: $(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT) -nqdp_vaxpy_chp -P$(bagelTargetName) -p > $@

bagel_vaxpy_ch_minus.S: $(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT) -nqdp_vaxpy_chm -P$(bagelTargetName) > $@

bagel_vaxmy_ch_plus.S: $(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT) -m -nqdp_vaxmy_chp -P$(bagelTargetName) -p > $@

bagel_vaxmy_ch_minus.S: $(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT) -m -nqdp_vaxmy_chm -P$(bagelTargetName) > $@

bagel_vxpay_ch_plus.S: $(bagelWriterDir)/write_vxpay_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vxpay_chproj$(EXEEXT) -nqdp_vxpay_chp -P$(bagelTargetName) -p > $@

bagel_vxpay_ch_minus.S: $(bagelWriterDir)/write_vxpay_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vxpay_chproj$(EXEEXT) -nqdp_vxpay_chm -P$(bagelTargetName) > $@

bagel_vxmay_ch_plus.S: $(bagelWriterDir)/write_vxpay_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vxpay_chproj$(EXEEXT) -m -nqdp_vxmay_chp -P$(bagelTargetName) -p > $@

bagel_vxmay_ch_minus.S: $(bagelWriterDir)/write_vxpay_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vxpay_chproj$(EXEEXT) -m -nqdp_vxmay_chm -P$(bagelTargetName) > $@


bagel_vaxpby_ch_plus.S: $(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT) -nqdp_vaxpby_chp -P$(bagelTargetName) -p > $@

bagel_vaxpby_ch_minus.S: $(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT) -nqdp_vaxpby_chm -P$(bagelTargetName) > $@

bagel_vaxmby_ch_plus.S: $(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT) -m -nqdp_vaxmby_chp -P$(bagelTargetName) -p > $@

bagel_vaxmby_ch_minus.S: $(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT) -m -nqdp_vaxmby_chm -P$(bagelTargetName) > $@

bagel_vscal_g5.S: $(bagelWriterDir)/write_vscal_g5$(EXEEXT)
	$(bagelWriterDir)/write_vscal_g5$(EXEEXT) -nqdp_vscal_g5 -P$(bagelTargetName) > $@

bagel_vaxpbg5y.S: $(bagelWriterDir)/write_vaxpbg5y$(EXEEXT)
	$(bagelWriterDir)/write_vaxpbg5y$(EXEEXT) -nqdp_vaxpbg5y -P$(bagelTargetName) > $@

bagel_vaxmbg5y.S: $(bagelWriterDir)/write_vaxpbg5y$(EXEEXT)
	$(bagelWriterDir)/write_vaxpbg5y$(EXEEXT) -m -nqdp_vaxmbg5y -P$(bagelTargetName) > $@

bagel_vg5axpby.S: $(bagelWriterDir)/write_vg5axpby$(EXEEXT)
	$(bagelWriterDir)/write_vg5axpby$(EXEEXT) -nqdp_vg5axpby -P$(bagelTargetName) > $@

bagel_vg5axmby.S: $(bagelWriterDir)/write_vg5axpby$(EXEEXT)
	$(bagelWriterDir)/write_vg5axpby$(EXEEXT) -m -nqdp_vg5axmby -P$(bagelTargetName) > $@

bagel_vag5xpby.S: $(bagelWriterDir)/write_vag5xpby$(EXEEXT)
	$(bagelWriterDir)/write_vag5xpby$(EXEEXT) -nqdp_vag5xpby -P$(bagelTargetName) > $@

bagel_vag5xmby.S: $(bagelWriterDir)/write_vag5xpby$(EXEEXT)
	$(bagelWriterDir)/write_vag5xpby$(EXEEXT) -m -nqdp_vag5xmby -P$(bagelTargetName) > $@

bagel_vxpag5y.S: $(bagelWriterDir)/write_vxpag5y$(EXEEXT)
	$(bagelWriterDir)/write_vxpag5y$(EXEEXT) -nqdp_vxpag5y -P$(bagelTargetName) > $@

bagel_vxmag5y.S: $(bagelWriterDir)/write_vxpag5y$(EXEEXT)
	$(bagelWriterDir)/write_vxpag5y$(EXEEXT) -m -nqdp_vxmag5y -P$(bagelTargetName) > $@

bagel_vaxpbg5iy.S: $(bagelWriterDir)/write_vaxpbg5iy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpbg5iy$(EXEEXT) -nqdp_vaxpbg5iy -P$(bagelTargetName) > $@

bagel_vaxmbg5iy.S: $(bagelWriterDir)/write_vaxpbg5iy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpbg5iy$(EXEEXT) -m -nqdp_vaxmbg5iy -P$(bagelTargetName) > $@

bagel_vxpag5iy.S: $(bagelWriterDir)/write_vxpag5iy$(EXEEXT)
	$(bagelWriterDir)/write_vxpag5iy$(EXEEXT) -nqdp_vxpag5iy -P$(bagelTargetName) > $@

bagel_vxmag5iy.S: $(bagelWriterDir)/write_vxpag5iy$(EXEEXT)
	$(bagelWriterDir)/write_vxpag5iy$(EXEEXT) -m -nqdp_vxmag5iy -P$(bagelTargetName) > $@

bagel_su3_mm.S: $(bagelWriterDir)/write_su3_mm$(EXEEXT)
	$(bagelWriterDir)/write_su3_mm$(EXEEXT) -nqdp_su3_mm -P$(bagelTargetName) > $@

########## NOARCH TARGETS ##############################

noarch_bagel_vaxpy3.cc: $(bagelWriterDir)/write_vaxpy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy$(EXEEXT) -nqdp_vaxpy3 -P$(bagelTargetName) > $@

noarch_bagel_vaxmy3.cc: $(bagelWriterDir)/write_vaxpy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy$(EXEEXT) -m -nqdp_vaxmy3 -P$(bagelTargetName) > $@
noarch_bagel_vaxpby3.cc: $(bagelWriterDir)/write_vaxpby$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby$(EXEEXT) -nqdp_vaxpby3 -P$(bagelTargetName) > $@

noarch_bagel_vaxmby3.cc: $(bagelWriterDir)/write_vaxpby$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby$(EXEEXT) -m -nqdp_vaxmby3 -P$(bagelTargetName) > $@

noarch_bagel_vadd3.cc: $(bagelWriterDir)/write_vadd$(EXEEXT)
	$(bagelWriterDir)/write_vadd$(EXEEXT) -nqdp_vadd3 -P$(bagelTargetName) > $@

noarch_bagel_vsub3.cc: $(bagelWriterDir)/write_vadd$(EXEEXT)
	$(bagelWriterDir)/write_vadd$(EXEEXT) -m -nqdp_vsub3 -P$(bagelTargetName) > $@

noarch_bagel_vscal3.cc: $(bagelWriterDir)/write_vscal$(EXEEXT)
	$(bagelWriterDir)/write_vscal$(EXEEXT) -nqdp_vscal3 -P$(bagelTargetName) > $@

noarch_bagel_local_sumsq.cc: $(bagelWriterDir)/write_local_sumsq$(EXEEXT)
	$(bagelWriterDir)/write_local_sumsq$(EXEEXT) -nqdp_lsum2 -P$(bagelTargetName) > $@

noarch_bagel_local_vcdot.cc: $(bagelWriterDir)/write_local_vcdot$(EXEEXT)
	$(bagelWriterDir)/write_local_vcdot$(EXEEXT) -nqdp_lcdot -P$(bagelTargetName) > $@


noarch_bagel_local_vcdot_real.cc: $(bagelWriterDir)/write_local_vcdot_real$(EXEEXT)
	$(bagelWriterDir)/write_local_vcdot_real$(EXEEXT) -nqdp_lcdotr -P$(bagelTargetName) > $@

noarch_bagel_vscal_ch_plus.cc: $(bagelWriterDir)/write_vscal_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vscal_chproj$(EXEEXT) -nqdp_vscal_chp -P$(bagelTargetName) -p > $@

noarch_bagel_vscal_ch_minus.cc: $(bagelWriterDir)/write_vscal_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vscal_chproj$(EXEEXT) -nqdp_vscal_chm -P$(bagelTargetName) > $@

noarch_bagel_vadd_ch_plus.cc: $(bagelWriterDir)/write_vadd_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vadd_chproj$(EXEEXT) -nqdp_vadd_chp -P$(bagelTargetName) -p > $@


noarch_bagel_vadd_ch_minus.cc: $(bagelWriterDir)/write_vadd_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vadd_chproj$(EXEEXT)  -nqdp_vadd_chm -P$(bagelTargetName) > $@

noarch_bagel_vsub_ch_plus.cc: $(bagelWriterDir)/write_vadd_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vadd_chproj$(EXEEXT) -m -nqdp_vsub_chp -P$(bagelTargetName) -p > $@

noarch_bagel_vsub_ch_minus.cc: $(bagelWriterDir)/write_vadd_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vadd_chproj$(EXEEXT) -m -nqdp_vsub_chm -P$(bagelTargetName) > $@

noarch_bagel_vaxpy_ch_plus.cc: $(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT) -nqdp_vaxpy_chp -P$(bagelTargetName) -p > $@

noarch_bagel_vaxpy_ch_minus.cc: $(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT) -nqdp_vaxpy_chm -P$(bagelTargetName) > $@

noarch_bagel_vaxmy_ch_plus.cc: $(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT) -m -nqdp_vaxmy_chp -P$(bagelTargetName) -p > $@

noarch_bagel_vaxmy_ch_minus.cc: $(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpy_chproj$(EXEEXT) -m -nqdp_vaxmy_chm -P$(bagelTargetName) > $@

noarch_bagel_vxpay_ch_plus.cc: $(bagelWriterDir)/write_vxpay_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vxpay_chproj$(EXEEXT) -nqdp_vxpay_chp -P$(bagelTargetName) -p > $@

noarch_bagel_vxpay_ch_minus.cc: $(bagelWriterDir)/write_vxpay_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vxpay_chproj$(EXEEXT) -nqdp_vxpay_chm -P$(bagelTargetName) > $@

noarch_bagel_vxmay_ch_plus.cc: $(bagelWriterDir)/write_vxpay_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vxpay_chproj$(EXEEXT) -m -nqdp_vxmay_chp -P$(bagelTargetName) -p > $@

noarch_bagel_vxmay_ch_minus.cc: $(bagelWriterDir)/write_vxpay_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vxpay_chproj$(EXEEXT) -m -nqdp_vxmay_chm -P$(bagelTargetName) > $@

noarch_bagel_vaxpby_ch_plus.cc: $(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT) -nqdp_vaxpby_chp -P$(bagelTargetName) -p > $@

noarch_bagel_vaxpby_ch_minus.cc: $(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT) -nqdp_vaxpby_chm -P$(bagelTargetName) > $@

noarch_bagel_vaxmby_ch_plus.cc: $(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT) -m -nqdp_vaxmby_chp -P$(bagelTargetName) -p > $@

noarch_bagel_vaxmby_ch_minus.cc: $(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT)
	$(bagelWriterDir)/write_vaxpby_chproj$(EXEEXT) -m -nqdp_vaxmby_chm -P$(bagelTargetName) > $@

noarch_bagel_vscal_g5.cc: $(bagelWriterDir)/write_vscal_g5$(EXEEXT)
	$(bagelWriterDir)/write_vscal_g5$(EXEEXT) -nqdp_vscal_g5 -P$(bagelTargetName) > $@

noarch_bagel_vaxpbg5y.cc: $(bagelWriterDir)/write_vaxpbg5y$(EXEEXT)
	$(bagelWriterDir)/write_vaxpbg5y$(EXEEXT) -nqdp_vaxpbg5y -P$(bagelTargetName) > $@

noarch_bagel_vaxmbg5y.cc: $(bagelWriterDir)/write_vaxpbg5y$(EXEEXT)
	$(bagelWriterDir)/write_vaxpbg5y$(EXEEXT) -m -nqdp_vaxmbg5y -P$(bagelTargetName) > $@

noarch_bagel_vxpag5y.cc: $(bagelWriterDir)/write_vxpag5y$(EXEEXT)
	$(bagelWriterDir)/write_vxpag5y$(EXEEXT) -nqdp_vxpag5y -P$(bagelTargetName) > $@

noarch_bagel_vxmag5y.cc: $(bagelWriterDir)/write_vxpag5y$(EXEEXT)
	$(bagelWriterDir)/write_vxpag5y$(EXEEXT) -m -nqdp_vxmag5y -P$(bagelTargetName) > $@

noarch_bagel_vag5xpby.cc: $(bagelWriterDir)/write_vag5xpby$(EXEEXT)
	$(bagelWriterDir)/write_vag5xpby$(EXEEXT) -nqdp_vag5xpby -P$(bagelTargetName) > $@

noarch_bagel_vag5xmby.cc: $(bagelWriterDir)/write_vag5xpby$(EXEEXT)
	$(bagelWriterDir)/write_vag5xpby$(EXEEXT) -m -nqdp_vag5xmby -P$(bagelTargetName) > $@

noarch_bagel_vg5axpby.cc: $(bagelWriterDir)/write_vg5axpby$(EXEEXT)
	$(bagelWriterDir)/write_vg5axpby$(EXEEXT) -nqdp_vg5axpby -P$(bagelTargetName) > $@

noarch_bagel_vg5axmby.cc: $(bagelWriterDir)/write_vg5axpby$(EXEEXT)
	$(bagelWriterDir)/write_vg5axpby$(EXEEXT) -m -nqdp_vg5axmby -P$(bagelTargetName) > $@

noarch_bagel_vaxpbg5iy.cc: $(bagelWriterDir)/write_vaxpbg5iy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpbg5iy$(EXEEXT) -nqdp_vaxpbg5iy -P$(bagelTargetName) > $@

noarch_bagel_vaxmbg5iy.cc: $(bagelWriterDir)/write_vaxpbg5iy$(EXEEXT)
	$(bagelWriterDir)/write_vaxpbg5iy$(EXEEXT) -m -nqdp_vaxmbg5iy -P$(bagelTargetName) > $@

noarch_bagel_vxpag5iy.cc: $(bagelWriterDir)/write_vxpag5iy$(EXEEXT)
	$(bagelWriterDir)/write_vxpag5iy$(EXEEXT) -nqdp_vxpag5iy -P$(bagelTargetName) > $@

noarch_bagel_vxmag5iy.cc: $(bagelWriterDir)/write_vxpag5iy$(EXEEXT)
	$(bagelWriterDir)/write_vxpag5iy$(EXEEXT) -m -nqdp_vxmag5iy -P$(bagelTargetName) > $@

noarch_bagel_su3_mm.cc: $(bagelWriterDir)/write_su3_mm$(EXEEXT)
	$(bagelWriterDir)/write_su3_mm$(EXEEXT) -nqdp_su3_mm -P$(bagelTargetName) > $@
