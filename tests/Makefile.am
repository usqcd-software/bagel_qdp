TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

AM_CXX = @QDPXX_CXX@
INCFLAGS= -I$(TOPDIR)/lib -I$(BUILDDIR)/lib
AM_CXXFLAGS= $(INCFLAGS) @CXXFLAGS@ @QDPXX_CXXFLAGS@
AM_LDFLAGS= -L$(BUILDDIR)/lib @LDFLAGS@ @QDPXX_LDFLAGS@
LDADD= @QDPXX_LIBS@ -lbagel_qdp


bin_PROGRAMS=t_blas

t_blas_SOURCES=t_blas.cc
t_blas_DEPENDENCIES=build_libs

.PHONY: build_libs
build_libs:
	cd $(BUILDDIR)/lib; $(MAKE)


