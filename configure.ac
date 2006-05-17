AC_INIT([bagel_qdp],[1.1.5],[bjoo@jlab.org])
AC_CONFIG_AUX_DIR(config)

dnl boilerplace
AC_LANG([C++])
AC_PROG_CXX([g++])
AC_PROG_RANLIB
AC_CHECK_TOOL(AR, ar, [ar])

dnl init automake
AM_INIT_AUTOMAKE

dnl set up assembler for automake
AM_PROG_AS

dnl user options go to wfm_options_internal.h
AC_CONFIG_HEADERS([lib/bagel_qdp_options_internal.h])

dnl user must specify where the bagel installation is 
AC_ARG_WITH(bagel,
   AC_HELP_STRING([--with_bagel=BAGELDIR],
                  [ Use BAGEL installed in BAGELDIR ] ),
   [ bagel_bindir=${with_bagel}/bin
     bagel_CXX=`${bagel_bindir}/bagel-config --cxx`
     bagel_CXXFLAGS=`${bagel_bindir}/bagel-config --cxxflags`
     bagel_LDFLAGS=`${bagel_bindir}/bagel-config --ldflags`
     bagel_LIBS=`${bagel_bindir}/bagel-config --libs`
     AC_SUBST(bagel_CXX)
     AC_SUBST(bagel_CXXFLAGS)
     AC_SUBST(bagel_LDFLAGS)
     AC_SUBST(bagel_LIBS)
    ],
)

AC_ARG_WITH(qdp,
AC_HELP_STRING(
    [--with-qdp=DIR],
    [Build tests with QDP++ where QDP++ is installed in DIR. Not specifying this will disable tests which you should do if you want to install a production version]
  ),
  [
   QDPXX_HOME="$with_qdp"
   QDP_TEST_ENABLED="yes"
  ],
  [ 
    AC_MSG_NOTICE( [ QDP test disabled ]) 
    QDP_TEST_ENABLED="no"
  ]
)


AC_MSG_NOTICE( [Bagel CXX:      ${bagel_CXX} ] )
AC_MSG_NOTICE( [Bagel CXXFLAGS: ${bagel_CXXFLAGS} ] )
AC_MSG_NOTICE( [Bagel LDFLAGS:  ${bagel_LDFLAGS} ] )
AC_MSG_NOTICE( [Bagel LIBS:     ${bagel_LIBS} ] )

dnl allow setting of precision
AC_ARG_ENABLE(precision, 
   AC_HELP_STRING([--enable-precision=PRECISION], 
                  [ Enable Precision (single|double). Default is double] ),
   [ ac_PRECISION=${enable_precision}],
   [ ac_PRECISION="double" ]
)

dnl set the target cpu
AC_ARG_ENABLE(target_cpu,
   AC_HELP_STRING([--enable-target-cpu=CPU],
       [ Create Assembler for target CPU (ppc440|powerIII|alpha064|alpha164|alpha264|usparcII|noarch) default is ppc440 ]),
       [ ac_TARGET_CPU=${enable_target_cpu} ],	
       [ ac_TARGET_CPU=ppc440 ]
)

dnl configuration
case ${ac_PRECISION} in
  single)
	AC_MSG_NOTICE([Using single precision])
	#  don't define this if not using double precision
	#AC_DEFINE([USE_DOUBLE],[0],[ Use Single Precision] )
        ;;
  double)
	AC_MSG_NOTICE([Using double precision])
	AC_DEFINE([USE_DOUBLE],[1],[ Use Double Precision] )
	;;
  *)
	AC_MSG_ERROR([Invalid value for precision])
	;;
esac

dnl get munged target name
case ${ac_TARGET_CPU} in
     ppc440|powerIII|alpha064|alpha164|alpha264|noarch|usparcII)
	AC_MSG_NOTICE([Producing code for ${ac_TARGET_CPU} ])
	if test "X${ac_PRECISION}X" == "XsingleX";
        then
           bagel_MUNGED_TARGET=${ac_TARGET_CPU}"s";
        else
	   bagel_MUNGED_TARGET=${ac_TARGET_CPU};
        fi
	AC_MSG_NOTICE([Bagel Munged Target Name: ${bagel_MUNGED_TARGET}])
	;;
      *)
	AC_MSG_ERROR([Unknown BAGEL Target CPU ${ac_TARGET_CPU}]);
	;;
esac
AC_SUBST(bagel_MUNGED_TARGET)

dnl decide whether the MADD, NOMADD, or NOARCH target list gets
dnl added to the library list
AM_CONDITIONAL(BUILD_ARCH_SOURCES,
	[ test  "X${ac_TARGET_CPU}X" == "Xppc440X" \
	     -o "X${ac_TARGET_CPU}X" == "XpowerIIIX" ])

AM_CONDITIONAL(BUILD_NOARCH_SOURCES,
	[ test  "X${ac_TARGET_CPU}X" == "XnoarchX" ])


if test "X${QDP_TEST_ENABLED}X" = "XyesX"; then 

  if test "X${QDPXX_HOME}X" = "XX" ; then	   
    AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [])
  else
    AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [], [${QDPXX_HOME}/bin:${PATH}])
  fi

  if test "X${QDPXX_CONFIG}X" = "XX" ; then
    AC_MSG_ERROR([QDP++ configuration program qdp++-config not found.])
  fi

  AC_MSG_NOTICE([Found QDP++ configuration program ${QDPXX_CONFIG}])
  AC_SUBST(QDPXX_CXX, "`${QDPXX_CONFIG} --cxx`")
  AC_MSG_NOTICE([QDP++ compiler: ${CXX}])
  AC_SUBST(QDPXX_CXXFLAGS, "`${QDPXX_CONFIG} --cxxflags`")
  AC_MSG_NOTICE([QDP++ compile flags: ${QDPXX_CXXFLAGS}])
  AC_SUBST(QDPXX_LDFLAGS,  "`${QDPXX_CONFIG} --ldflags`")
  AC_MSG_NOTICE([QDP++ linking flags: ${QDPXX_LDFLAGS}])
  AC_SUBST(QDPXX_LIBS,     "`${QDPXX_CONFIG} --libs`")
  AC_MSG_NOTICE([QDP++ libraries flags: ${QDPXX_LIBS}])
  AC_SUBST(QDPXX_ND,       "`${QDPXX_CONFIG} --Nd`")
  if test "X${QDPXX_ND}X" = "XX" ; then
    AC_MSG_ERROR([Cannot determine QDP++ spacetime dimension])
  else
    AC_MSG_NOTICE([QDP++ spacetime dimension: ${QDPXX_ND}])
  fi
  AC_SUBST(QDPXX_NC,       "`${QDPXX_CONFIG} --Nc`")
  AC_MSG_NOTICE([QDP++ number of colors: ${QDPXX_NC}])
  AC_SUBST(QDPXX_NS,       "`${QDPXX_CONFIG} --Ns`")
  if test "X${QDPXX_NS}X" = "XX" ; then
    AC_MSG_ERROR([Cannot determine number of spin components in QDP++])
  else
    AC_MSG_NOTICE([QDP++ number of spin components: ${QDPXX_NS}])
  fi
  AC_SUBST(QDPXX_PARALLEL_ARCH, "`${QDPXX_CONFIG} --parallel-arch`")
  AC_MSG_NOTICE([QDP++ parallel arch: ${QDPXX_PARALLEL_ARCH}])

  AC_SUBST(QDPXX_PRECISION, "`${QDPXX_CONFIG} --precision`")
  AC_MSG_NOTICE([QDP++ precision: ${QDPXX_PRECISION}])
  dnl This function is defined in acinclude.m4
  dnl it tries to try to compile the program
  AC_MSG_CHECKING([if we can compile/link a simple QDP++ program])
  PAC_QDPXX_LINK_CXX_FUNC(
    ${QDPXX_CXX},
    ${QDPXX_CXXFLAGS},
    ${QDPXX_LDFLAGS},
    ${QDPXX_LIBS},
    ,
    ,
    [qdpxx_link_ok=yes],
    [qdpxx_link_ok=no]
  )
  if test "X${qdpxx_link_ok}X" = "XyesX" ; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Cannot compile/link a program with QDP++.
      Use --with-qdp++=<dir> to select a working version.])
  fi
else 
     dnl QDP is disabled
     AC_SUBST(QDPXX_CXX, "")
     AC_SUBST(QDPXX_CXXFLAGS,"")
     AC_SUBST(QDPXX_LDFLAGS,"")
     AC_SUBST(QDPXX_LIBS,"")
fi

AM_CONDITIONAL(BUILD_QDP_TEST, [ test "X${QDP_TEST_ENABLED}X" = "XyesX" ])

AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(code_gen/Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(tests/Makefile)
AC_CONFIG_FILES(bagel-qdp-config)
AC_OUTPUT
