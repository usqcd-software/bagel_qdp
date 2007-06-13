#ifndef BAGEL_QDP_OPTIONS_H
#define BAGEL_QDP_OPTIONS_H

/* Undef the unwanted from the environment -- eg the compiler command line */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#include "bagel_qdp_options_internal.h"

/* Prefix everything with BAGEL_QDP_ */
static const char* const BAGEL_QDP_PACKAGE(PACKAGE);
static const char* const BAGEL_QDP_PACKAGE_BUGREPORT(PACKAGE_BUGREPORT);
static const char* const BAGEL_QDP_PACKAGE_NAME(PACKAGE_NAME);
static const char* const BAGEL_QDP_PACKAGE_STRING(PACKAGE_STRING);
static const char* const BAGEL_QDP_PACKAGE_TARNAME(PACKAGE_TARNAME);
static const char* const BAGEL_QDP_PACKAGE_VERSION(PACKAGE_VERSION);
static const char* const BAGEL_QDP_VERSION(VERSION);
                                                                                
/* Undef the unwanted */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

/* User defined */
#ifdef USE_DOUBLE
typedef double BAGELQDPFloat;
#else
typedef float BAGELQDPFloat;
#endif

#endif

