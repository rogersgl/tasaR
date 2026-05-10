#ifndef PANDASEQ_CONFIG_H
#define PANDASEQ_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* Package metadata expected by args.c/help/version output */
#define PACKAGE_NAME "pandaseq"
#define PACKAGE_TARNAME "pandaseq"
#define PACKAGE_VERSION "2.11"
#define PACKAGE_STRING "pandaseq 2.11"
#define PACKAGE_BUGREPORT "andre@masella.name"
#define PACKAGE_URL "https://github.com/neufeld/pandaseq"

/* Build-time values from upstream configure.ac */
#define MAX_LEN 450
#define VERSION_MAJOR 2
#define VERSION_MINOR 11
#define LIB_MAJOR 7
#define LIB_MINOR 0

/* Feature flags used by the vendored source on macOS */
#define HAVE_PTHREAD 1
#define HAVE_SYS_PARAM_H 1
#define HAVE_SYS_SYSCTL_H 1
#define HAVE_BZLIB_H 1
#define HAVE_LTDL_H 1
#define HAVE_LIBBZ2 1
#define HAVE_LIBLTDL 1
#define HAVE_ZLIB_H 1
#define HAVE_LIBZ 1

#endif