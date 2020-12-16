/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/* Defines debugging . */
/* #undef DEBUG */

/* No debugging. */
#define DEBUGLEVEL 0

/* Define to 1 to support Advanced Bit Manipulation */
#define HAVE_ABM 1

/* Define to 1 to support Multi-Precision Add-Carry Instruction Extensions */
#define HAVE_ADX 1

/* Define to 1 to support Advanced Encryption Standard New Instruction Set
   (AES-NI) */
#define HAVE_AES 1

/* Support Altivec instructions */
/* #undef HAVE_ALTIVEC */

/* Define to 1 to support Advanced Vector Extensions */
#define HAVE_AVX 1

/* Define to 1 to support Advanced Vector Extensions 2 */
#define HAVE_AVX2 1

/* Define to 1 to support AVX-512 Byte and Word Instructions */
#define HAVE_AVX512_BW 1

/* Define to 1 to support AVX-512 Conflict Detection Instructions */
#define HAVE_AVX512_CD 1

/* Define to 1 to support AVX-512 Doubleword and Quadword Instructions */
#define HAVE_AVX512_DQ 1

/* Define to 1 to support AVX-512 Exponential & Reciprocal Instructions */
/* #undef HAVE_AVX512_ER */

/* Define to 1 to support AVX-512 Foundation Extensions */
#define HAVE_AVX512_F 1

/* Define to 1 to support AVX-512 Integer Fused Multiply Add Instructions */
/* #undef HAVE_AVX512_IFMA */

/* Define to 1 to support AVX-512 Conflict Prefetch Instructions */
/* #undef HAVE_AVX512_PF */

/* Define to 1 to support AVX-512 Vector Byte Manipulation Instructions */
/* #undef HAVE_AVX512_VBMI */

/* Define to 1 to support AVX-512 Vector Length Extensions */
#define HAVE_AVX512_VL 1

/* Define to 1 to support Bit Manipulation Instruction Set 1 */
#define HAVE_BMI1 1

/* Define to 1 to support Bit Manipulation Instruction Set 2 */
#define HAVE_BMI2 1

/* Define to 1 if you have the `clock_gettime' function. */
#define HAVE_CLOCK_GETTIME 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <float.h> header file. */
#define HAVE_FLOAT_H 1

/* Define to 1 to support Fused Multiply-Add Extensions 3 */
#define HAVE_FMA3 1

/* Define to 1 to support Fused Multiply-Add Extensions 4 */
/* #undef HAVE_FMA4 */

/* Define to 1 if you have the `gethrtime' function. */
/* #undef HAVE_GETHRTIME */

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `pthread' library (-lpthread). */
#define HAVE_LIBPTHREAD 1

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if you have the `localtime_r' function. */
#define HAVE_LOCALTIME_R 1

/* Define to 1 if you have the <mach/mach.h> header file. */
/* #undef HAVE_MACH_MACH_H */

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define to 1 to support Multimedia Extensions */
#define HAVE_MMX 1

/* Define to 1 to support Memory Protection Extensions */
#define HAVE_MPX 1

/* Defined if OpenMP should and can be used */
#define HAVE_OPENMP 1

/* Define to 1 if you have the `pow' function. */
#define HAVE_POW 1

/* Define to 1 to support Prefetch Vector Data Into Caches WT1 */
/* #undef HAVE_PREFETCHWT1 */

/* Define to 1 to support Digital Random Number Generator */
#define HAVE_RDRND 1

/* Define to 1 if your system has a GNU libc compatible `realloc' function,
   and to 0 otherwise. */
#define HAVE_REALLOC 1

/* Define to 1 if you have the `realpath' function. */
#define HAVE_REALPATH 1

/* Define to 1 to support Secure Hash Algorithm Extension */
/* #undef HAVE_SHA */

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 to support Streaming SIMD Extensions */
#define HAVE_SSE 1

/* Define to 1 to support Streaming SIMD Extensions */
#define HAVE_SSE2 1

/* Define to 1 to support Streaming SIMD Extensions 3 */
#define HAVE_SSE3 1

/* Define to 1 to support Streaming SIMD Extensions 4.1 */
#define HAVE_SSE4_1 1

/* Define to 1 to support Streaming SIMD Extensions 4.2 */
#define HAVE_SSE4_2 1

/* Define to 1 to support AMD Streaming SIMD Extensions 4a */
/* #undef HAVE_SSE4a */

/* Define to 1 to support Supplemental Streaming SIMD Extensions 3 */
#define HAVE_SSSE3 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strcspn' function. */
#define HAVE_STRCSPN 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strstr' function. */
#define HAVE_STRSTR 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Support VSX instructions */
/* #undef HAVE_VSX */

/* Define to 1 to support eXtended Operations Extensions */
/* #undef HAVE_XOP */

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Name of package */
#define PACKAGE "kalign"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME "kalign"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "kalign 3.3"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "kalign"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "3.3"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Enable extensions on AIX 3, Interix.  */
#ifndef _ALL_SOURCE
# define _ALL_SOURCE 1
#endif
/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif
/* Enable threading extensions on Solaris.  */
#ifndef _POSIX_PTHREAD_SEMANTICS
# define _POSIX_PTHREAD_SEMANTICS 1
#endif
/* Enable extensions on HP NonStop.  */
#ifndef _TANDEM_SOURCE
# define _TANDEM_SOURCE 1
#endif
/* Enable general extensions on Solaris.  */
#ifndef __EXTENSIONS__
# define __EXTENSIONS__ 1
#endif


/* Version number of package */
#define VERSION "3.3"

/* Define to 1 if on MINIX. */
/* #undef _MINIX */

/* Define to 2 if the system does not provide POSIX.1 features except with
   this defined. */
/* #undef _POSIX_1_SOURCE */

/* Define to 1 if you need to in order for `stat' and other things to work. */
/* #undef _POSIX_SOURCE */

/* Define for Solaris 2.5.1 so the uint32_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT32_T */

/* Define for Solaris 2.5.1 so the uint64_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT64_T */

/* Define for Solaris 2.5.1 so the uint8_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT8_T */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to the type of a signed integer type of width exactly 16 bits if
   such a type exists and the standard includes do not define it. */
/* #undef int16_t */

/* Define to the type of a signed integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
/* #undef int32_t */

/* Define to the type of a signed integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
/* #undef int64_t */

/* Define to the type of a signed integer type of width exactly 8 bits if such
   a type exists and the standard includes do not define it. */
/* #undef int8_t */

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to rpl_realloc if the replacement function should be used. */
/* #undef realloc */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef ssize_t */

/* Define to the type of an unsigned integer type of width exactly 16 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint16_t */

/* Define to the type of an unsigned integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint32_t */

/* Define to the type of an unsigned integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint64_t */

/* Define to the type of an unsigned integer type of width exactly 8 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint8_t */
