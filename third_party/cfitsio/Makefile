#
# Makefile for cfitsio library:
#       libcfits.a
#
# Oct-96 : original version by 
#
#       JDD/WDP
#       NASA GSFC
#       Oct 1996
#
# 25-Jan-01 : removed conditional drvrsmem.c compilation because this
#             is now handled within the source file itself.
# 09-Mar-98 : modified to conditionally compile drvrsmem.c. Also
# changes to target all (deleted clean), added DEFS, LIBS, added
# DEFS to .c.o, added SOURCES_SHMEM and MY_SHMEM, expanded getcol*
# and putcol* in SOURCES, modified OBJECTS, mv changed to /bin/mv
# (to bypass aliasing), cp changed to /bin/cp, add smem and
# testprog targets. See also changes and comments in configure.in


# CFITSIO version numbers:
CFITSIO_MAJOR	= 3
CFITSIO_MINOR	= 37
CFITSIO_SONAME	= 2

prefix		= /Users/nate/Dropbox/ThreeDAsymmetryTests/ThreeDAsymmetryPackage_V2/third_party/cfitsio
exec_prefix	= ${prefix}
DESTDIR		= 
CFITSIO_LIB	= ${DESTDIR}${exec_prefix}/lib
CFITSIO_INCLUDE	= ${DESTDIR}${prefix}/include
INSTALL_DIRS	= ${prefix} ${CFITSIO_INCLUDE} ${CFITSIO_LIB} ${CFITSIO_LIB}/pkgconfig


SHELL =		/bin/sh
ARCHIVE =	ar rv
RANLIB =	ranlib
CC =		gcc
CFLAGS =	-g -O2 -arch i386 -arch x86_64 -Dg77Fortran -fPIC -fno-common
SSE_FLAGS =	
FC =		gfortran
LDFLAGS =	$(CFLAGS)
DEFS =		-DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_BZLIB_H=1 -DHAVE_BZIP2=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MATH_H=1 -DHAVE_LIMITS_H=1 -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -DHAVE_FTRUNCATE=1 -DHAVE_LONGLONG=1 -DHAVE_SHMEM_SERVICES=1 -DHAVE_UNION_SEMUN=1 -DHAVE_NET_SERVICES=1
LIBS =		
FLEX =		flex
BISON =		bison

SHLIB_LD =	gcc -dynamiclib -arch i386 -arch x86_64 -install_name libcfitsio.${CFITSIO_SONAME}${SHLIB_SUFFIX} -compatibility_version ${CFITSIO_SONAME} -current_version ${CFITSIO_SONAME}.${CFITSIO_MAJOR}.${CFITSIO_MINOR}
SHLIB_SUFFIX =	.dylib
CFITSIO_SHLIB =	libcfitsio.${CFITSIO_SONAME}.${CFITSIO_MAJOR}.${CFITSIO_MINOR}${SHLIB_SUFFIX}
CFITSIO_SHLIB_SONAME = libcfitsio.${CFITSIO_SONAME}${SHLIB_SUFFIX}


CORE_SOURCES = 	buffers.c cfileio.c checksum.c drvrfile.c drvrmem.c \
		drvrnet.c drvrsmem.c drvrgsiftp.c editcol.c edithdu.c eval_l.c \
		eval_y.c eval_f.c fitscore.c getcol.c getcolb.c getcold.c getcole.c \
		getcoli.c getcolj.c getcolk.c getcoll.c getcols.c getcolsb.c \
		getcoluk.c getcolui.c getcoluj.c getkey.c group.c grparser.c \
		histo.c iraffits.c \
		modkey.c putcol.c putcolb.c putcold.c putcole.c putcoli.c \
		putcolj.c putcolk.c putcoluk.c putcoll.c putcols.c putcolsb.c \
		putcolu.c putcolui.c putcoluj.c putkey.c region.c scalnull.c \
		swapproc.c wcssub.c wcsutil.c imcompress.c quantize.c ricecomp.c \
		pliocomp.c fits_hcompress.c fits_hdecompress.c \
		simplerng.c

ZLIB_SOURCES =	zlib/adler32.c zlib/crc32.c zlib/deflate.c zlib/infback.c \
		zlib/inffast.c zlib/inflate.c zlib/inftrees.c zlib/trees.c \
		zlib/uncompr.c zlib/zcompress.c zlib/zuncompress.c zlib/zutil.c

SOURCES = ${CORE_SOURCES} ${ZLIB_SOURCES} ${FITSIO_SRC}

OBJECTS = 	${SOURCES:.c=.o}

CORE_OBJECTS = 	${CORE_SOURCES:.c=.o} ${ZLIB_SOURCES:.c=.o}


FITSIO_SRC =	f77_wrap1.c f77_wrap2.c f77_wrap3.c f77_wrap4.c

# ============ description of all targets =============
#       -  <<-- ignore error code

all:
		@if [ "x${FC}" = x ]; then \
			${MAKE} all-nofitsio; \
		else \
			${MAKE} stand_alone; \
		fi

all-nofitsio:
		${MAKE} stand_alone "FITSIO_SRC="

stand_alone:	libcfitsio.a

libcfitsio.a:	${OBJECTS}
		${ARCHIVE} libcfitsio.a ${OBJECTS}; \
		${RANLIB} libcfitsio.a;

shared: libcfitsio${SHLIB_SUFFIX}

libcfitsio${SHLIB_SUFFIX}: ${OBJECTS}
		${SHLIB_LD} ${LDFLAGS} -o ${CFITSIO_SHLIB} ${OBJECTS} -lm ${LIBS}
		@if [ "x${CFITSIO_SHLIB_SONAME}" != x ]; then \
			ln -sf ${CFITSIO_SHLIB} ${CFITSIO_SHLIB_SONAME}; \
			ln -sf ${CFITSIO_SHLIB_SONAME} $@; \
		fi

install:	libcfitsio.a $(INSTALL_DIRS)
		/bin/cp -a libcfitsio* ${CFITSIO_LIB}
		/bin/cp fitsio.h fitsio2.h longnam.h drvrsmem.h ${CFITSIO_INCLUDE}
		/bin/cp cfitsio.pc ${CFITSIO_LIB}/pkgconfig

.c.o:
		$(CC) -c -o ${<D}/${@F} $(CFLAGS) $(DEFS) $<

swapproc.o:	swapproc.c
		$(CC) -c $(CFLAGS) $(SSE_FLAGS) $(DEFS) $<

smem:		smem.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o smem smem.o -L. -lcfitsio -lm

testprog:	testprog.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o testprog testprog.o -L. -lcfitsio -lm ${LIBS}

fpack:		fpack.o fpackutil.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o fpack fpack.o fpackutil.o libcfitsio.a -lm ${LIBS}

funpack:	funpack.o fpackutil.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o funpack funpack.o fpackutil.o libcfitsio.a -lm ${LIBS}

fitscopy:	fitscopy.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o fitscopy fitscopy.o -L. -lcfitsio -lm ${LIBS}

speed:		speed.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o speed speed.o -L. -lcfitsio -lm ${LIBS}

imcopy:		imcopy.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o imcopy imcopy.o -L. -lcfitsio -lm ${LIBS}

listhead:	listhead.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o listhead listhead.o -L. -lcfitsio -lm ${LIBS}

cookbook:	cookbook.o libcfitsio.a ${OBJECTS}
		${CC} $(CFLAGS) $(DEFS) -o cookbook cookbook.o -L. -lcfitsio -lm ${LIBS}

eval:		# Rebuild eval_* files from flex/bison source
		$(FLEX) -t eval.l > eval_l.c1
		# Note workaround for yyfree=fffree conflict
		/bin/sed -e 's/yy/ff/g' -e 's/YY/FF/g' eval_l.c1 -e 's/fffree/yyfffree/g' > eval_l.c
		/bin/rm -f eval_l.c1
		$(BISON) -d -v -y eval.y
		/bin/sed -e 's/yy/ff/g' -e 's/YY/FF/g' y.tab.c > eval_y.c
		/bin/sed -e 's/yy/ff/g' -e 's/YY/FF/g' y.tab.h > eval_tab.h
		/bin/rm -f y.tab.c y.tab.h

clean:
	-	/bin/rm -f *.o zlib/*.o libcfitsio* fpack funpack \
			smem testprog y.output so_locations

distclean:	clean
	-	/bin/rm -f Makefile cfitsio.pc config.log config.status configure.lineno

# Make target which outputs the list of the .o contained in the cfitsio lib
# usefull to build a single big shared library containing Tcl/Tk and other
# extensions.  used for the Tcl Plugin. 

cfitsioLibObjs:
	@echo ${CORE_OBJECTS}

cfitsioLibSrcs:
	@echo ${SOURCES}

# This target actually builds the objects needed for the lib in the above
# case
objs: ${CORE_OBJECTS}

$(INSTALL_DIRS):
	@if [ ! -d $@ ]; then mkdir -p $@; fi
