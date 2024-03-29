PROJECT(CFITSIO)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)

# Allow the developer to select if Dynamic or Static libraries are built
OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" ON)
OPTION (USE_PTHREADS "Thread-safe build (using pthreads)" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set (LIB_DESTINATION "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
set (INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/")

# Define project version
SET(${PROJECT_NAME}_MAJOR_VERSION 3)
SET(${PROJECT_NAME}_MINOR_VERSION 36)
SET(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION})
SET(LIB_NAME cfitsio)

# Define
IF(MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
ENDIF()

#add_subdirectory (src)

SET (LIB_TYPE STATIC)
IF (BUILD_SHARED_LIBS)
  SET (LIB_TYPE SHARED)
ENDIF (BUILD_SHARED_LIBS)


FILE(GLOB H_FILES "*.h")

IF (USE_PTHREADS)
    FIND_PACKAGE(pthreads REQUIRED)
    INCLUDE_DIRECTORIES(${PTHREADS_INCLUDE_DIR})
    ADD_DEFINITIONS(-D_REENTRANT)
ENDIF()

SET(SRC_FILES
    buffers.c cfileio.c checksum.c drvrfile.c drvrmem.c
    drvrnet.c drvrsmem.c drvrgsiftp.c editcol.c edithdu.c eval_l.c
    eval_y.c eval_f.c fitscore.c getcol.c getcolb.c getcold.c getcole.c
    getcoli.c getcolj.c getcolk.c getcoll.c getcols.c getcolsb.c
    getcoluk.c getcolui.c getcoluj.c getkey.c group.c grparser.c
    histo.c iraffits.c
    modkey.c putcol.c putcolb.c putcold.c putcole.c putcoli.c
    putcolj.c putcolk.c putcoluk.c putcoll.c putcols.c putcolsb.c
    putcolu.c putcolui.c putcoluj.c putkey.c region.c scalnull.c
    swapproc.c wcssub.c wcsutil.c imcompress.c quantize.c ricecomp.c
    pliocomp.c fits_hcompress.c fits_hdecompress.c zlib/zuncompress.c
    zlib/zcompress.c zlib/adler32.c zlib/crc32.c zlib/inffast.c
    zlib/inftrees.c zlib/trees.c zlib/zutil.c zlib/deflate.c
    zlib/infback.c zlib/inflate.c zlib/uncompr.c simplerng.c
    f77_wrap1.c f77_wrap2.c f77_wrap3.c f77_wrap4.c
)

ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${H_FILES} ${SRC_FILES})
TARGET_LINK_LIBRARIES(${LIB_NAME} m ${PTHREADS_LIBRARY}
)

SET_TARGET_PROPERTIES(${LIB_NAME} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_MAJOR_VERSION})
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION ${LIB_DESTINATION})
install(FILES ${H_FILES} DESTINATION ${INCLUDE_INSTALL_DIR} COMPONENT Devel)

ENABLE_TESTING()

ADD_EXECUTABLE(TestProg testprog.c)
TARGET_LINK_LIBRARIES(TestProg ${LIB_NAME})
ADD_TEST(TestProg TestProg)

ADD_EXECUTABLE(cookbook cookbook.c)
TARGET_LINK_LIBRARIES(cookbook ${LIB_NAME})
ADD_TEST(cookbook cookbook)

ADD_EXECUTABLE(FPack fpack.c fpackutil.c)
TARGET_LINK_LIBRARIES(FPack ${LIB_NAME})

ADD_EXECUTABLE(Funpack funpack.c fpackutil.c)
TARGET_LINK_LIBRARIES(Funpack ${LIB_NAME})

# To expands the command line arguments in Windows, see:
# http://msdn.microsoft.com/en-us/library/8bch7bkk.aspx
if(MSVC)
  set_target_properties(FPack Funpack PROPERTIES 
    LINK_FLAGS "setargv.obj"
  )
endif(MSVC)
