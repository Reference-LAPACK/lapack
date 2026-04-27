# This Makefile compiles the CBLAS routines

# Error handling routines for level 2 & 3
set(ERRHAND cblas_globals.c cblas_xerbla.c xerbla.c)

#
#
# CBLAS routines
#
# Level 1
#
#

# Files for level 1 single precision real
set(SLEV1_C
  cblas_srotg.c cblas_srotmg.c cblas_srot.c cblas_srotm.c cblas_sswap.c
  cblas_sscal.c cblas_scopy.c cblas_saxpy.c cblas_sdot.c cblas_sdsdot.c
  cblas_snrm2.c cblas_sasum.c cblas_isamax.c cblas_saxpby.c)

set(SLEV1_F
  sdotsub.f sdsdotsub.f snrm2sub.f sasumsub.f isamaxsub.f)

# Files for level 1 double precision real
set(DLEV1_C
  cblas_drotg.c cblas_drotmg.c cblas_drot.c  cblas_drotm.c cblas_dswap.c
  cblas_dscal.c cblas_dcopy.c cblas_daxpy.c cblas_ddot.c cblas_dsdot.c
  cblas_dnrm2.c cblas_dasum.c cblas_idamax.c cblas_daxpby.c)

set(DLEV1_F
  ddotsub.f dsdotsub.f dnrm2sub.f dasumsub.f idamaxsub.f)

# Files for level 1 single precision complex
set(CLEV1_C
  cblas_crotg.c cblas_csrot.c cblas_cswap.c cblas_cscal.c cblas_csscal.c
  cblas_ccopy.c cblas_caxpy.c cblas_cdotu_sub.c cblas_cdotc_sub.c cblas_scnrm2.c
  cblas_scasum.c cblas_icamax.c cblas_scabs1.c cblas_caxpby.c)

set(CLEV1_F
  cdotcsub.f cdotusub.f icamaxsub.f scasumsub.f scnrm2sub.f scabs1sub.f)

# Files for level 1 double precision complex
set(ZLEV1_C
  cblas_zrotg.c cblas_zdrot.c cblas_zswap.c cblas_zscal.c cblas_zdscal.c
  cblas_zcopy.c cblas_zaxpy.c cblas_zdotu_sub.c cblas_zdotc_sub.c cblas_dznrm2.c
  cblas_dzasum.c cblas_izamax.c cblas_dcabs1.c cblas_zaxpby.c)

set(ZLEV1_F
  zdotcsub.f zdotusub.f izamaxsub.f dzasumsub.f dznrm2sub.f dcabs1sub.f)

#
#
# CBLAS routines
#
# Level 2
#
#

# Files for level 2 single precision real
set(SLEV2
  cblas_sgemv.c cblas_sgbmv.c cblas_sger.c cblas_ssbmv.c cblas_sspmv.c
  cblas_sspr.c cblas_sspr2.c cblas_ssymv.c cblas_ssyr.c cblas_ssyr2.c
  cblas_stbmv.c cblas_stbsv.c cblas_stpmv.c cblas_stpsv.c cblas_strmv.c
  cblas_strsv.c cblas_sskewsymv.c cblas_sskewsyr2.c)

# Files for level 2 double precision real
set(DLEV2
  cblas_dgemv.c cblas_dgbmv.c cblas_dger.c cblas_dsbmv.c cblas_dspmv.c
  cblas_dspr.c cblas_dspr2.c cblas_dsymv.c cblas_dsyr.c cblas_dsyr2.c
  cblas_dtbmv.c cblas_dtbsv.c cblas_dtpmv.c cblas_dtpsv.c cblas_dtrmv.c
  cblas_dtrsv.c cblas_dskewsymv.c cblas_dskewsyr2.c)

# Files for level 2 single precision complex
set(CLEV2
  cblas_cgemv.c cblas_cgbmv.c cblas_chemv.c cblas_chbmv.c cblas_chpmv.c
  cblas_ctrmv.c cblas_ctbmv.c cblas_ctpmv.c cblas_ctrsv.c cblas_ctbsv.c
  cblas_ctpsv.c cblas_cgeru.c cblas_cgerc.c cblas_cher.c cblas_cher2.c
  cblas_chpr.c cblas_chpr2.c)

# Files for level 2 double precision complex
set(ZLEV2
  cblas_zgemv.c cblas_zgbmv.c cblas_zhemv.c cblas_zhbmv.c cblas_zhpmv.c
  cblas_ztrmv.c cblas_ztbmv.c cblas_ztpmv.c cblas_ztrsv.c cblas_ztbsv.c
  cblas_ztpsv.c cblas_zgeru.c cblas_zgerc.c cblas_zher.c cblas_zher2.c
  cblas_zhpr.c cblas_zhpr2.c)

#
#
# CBLAS routines
#
# Level 3
#
#

# Files for level 3 single precision real
set(SLEV3
  cblas_sgemm.c cblas_ssymm.c cblas_ssyrk.c cblas_ssyr2k.c cblas_strmm.c
  cblas_strsm.c cblas_sgemmtr.c cblas_sskewsymm.c cblas_sskewsyr2k.c)

# Files for level 3 double precision real
set(DLEV3
  cblas_dgemm.c cblas_dsymm.c cblas_dsyrk.c cblas_dsyr2k.c cblas_dtrmm.c
  cblas_dtrsm.c cblas_dgemmtr.c cblas_dskewsymm.c cblas_dskewsyr2k.c)

# Files for level 3 single precision complex
set(CLEV3
  cblas_cgemm.c cblas_csymm.c cblas_chemm.c cblas_cherk.c cblas_cher2k.c
  cblas_ctrmm.c cblas_ctrsm.c cblas_csyrk.c cblas_csyr2k.c cblas_cgemmtr.c)

# Files for level 3 double precision complex
set(ZLEV3
  cblas_zgemm.c cblas_zsymm.c cblas_zhemm.c cblas_zherk.c cblas_zher2k.c
  cblas_ztrmm.c cblas_ztrsm.c cblas_zsyrk.c cblas_zsyr2k.c cblas_zgemmtr.c)

set(SOURCES_C)
set(SOURCES_F)
if(BUILD_SINGLE)
  list(APPEND SOURCES_C ${SLEV1_C} ${SLEV2} ${SLEV3} ${ERRHAND})
  list(APPEND SOURCES_F ${SLEV1_F})
endif()
if(BUILD_DOUBLE)
  list(APPEND SOURCES_C ${DLEV1_C} ${DLEV2} ${DLEV3} ${ERRHAND})
  list(APPEND SOURCES_F ${DLEV1_F})
endif()
if(BUILD_COMPLEX)
  list(APPEND SOURCES_C ${CLEV1_C} ${CLEV2} ${CLEV3} ${ERRHAND})
  list(APPEND SOURCES_F ${CLEV1_F})
endif()
if(BUILD_COMPLEX16)
  list(APPEND SOURCES_C ${ZLEV1_C} ${ZLEV2} ${ZLEV3} ${ERRHAND})
  list(APPEND SOURCES_F ${ZLEV1_F})
endif()
list(REMOVE_DUPLICATES SOURCES_C)
list(REMOVE_DUPLICATES SOURCES_F)

add_library(${CBLASLIB}_obj OBJECT ${SOURCES_C} ${SOURCES_F})
if(HAS_ATTRIBUTE_WEAK_SUPPORT)
  target_compile_definitions(${CBLASLIB}_obj PRIVATE
    "$<$<COMPILE_LANGUAGE:C>:HAS_ATTRIBUTE_WEAK_SUPPORT>")
endif()

if(BUILD_INDEX64_EXT_API)
  include(ExtendedAPIHelpers)

  # Define list of C files for 64bit interface
  set(SOURCES_64_C ${SOURCES_C})
  list(REMOVE_ITEM SOURCES_64_C cblas_globals.c)

  # Generate list of Fortran files for 64bit interface
  generate_64bit_suffixed_sources(${CBLASLIB} SOURCES_F SOURCES_64_F)

  add_library(${CBLASLIB}_64_obj OBJECT ${SOURCES_64_C} ${SOURCES_64_F})
  target_compile_definitions(${CBLASLIB}_64_obj PRIVATE
    "$<$<COMPILE_LANGUAGE:C>:WeirdNEC>" "$<$<COMPILE_LANGUAGE:C>:CBLAS_API64>")
  target_compile_options(${CBLASLIB}_64_obj PRIVATE
    "$<$<COMPILE_LANGUAGE:Fortran>:${FOPT_ILP64}>")
  if(HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(${CBLASLIB}_64_obj PRIVATE
      "$<$<COMPILE_LANGUAGE:C>:HAS_ATTRIBUTE_WEAK_SUPPORT>")
  endif()
endif()

add_library(${CBLASLIB}
  $<TARGET_OBJECTS:${CBLASLIB}_obj>
  $<$<BOOL:${BUILD_INDEX64_EXT_API}>: $<TARGET_OBJECTS:${CBLASLIB}_64_obj>>)

set_target_properties(
  ${CBLASLIB} PROPERTIES
  LINKER_LANGUAGE C
  VERSION ${LAPACK_VERSION}
  SOVERSION ${LAPACK_MAJOR_VERSION}
  )

target_include_directories(${CBLASLIB} PUBLIC
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(${CBLASLIB} PUBLIC ${BLAS_LIBRARIES})
lapack_install_library(${CBLASLIB})
