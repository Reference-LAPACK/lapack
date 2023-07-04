# THE CBLAS C INTERFACE TO BLAS

## Contents
[1. Introduction](#1-introduction)

[1.1 Naming Schemes](#11-naming-schemes)

[1.2 Integers](#12-integers)

[2. Function List](#2-function-list)

[2.1 BLAS Level 1](#21-blas-level-1)

[2.2 BLAS Level 2](#22-blas-level-2)

[2.3 BLAS Level 3](#23-blas-level-3)

[3. Examples](#3-examples)

[3.1 Calling DGEMV](#31-calling-dgemv)

[3.2 Calling DGEMV_64](#32-calling-dgemv_64)

## 1. Introduction
This document describes CBLAS, the C language interface to the Basic Linear Algebra Subprograms (BLAS).
In comparison to BLAS Fortran interfaces CBLAS interfaces support both row-major and column-major matrix
ordering with the `layout` parameter.
The prototypes for CBLAS interfaces, associated macros and type definitions are contained in the header
file [cblas.h](../CBLAS/include/cblas.h)

### 1.1 Naming Schemes
The naming scheme for the CBLAS interface is to take the Fortran BLAS routine name, make it lower case,
and add the prefix `cblas_`. For example, the BLAS routine `DGEMM` becomes `cblas_dgemm`.

CBLAS routines also support `_64` suffix that enables large data arrays support in the LP64 interface library
(default build configuration). This suffix allows mixing LP64 and ILP64 programming models in one application.
For example, `cblas_dgemm` with 32-bit integer type support can be mixed with `cblas_dgemm_64`
that supports 64-bit integer type. 

### 1.2 Integers
Variables with the Fortran type integer are converted to `CBLAS_INT` in CBLAS. By default
the CBLAS interface is built with 32-bit integer type, but it can be re-defined to 64-bit integer type.

## 2. Function List
This section contains the list of the currently available CBLAS interfaces.

### 2.1 BLAS Level 1
* Single Precision Real:
  ```
  SROTG SROTMG SROT SROTM SSWAP SSCAL
  SCOPY SAXPY SDOT  SDSDOT SNRM2 SASUM 
  ISAMAX
  ```
* Double Precision Real:
  ```
  DROTG DROTMG DROT DROTM DSWAP DSCAL
  DCOPY DAXPY DDOT  DSDOT DNRM2 DASUM 
  IDAMAX
  ```
* Single Precision Complex:
  ```
  CROTG CSROT CSWAP CSCAL CSSCAL CCOPY 
  CAXPY CDOTU_SUB CDOTC_SUB ICAMAX SCABS1
  ```
* Double Precision Complex:
  ```
  ZROTG ZDROT ZSWAP ZSCAL ZDSCAL ZCOPY 
  ZAXPY ZDOTU_SUB ZDOTC_SUB IZAMAX DCABS1
  DZNRM2 DZASUM
  ```
### 2.2 BLAS Level 2
* Single Precision Real:
  ```
  SGEMV SGBMV SGER  SSBMV SSPMV SSPR
  SSPR2 SSYMV SSYR  SSYR2 STBMV STBSV
  STPMV STPSV STRMV STRSV
  ```
* Double Precision Real:
  ```
  DGEMV DGBMV DGER  DSBMV DSPMV DSPR
  DSPR2 DSYMV DSYR  DSYR2 DTBMV DTBSV
  DTPMV DTPSV DTRMV DTRSV
  ```
* Single Precision Complex:
  ```
  CGEMV CGBMV CHEMV CHBMV CHPMV CTRMV
  CTBMV CTPMV CTRSV CTBSV CTPSV CGERU
  CGERC CHER  CHER2 CHPR  CHPR2
  ```
* Double Precision Complex:
  ```
  ZGEMV ZGBMV ZHEMV ZHBMV ZHPMV ZTRMV
  ZTBMV ZTPMV ZTRSV ZTBSV ZTPSV ZGERU
  ZGERC ZHER  ZHER2 ZHPR  ZHPR2
  ```
### 2.3 BLAS Level 3
* Single Precision Real:
  ```
  SGEMM SSYMM SSYRK SSERK2K STRMM STRSM
  ```
* Double Precision Real:
  ```
  DGEMM DSYMM DSYRK DSERK2K DTRMM DTRSM
  ```
* Single Precision Complex:
  ```
  CGEMM CSYMM CHEMM CHERK CHER2K CTRMM
  CTRSM CSYRK CSYR2K
  ```
* Double Precision Complex:
  ```
  ZGEMM ZSYMM ZHEMM ZHERK ZHER2K ZTRMM 
  ZTRSM ZSYRK ZSYR2K
  ```

## 3. Examples
This section contains examples of calling CBLAS functions from a C program.

### 3.1 Calling DGEMV
The variable declarations should be as follows:
```
   double *a, *x, *y;
   double alpha, beta;
   CBLAS_INT m, n, lda, incx, incy;
```
The CBLAS function call is then:
```
cblas_dgemv( CblasColMajor, CblasNoTrans, m, n, alpha, a, lda, x, incx, beta,
                y, incy );
```

### 3.2 Calling DGEMV_64
The variable declarations should be as follows:
```
   double *a, *x, *y;
   double alpha, beta;
   int64_t m, n, lda, incx, incy;
```
The CBLAS function call is then:
```
cblas_dgemv_64( CblasColMajor, CblasNoTrans, m, n, alpha, a, lda, x, incx, beta,
                y, incy );
```
