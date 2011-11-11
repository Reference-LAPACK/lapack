##############################################################################
# Copyright (c) 2010, Intel Corp.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#   * Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name of Intel Corporation nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.
##############################################################################
# Contents: Native C interface to LAPACK
# Author: Intel Corporation
# September, 2010
##############################################################################
# This is the make.inc example. The following settings are used:
#
# Compiler: gcc
# Configuration file: turned off (default)
# Complex types: C99 (default)
# Name pattern: mixed case (default)
# (64-bit) Data model: LP64 (default)
#
# Basic include options.
# CC is the C compiler, normally invoked with options CFLAGS.
# LINKER is the linker, invoked with LDFLAGS.
#
# If libraries lapack.a and blas.a are built with
# - ifort, set:    LINKER = ifort
#                  LDFLAGS = -nofor-main
# - gfortran, set: LINKER = gfortran
#
CC = gcc
CFLAGS =
LINKER = $(CC)
LDFLAGS =
#
# The name of the libraries to be created/linked to
# Ensure that the libraries have the same data model (LP64/ILP64).
#
LAPACKE = lapacke.a
LIBS = ../../../lapack-3.2.1/lapack.a ../../../lapack-3.2.1/blas.a -lm
#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib
