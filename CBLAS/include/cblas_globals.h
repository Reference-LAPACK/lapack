#ifndef CBLAS_GLOBALS_H
#define CBLAS_GLOBALS_H

#if defined(CBLAS_DLL_IMPORTS)
  #define CBLAS_GLOBAL_SYMBOL __declspec(dllimport)
#else
  #define CBLAS_GLOBAL_SYMBOL
#endif

extern CBLAS_GLOBAL_SYMBOL int CBLAS_CallFromC;
extern CBLAS_GLOBAL_SYMBOL int RowMajorStrg;

#endif
