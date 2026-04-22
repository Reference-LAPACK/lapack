#ifndef CBLAS_GLOBALS_H
#define CBLAS_GLOBALS_H

#if defined(_WIN32)
  #if defined(CBLAS_DLL_EXPORTS)
    #define CBLAS_GLOBAL_SYMBOL __declspec(dllexport)
  #elif defined(CBLAS_DLL_IMPORTS)
    #define CBLAS_GLOBAL_SYMBOL __declspec(dllimport)
  #else
    #define CBLAS_GLOBAL_SYMBOL
  #endif
#else
  #define CBLAS_GLOBAL_SYMBOL
#endif

extern CBLAS_GLOBAL_SYMBOL int CBLAS_CallFromC;
extern CBLAS_GLOBAL_SYMBOL int RowMajorStrg;

#endif
