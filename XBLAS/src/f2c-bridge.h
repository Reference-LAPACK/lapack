#if !defined(XBLAS_F2C_BRIDGE_H_)
#define XBLAS_F2C_BRIDGE_H_

/*
  Adapting to a specific Fortran name mangling requires adding
  some of the following CPP flags in make.inc:

    CONFIG_FC_UNDERSCORE       : Add one underscore to each name.

    CONFIG_FC_DBL_UNDERSCORE   : If the name has an underscore in it already,
                                 add two; otherwise add one.

    CONFIG_FC_UCASE            : The name is converted to upper-case; default
                                 is lower-case.

    CONFIG_FC_RETURNS_DBL_REAL : REAL values are returned in C doubles as in
                                 f2c.  (currently unused)

  The first two flags are mutually exclusive.
*/

#if defined(CONFIG_FC_UNDERSCORE)
#define FC_UNDER(x) x##_
#define FC_UNDER2(x) x##_
#elif defined(CONFIG_FC_DBL_UNDERSCORE)
#define FC_UNDER(x) x##_
#define FC_UNDER2(x) x##__
#else
#define FC_UNDER(x) x
#define FC_UNDER2(x) x
#endif

#if defined(CONFIG_FC_UCASE)
#define FC_FUNC(x,X) FC_UNDER(X)
#define FC_FUNC_(x,X) FC_UNDER2(X)
#else
#define FC_FUNC(x,X) FC_UNDER(x)
#define FC_FUNC_(x,X) FC_UNDER2(x)
#endif

#if defined(CONFIG_FC_RETURNS_DBL_REAL)
typedef double fc_real_result_t;
#else
typedef float fc_real_result_t;
#endif

/*
  Possible future config options:

  CONFIG_FC_STUPID_TYPES
    When the Fortran compiler defines 64-bit REALs and INTEGERs
    as the default KIND.
*/

#endif /* XBLAS_F2C_BRIDGE_H_ */
