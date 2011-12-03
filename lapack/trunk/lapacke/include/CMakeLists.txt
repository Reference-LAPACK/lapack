ADD_CUSTOM_TARGET( include ALL DEPENDS 
    lapacke.h lapacke_config.h lapacke_utils.h )

INSTALL(
  FILES lapacke.h lapacke_config.h
  DESTINATION include )
