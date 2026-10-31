find_package(Doxygen REQUIRED)

set(DOXYGEN_PROJECT_BRIEF "LAPACK: Linear Algebra PACKage")
set(DOXYGEN_PROJECT_NUMBER ${LAPACK_VERSION})
set(DOXYGEN_OUTPUT_DIRECTORY "${LAPACK_BINARY_DIR}/DOCS")
set(DOXYGEN_PROJECT_LOGO "${CMAKE_CURRENT_SOURCE_DIR}/lapack.png")
set(DOXYGEN_OPTIMIZE_FOR_FORTRAN YES)
set(DOXYGEN_SOURCE_BROWSER YES)
set(DOXYGEN_CREATE_SUBDIRS YES)
set(DOXYGEN_SEPARATE_MEMBER_PAGES YES)
set(DOXYGEN_TAB_SIZE 8)
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_FILE_PATTERNS *.f *.f90 *.c *.h )
set(DOXYGEN_RECURSIVE YES)
set(DOXYGEN_GENERATE_TREEVIEW YES)
set(DOXYGEN_DOT_IMAGE_FORMAT svg)
set(DOXYGEN_INTERACTIVE_SVG YES)
set(DOXYGEN_WARN_NO_PARAMDOC YES)
set(DOXYGEN_WARN_LOGFILE doxygen_error)
set(DOXYGEN_LAYOUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/DoxygenLayout.xml")

# Exclude functions that are duplicated, creating conflicts.
set(DOXYGEN_EXCLUDE .git
                    .github
                    SRC/VARIANTS
                    BLAS/SRC/lsame.f
                    BLAS/SRC/xerbla.f
                    BLAS/SRC/xerbla_array.f
                    INSTALL/slamchf77.f
                    INSTALL/dlamchf77.f )

if (BUILD_HTML_DOCUMENTATION)
    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_GENERATE_MAN NO)
    set(DOXYGEN_INLINE_SOURCES YES)
    set(DOXYGEN_CALL_GRAPH YES)
    set(DOXYGEN_CALLER_GRAPH YES)

    set(DOXYGEN_HTML_OUTPUT explore-html)
    set(DOXYGEN_HTML_TIMESTAMP YES)
    doxygen_add_docs(
        html

        # Doxygen INPUT =
        BLAS
        CBLAS
        SRC
        INSTALL
        TESTING
        DOCS/groups-usr.dox
        README.md

        WORKING_DIRECTORY "${LAPACK_SOURCE_DIR}"
        COMMENT "Generating html LAPACK documentation (it will take some time... time to grab a coffee)"
    )
    unset(DOXYGEN_HTML_OUTPUT)
    unset(DOXYGEN_HTML_TIMESTAMP)
endif()
if (BUILD_MAN_DOCUMENTATION)
    set(DOXYGEN_GENERATE_HTML NO)
    set(DOXYGEN_GENERATE_MAN YES)
    set(DOXYGEN_INLINE_SOURCES NO)
    set(DOXYGEN_CALL_GRAPH NO)
    set(DOXYGEN_CALLER_GRAPH NO)

    set(DOXYGEN_MAN_LINKS YES)
    doxygen_add_docs(
        man

        # Doxygen INPUT =
        BLAS
        CBLAS
        SRC
        INSTALL
        TESTING
        DOCS/groups-usr.dox

        WORKING_DIRECTORY "${LAPACK_SOURCE_DIR}"
        COMMENT "Generating man LAPACK documentation"
    )
    unset(DOXYGEN_MAN_LINKS)
endif()
