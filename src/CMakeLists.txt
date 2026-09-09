# Each .cpp below is its own translation unit (no more #include-ing .cpp
# files into main.cpp) and has a matching .hpp declaring what it exposes.
set(NLOSINGLEINCLUSIVE_SOURCES
    main.cpp
    params.cpp
    utils.cpp
    spline_wrappers.cpp
    dipole_amplitude.cpp
    Sr_interp_1D.cpp
    nlo_coefficients.cpp
    point_tables.cpp
    intde1.cpp
    sigma_NLO.cpp
    sigma_LO.cpp
)
set(NLOSINGLEINCLUSIVE_HEADERS
    common.hpp
    params.hpp
    utils.hpp
    spline_wrappers.hpp
    dipole_amplitude.hpp
    Sr_interp_1D.hpp
    nlo_coefficients.hpp
    point_tables.hpp
    intde1.hpp
    sigma_NLO.hpp
    sigma_LO.hpp
)

add_executable(nlosingleinclusive
    ${NLOSINGLEINCLUSIVE_SOURCES}
    ${NLOSINGLEINCLUSIVE_HEADERS}
)
# Headers carry no build step of their own; list them only so they show up
# and are tracked as dependencies in IDEs.
set_source_files_properties(
    ${NLOSINGLEINCLUSIVE_HEADERS}
    PROPERTIES HEADER_FILE_ONLY TRUE
)
