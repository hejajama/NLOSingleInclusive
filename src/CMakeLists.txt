# Each .cpp below is its own translation unit (no more #include-ing .cpp
# files into main.cpp) and has a matching .hpp declaring what it exposes.
set(NLOSINGLEINCLUSIVE_SOURCES
    main.cpp
    params.cpp
    utils.cpp
    bksol_nlodisfit.cpp
    Sr_interp_2D.cpp
    Sr_interp_1D.cpp
    func.cpp
    xi_interp.cpp
    xi_int.cpp
    init_interp.cpp
    intde1.cpp
    sigma_NLO.cpp
    sigma_LO.cpp
    init.cpp
)
set(NLOSINGLEINCLUSIVE_HEADERS
    common.hpp
    params.hpp
    utils.hpp
    bksol_nlodisfit.hpp
    Sr_interp_2D.hpp
    Sr_interp_1D.hpp
    func.hpp
    xi_interp.hpp
    xi_int.hpp
    init_interp.hpp
    intde1.hpp
    sigma_NLO.hpp
    sigma_LO.hpp
    init.hpp
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
