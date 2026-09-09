# main.cpp #includes all the other .cpp files directly (unity build), so it
# is the only file actually passed to the compiler. The others are listed
# and added as sources too (marked HEADER_FILE_ONLY) purely so they show up
# in IDEs and are tracked as dependencies of the target.
set(NLOSINGLEINCLUSIVE_MAIN_SOURCE main.cpp)
set(NLOSINGLEINCLUSIVE_INCLUDED_SOURCES
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

add_executable(nlosingleinclusive
    ${NLOSINGLEINCLUSIVE_MAIN_SOURCE}
    ${NLOSINGLEINCLUSIVE_INCLUDED_SOURCES}
)
set_source_files_properties(
    ${NLOSINGLEINCLUSIVE_INCLUDED_SOURCES}
    PROPERTIES HEADER_FILE_ONLY TRUE
)
