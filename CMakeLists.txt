cmake_minimum_required(VERSION 3.1)

project(alvin)

add_definitions("-std=c++11")
if( CYGWIN )
	add_definitions("-std=gnu++11")
endif()

add_definitions("-Wall")
add_definitions("-Wno-reorder")
add_definitions("-Wno-unknown-pragmas")
add_definitions("-msse3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

set( CMAKE_BUILD_TYPE Release )

add_executable(alvin rationalinteger_invariantsqf.cpp quadraticinteger_alvin.cpp rationalinteger_alvin.cpp app.cpp CoxIter/lib/string.cpp alvin.cpp main.cpp CoxIter/lib/regexp.cpp )
find_library(GAOL_LIBRARY gaol)
if( GAOL_LIBRARY )
	target_sources(alvin PUBLIC "${CMAKE_CURRENT_LIST_DIR}/rcyclotomic7integer_alvin.cpp")
endif()

find_library(PARI_LIBRARY pari)
target_link_libraries(alvin ${PARI_LIBRARY})

# ------------------------------------------------------------------
# CoxIter

# CoxIter: Graphs: graphs, products of graphs and iterators
add_library( coxiter_graphs STATIC CoxIter/graphs.product.set.cpp CoxIter/graphs.product.cpp CoxIter/graphs.list.n.cpp CoxIter/graphs.list.iterator.cpp CoxIter/graphs.list.cpp CoxIter/graph.cpp )

# CoxIter
add_library( coxiter STATIC CoxIter/coxiter.cpp)
target_link_libraries( coxiter coxiter_graphs )

add_library( coxiter_growthrate STATIC CoxIter/lib/paripolynomials.cpp CoxIter/growthrate.cpp )
target_link_libraries( coxiter_growthrate ${PARI_LIBRARY} )
target_link_libraries( coxiter coxiter_growthrate )

target_link_libraries( alvin coxiter )

# ------------------------------------------------------------------
# alvin

# Maths: polynomials, fractions, varia functions
add_library( alvin_maths STATIC CoxIter/lib/math_tools.cpp CoxIter/lib/polynomials.cpp CoxIter/lib/numbers/number_template.cpp CoxIter/lib/numbers/mpz_rational.cpp )
target_link_libraries( alvin alvin_maths )

# Numbers: rationals, quadratic and RC7
add_library( alvin_numbers STATIC CoxIter/lib/math_tools.cpp quadraticinteger.cpp rationalinteger.cpp algebraicinteger.cpp )
if( GAOL_LIBRARY )
	target_sources( alvin_numbers PUBLIC "${CMAKE_CURRENT_LIST_DIR}/rcyclotomic7integer.cpp" )
endif()
target_link_libraries( alvin alvin_numbers )

# Fractions
add_library( alvin_fractions STATIC alvinfraction.cpp alvinfractions.cpp rationalinteger_alvinfractions.cpp quadraticinteger_alvinfractions.cpp )
if( GAOL_LIBRARY )
	target_sources( alvin_fractions PUBLIC "${CMAKE_CURRENT_LIST_DIR}/rcyclotomic7integer_alvinfractions.cpp" )
endif()
target_link_libraries( alvin_fractions alvin_numbers )
target_link_libraries( alvin alvin_fractions )

if( GAOL_LIBRARY )
	add_library( RC7alvinConstants STATIC rcyclotomic7integer_constants.cpp )
	target_link_libraries( alvin RC7alvinConstants )
endif()

add_library( alvin_nr rationalinteger_notreflective.cpp notreflective.cpp infinitensymetries.cpp rationalinteger_infinitensymetries.cpp quadraticinteger_infinitensymetries.cpp quadraticinteger_big.cpp)
if( GAOL_LIBRARY )
	target_sources( alvin_nr PUBLIC "${CMAKE_CURRENT_LIST_DIR}/rcyclotomic7integer_infinitensymetries.cpp" )
endif()

target_link_libraries( alvin alvin_nr )

# ------------------------------------------------------------------
# External libraries
find_package( OpenMP )
if( OPENMP_FOUND )
	set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

	add_definitions("-fopenmp")
	add_definitions("-lomp")
	add_definitions(-D_COMPILE_WITH_OPENMP_)
endif()

find_library(PCRE_LIBRARY pcre)
target_link_libraries(alvin ${PCRE_LIBRARY})

find_library(IGRAPH_LIBRARY igraph)
target_link_libraries(alvin ${IGRAPH_LIBRARY})

find_library(GMP_LIBRARY gmp)
target_link_libraries(alvin ${GMP_LIBRARY})

find_library(GMP_LIBRARYXX gmpxx)
target_link_libraries(alvin ${GMP_LIBRARYXX})

if( GAOL_LIBRARY )
	target_link_libraries(alvin ${GAOL_LIBRARY})

	find_library(CRLIBM_LIBRARY crlibm)
	target_link_libraries(alvin ${CRLIBM_LIBRARY})

	find_library(GDTOA_LIBRARY gdtoa)
	target_link_libraries(alvin ${GDTOA_LIBRARY})

	message("RC7 is there")

	add_definitions(-D_RC7AVAILABLE_)
else()
	message("No RC7")
endif()

# ------------------------------------------------------------------
# If the dot program is present, we will be able to create images
find_program( PROGRAM_DOT dot )
if( PROGRAM_DOT )
	add_definitions(-D_DOT_PROGRAM_FOUND_)
endif()

install(TARGETS alvin RUNTIME DESTINATION bin)
