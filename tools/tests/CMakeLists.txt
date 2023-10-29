project(tests)

add_executable(tests tests.cpp main.cpp  ../../quadraticinteger_alvinfractions.cpp ../../quadraticinteger_alvin.cpp ../../rcyclotomic7integer_alvin.cpp ../../quadraticinteger.cpp ../../quadraticinteger_big.cpp ../../alvinfraction.cpp ../../rationalinteger_alvinfractions.cpp ../../rcyclotomic7integer_alvinfractions.cpp ../../alvinfractions.cpp ../../rationalinteger_alvin.cpp ../../rationalinteger.cpp ../../rcyclotomic7integer.cpp ../../rcyclotomic7integer_constants.cpp ../../algebraicinteger.cpp  ../../CoxIter/lib/math_tools.cpp ../../alvin.cpp ../../CoxIter/lib/numbers/number_template.cpp ../../CoxIter/lib/numbers/mpz_rational.cpp ../../CoxIter/lib/string.cpp  ../../CoxIter/lib/regexp.cpp ../../CoxIter/lib/polynomials.cpp  ../../CoxIter/graphs.product.set.cpp ../../CoxIter/graphs.product.cpp ../../CoxIter/graphs.list.n.cpp ../../CoxIter/graphs.list.iterator.cpp ../../CoxIter/graphs.list.cpp ../../CoxIter/graph.cpp ../../CoxIter/coxiter.cpp ../../CoxIter/arithmeticity.cpp)

add_definitions("-std=c++11")
if( CYGWIN )
	add_definitions("-std=gnu++11")
endif()

add_definitions("-Wall")
add_definitions("-Wno-reorder")
add_definitions("-Wno-unknown-pragmas")
add_definitions("-msse3")

set( CMAKE_BUILD_TYPE Release )

#find_package( OpenMP )
#if( OPENMP_FOUND )
#	set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
#	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

#	add_definitions("-fopenmp")
#	add_definitions("-lomp")
#	add_definitions(-D_COMPILE_WITH_OPENMP_)
#endif()

find_library(PCRE_LIBRARY pcre)
target_link_libraries(tests ${PCRE_LIBRARY})

find_library(PARI_LIBRARY pari)
target_link_libraries(tests ${PARI_LIBRARY})

find_library(GMP_LIBRARY gmp)
target_link_libraries(tests ${GMP_LIBRARY})

find_library(GMP_LIBRARYXX gmpxx)
target_link_libraries(tests ${GMP_LIBRARYXX})

find_library(GAOL_LIBRARY gaol)
target_link_libraries(tests ${GAOL_LIBRARY})

find_library(CRLIBM_LIBRARY crlibm)
target_link_libraries(tests ${CRLIBM_LIBRARY})

find_library(GDTOA_LIBRARY gdtoa)
target_link_libraries(tests ${GDTOA_LIBRARY})

install(TARGETS tests RUNTIME DESTINATION bin)
