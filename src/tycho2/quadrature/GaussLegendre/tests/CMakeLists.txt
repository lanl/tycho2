add_executable(tycho2.Quadrature.test Quadrature.test.cpp)
target_link_libraries(tycho2.Quadrature.test PRIVATE tycho2 Catch2::Catch2 range-v3)
add_test(NAME tycho2.Quadrature.test COMMAND tycho2.Quadrature.test)
  
