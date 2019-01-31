add_executable(tycho2.Quadtrature.test Quadrature.test.cpp)
target_link_libraries(tycho2.Quadtrature.test PRIVATE tycho2 Catch2::Catch2)

add_test(NAME tycho2.Quadtrature.test COMMAND tycho2.Quadtrature.test)
  
