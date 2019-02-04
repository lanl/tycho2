add_executable(tycho2.mesh.Tetrahedral.test Tetrahedral.test.cpp)
target_link_libraries(tycho2.mesh.Tetrahedral.test PRIVATE tycho2 Catch2::Catch2)
add_test(NAME tycho2.mesh.Tetrahedral.test COMMAND tycho2.mesh.Tetrahedral.test)
