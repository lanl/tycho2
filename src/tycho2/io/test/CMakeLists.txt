add_executable(tycho2.json.test json.test.cpp)
target_link_libraries(tycho2.json.test PRIVATE Catch2::Catch2 nlohmann_json::nlohmann_json)
add_test(NAME tycho2.json.test COMMAND tycho2.json.test)

file( GLOB resources "resources/*" )
foreach( resource ${resources})
    file( COPY "${resource}" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" )
endforeach()
