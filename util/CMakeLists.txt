set(SRC_MESH
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/SerialMesh.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/ParallelMesh.cc"
   )
set(SRC_MESH2
    ${SRC_MESH}
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/Comm.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/Assert.cc"
   )
set(INC "${CMAKE_CURRENT_SOURCE_DIR}/../src")

add_executable(SerialMeshInfo.x ${SRC_MESH} "SerialMeshInfo.cc")
target_include_directories(SerialMeshInfo.x PUBLIC "${INC}")
target_compile_definitions(SerialMeshInfo.x PUBLIC "-DPARALLEL_MESH_READ_SERIAL_ONLY")
add_executable(ParallelMeshInfo.x ${SRC_MESH2} "ParallelMeshInfo.cc")
target_include_directories(ParallelMeshInfo.x PUBLIC "${INC}")
add_executable(PartitionColumns.x ${SRC_MESH} "PartitionColumns.cc")
target_include_directories(PartitionColumns.x PUBLIC "${INC}")
target_compile_definitions(PartitionColumns.x PUBLIC "-DPARALLEL_MESH_READ_SERIAL_ONLY")
add_executable(RefineSerialMesh.x ${SRC_MESH} "RefineSerialMesh.cc")
target_include_directories(RefineSerialMesh.x PUBLIC "${INC}")
target_compile_definitions(RefineSerialMesh.x PUBLIC "-DPARALLEL_MESH_READ_SERIAL_ONLY")

option(USE_METIS "Use METIS" OFF)
if (USE_METIS)
  set(METIS_DIR "" CACHE PATH "METIS install directory")
  find_library(METIS_LIBRARY NAMES metis PATHS "${METIS_DIR}/lib" NO_DEFAULT_PATH)
  set(METIS_INC "${METIS_DIR}/include")
  add_executable(PartitionMetis.x ${SRC_MESH} "PartitionMetis.cc")
  target_include_directories(PartitionMetis.x PUBLIC "${INC}" "${METIS_INC}")
  target_link_libraries(PartitionMetis.x PUBLIC "${METIS_LIBRARY}")
  target_compile_definitions(PartitionMetis.x PUBLIC "-DPARALLEL_MESH_READ_SERIAL_ONLY")
endif()

option(USE_MOAB "Use MOAB" OFF)
if (USE_MOAB)
  set(MOAB_DIR "" CACHE PATH "MOAB install directory")
  set(CGNS_DIR "" CACHE PATH "CGNS install directory")
  set(HDF5_DIR "" CACHE PATH "HDF5 install directory")
  set(MOAB_INC "${MOAB_DIR}/include")
  find_library(MOAB_LIBRARY NAMES MOAB PATHS "${MOAB_DIR}/lib" NO_DEFAULT_PATH)
  find_library(IMESH_LIBRARY NAMES iMesh PATHS "${MOAB_DIR}/lib" NO_DEFAULT_PATH)
  find_library(CGNS_LIBRARY NAMES cgns PATHS "${CGNS_DIR}/lib" NO_DEFAULT_PATH)
  find_library(HDF5_LIBRARY NAMES hdf5 PATHS "${HDF5_DIR}/lib" NO_DEFAULT_PATH)
  find_package(ZLIB REQUIRED)
  set(MOAB_LIBRARIES
      "${MOAB_LIBRARY}"
      "${IMESH_LIBRARY}"
      "${CGNS_LIBRARY}"
      "${HDF5_LIBRARY}"
      dl
      ZLIB::ZLIB)
  add_executable(MoabToSerialMesh.x ${SRC_MESH} "MoabToSerialMesh.cc")
  target_include_directories(MoabToSerialMesh.x PUBLIC "${INC}" "${MOAB_INC}")
  target_link_libraries(MoabToSerialMesh.x PUBLIC "${MOAB_LIBRARIES}")
  target_compile_definitions(MoabToSerialMesh.x PUBLIC "-DPARALLEL_MESH_READ_SERIAL_ONLY")
  add_executable(SerialMeshToMoab.x ${SRC_MESH} "SerialMeshToMoab.cc")
  target_include_directories(SerialMeshToMoab.x PUBLIC "${INC}" "${MOAB_INC}")
  target_link_libraries(SerialMeshToMoab.x PUBLIC "${MOAB_LIBRARIES}")
  target_compile_definitions(SerialMeshToMoab.x PUBLIC "-DPARALLEL_MESH_READ_SERIAL_ONLY")
  add_executable(ParallelMeshToMoab.x ${SRC_MESH} "ParallelMeshToMoab.cc")
  target_include_directories(ParallelMeshToMoab.x PUBLIC "${INC}" "${MOAB_INC}")
  target_link_libraries(ParallelMeshToMoab.x PUBLIC "${MOAB_LIBRARIES}")
  target_compile_definitions(ParallelMeshToMoab.x PUBLIC "-DPARALLEL_MESH_READ_SERIAL_ONLY")
endif()
