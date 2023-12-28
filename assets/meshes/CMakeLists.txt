# file(GLOB mesh_files *.obj *.ply)
file(GLOB mesh_files *)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/meshes)
file(COPY ${mesh_files} DESTINATION ${CMAKE_BINARY_DIR}/meshes)
