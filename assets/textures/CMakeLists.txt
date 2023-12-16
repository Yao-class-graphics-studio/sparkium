file(GLOB texture_files *.jpg *.png *.bmp *.hdr *.txt)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/textures)
file(COPY ${texture_files} DESTINATION ${CMAKE_BINARY_DIR}/textures)
