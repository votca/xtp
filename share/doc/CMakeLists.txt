if(NOT VOTCA_SPHINX_DIR)
  set(VOTCA_SPHINX_DIR "${CMAKE_CURRENT_BINARY_DIR}")
endif()
find_package(Python REQUIRED COMPONENTS Interpreter)
file(GLOB XMLS ${CMAKE_CURRENT_SOURCE_DIR}/../xtp/xml/*.xml)
set(XML_RST_FILES)
message(status "VOTCA_SPHINX_DIR: ${VOTCA_SPHINX_DIR}")
foreach(_XML ${XMLS})
  get_filename_component(_FILE "${_XML}" NAME_WE)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/${_FILE}.rst COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/extract_metadata.py -i ${_XML} -o ${VOTCA_SPHINX_DIR}/${_FILE}.rst
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/extract_metadata.py ${_XML})
  list(APPEND XML_RST_FILES ${VOTCA_SPHINX_DIR}/${_FILE}.rst)
endforeach()
add_custom_target(doc-xtp DEPENDS ${XML_RST_FILES})