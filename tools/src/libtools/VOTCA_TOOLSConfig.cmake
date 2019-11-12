include(CMakeFindDependencyMacro)
find_dependency(Eigen3 NO_MODULE)
find_dependency(Boost 1.53.0 REQUIRED COMPONENTS program_options)
include("${CMAKE_CURRENT_LIST_DIR}/VOTCA_TOOLS_Targets.cmake")
add_executable(VOTCA::votca_compare IMPORTED)
set_property(TARGET VOTCA::votca_compare PROPERTY IMPORTED_LOCATION "/home/yoar/votca/xtp/votca/bin/votca_compare")
if(FALSE)
  if(NOT TARGET MKL::MKL)
    add_library(MKL::MKL UNKNOWN IMPORTED)
    set_target_properties(MKL::MKL
      PROPERTIES
        IMPORTED_LOCATION ""
        INTERFACE_LINK_LIBRARIES ""
        INTERFACE_INCLUDE_DIRECTORIES ""
        IMPORTED_NO_SONAME )
  endif()
  if(NOT TARGET MKL::Interface)
    add_library(MKL::Interface UNKNOWN IMPORTED)
    set_target_properties(MKL::Interface
      PROPERTIES
        IMPORTED_LOCATION "")
  endif()
  if(NOT TARGET MKL::ThreadLayer)
    add_library(MKL::ThreadLayer UNKNOWN IMPORTED)
    set_target_properties(MKL::ThreadLayer
      PROPERTIES
        IMPORTED_LOCATION "")
  endif()
endif()
