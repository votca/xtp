# CMake generated Testfile for 
# Source directory: /home/yoar/votca/csg-tutorials
# Build directory: /home/yoar/votca/xtp/csg-tutorials
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(regression_UninstallExists "/usr/bin/cmake" "-DFileToCheck=/home/yoar/votca/xtp/csg-tutorials/cmake_uninstall.cmake" "-P" "/home/yoar/votca/csg-tutorials/CMakeModules/FileExists.cmake")
set_tests_properties(regression_UninstallExists PROPERTIES  LABELS "regression")
