# CMake generated Testfile for 
# Source directory: /home/yoar/votca/tools/src/tools
# Build directory: /home/yoar/votca/xtp/tools/src/tools
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(integration_votca_propertyHelp "votca_property" "--help")
set_tests_properties(integration_votca_propertyHelp PROPERTIES  LABELS "tools;votca;integration")
