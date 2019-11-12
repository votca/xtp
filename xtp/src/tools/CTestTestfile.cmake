# CMake generated Testfile for 
# Source directory: /home/yoar/votca/xtp/src/tools
# Build directory: /home/yoar/votca/xtp/xtp/src/tools
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(integration_xtp_mapHelp "xtp_map" "--help")
set_tests_properties(integration_xtp_mapHelp PROPERTIES  LABELS "xtp;tools;csg;ctp;votca;integration")
add_test(integration_xtp_runHelp "xtp_run" "--help")
set_tests_properties(integration_xtp_runHelp PROPERTIES  LABELS "xtp;tools;csg;ctp;votca;integration")
add_test(integration_xtp_toolsHelp "xtp_tools" "--help")
set_tests_properties(integration_xtp_toolsHelp PROPERTIES  LABELS "xtp;tools;csg;ctp;votca;integration")
add_test(integration_xtp_parallelHelp "xtp_parallel" "--help")
set_tests_properties(integration_xtp_parallelHelp PROPERTIES  LABELS "xtp;tools;csg;ctp;votca;integration")
add_test(integration_xtp_dumpHelp "xtp_dump" "--help")
set_tests_properties(integration_xtp_dumpHelp PROPERTIES  LABELS "xtp;tools;csg;ctp;votca;integration")
