# CMake generated Testfile for 
# Source directory: /home/yoar/votca/csg/share/template
# Build directory: /home/yoar/votca/xtp/csg/share/template
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(template_serialHelp "template_serial" "--help")
set_tests_properties(template_serialHelp PROPERTIES  LABELS "csg;tools;votca")
add_test(template_threadedHelp "template_threaded" "--help")
set_tests_properties(template_threadedHelp PROPERTIES  LABELS "csg;tools;votca")
