# CMake generated Testfile for 
# Source directory: /home/yoar/votca/csgapps
# Build directory: /home/yoar/votca/xtp/csgapps
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(fluctuationsHelp "fluctuations/csg_fluctuations" "--help")
set_tests_properties(fluctuationsHelp PROPERTIES  LABELS "csgapps;csg;tools;votca")
add_test(orientcorrHelp "orientcorr/csg_orientcorr" "--help")
set_tests_properties(orientcorrHelp PROPERTIES  LABELS "csgapps;csg;tools;votca")
add_test(part_distHelp "part_dist/csg_part_dist" "--help")
set_tests_properties(part_distHelp PROPERTIES  LABELS "csgapps;csg;tools;votca")
add_test(partial_rdfHelp "partial_rdf/csg_partial_rdf" "--help")
set_tests_properties(partial_rdfHelp PROPERTIES  LABELS "csgapps;csg;tools;votca")
add_test(radiiHelp "radii/csg_radii" "--help")
set_tests_properties(radiiHelp PROPERTIES  LABELS "csgapps;csg;tools;votca")
add_test(sphericalorderHelp "sphericalorder/csg_sphericalorder" "--help")
set_tests_properties(sphericalorderHelp PROPERTIES  LABELS "csgapps;csg;tools;votca")
add_test(traj_forceHelp "traj_force/csg_traj_force" "--help")
set_tests_properties(traj_forceHelp PROPERTIES  LABELS "csgapps;csg;tools;votca")
