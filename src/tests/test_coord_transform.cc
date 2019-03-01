#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE coord_transform
#include <boost/test/unit_test.hpp>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <votca/tools/constants.h>
#include <votca/xtp/coordtransform.h>
#include <votca/xtp/internalcoords.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmatom.h>

BOOST_AUTO_TEST_SUITE(test_coord_transform)

using namespace votca::xtp;
using namespace votca::tools;

const double tol = 1e-6;

BOOST_AUTO_TEST_CASE(single_bond) {
  std::cout << "Single bond" << std::endl;
  std::ofstream single_bondXYZ("single_bond.xyz");

  single_bondXYZ << "2" << std::endl;
  single_bondXYZ << " " << std::endl;
  single_bondXYZ << "C          0.00000        0.00000        0.00000"
                 << std::endl;
  single_bondXYZ << "C          1.00000        0.00000        0.00000"
                 << std::endl;

  Orbitals single_bond;
  single_bond.LoadFromXYZ("single_bond.xyz");
  InternalCoords ic(single_bond);

  BOOST_CHECK_EQUAL(ic.getNumBonds(), 1);
  BOOST_CHECK_EQUAL(ic.getNumAtoms(), 2);
}

BOOST_AUTO_TEST_CASE(single_angle) {

  std::cout << "Single Angle - H2O" << std::endl;
  std::ofstream single_angleXYZ("single_angle.xyz");

  single_angleXYZ << "3" << std::endl;
  single_angleXYZ << " " << std::endl;
  single_angleXYZ << "H          0.00000        0.00000        0.00000"
                  << std::endl;
  single_angleXYZ << "O          1.00000        0.00000        0.00000"
                  << std::endl;
  single_angleXYZ << "H          1.00000        1.00000        0.00000"
                  << std::endl;

  Orbitals single_angle;
  single_angle.LoadFromXYZ("single_angle.xyz");
  InternalCoords ic(single_angle);
  BOOST_CHECK_EQUAL(ic.getNumAtoms(), 3);
  BOOST_CHECK_EQUAL(ic.getNumAngles(), 1);
  BOOST_CHECK_EQUAL(ic.getNumBonds(), 2);
  BOOST_CHECK_EQUAL(ic.getNumDihedrals(), 0);
  double tol = 1e-6;

  std::cout << ic << std::endl;
}

BOOST_AUTO_TEST_CASE(single_dihedral_two_angles_ic_test) {
  std::cout << "single dihedral" << std::endl;
  std::ofstream single_dihedralXYZ("single_dihedral.xyz");

  single_dihedralXYZ << "4" << std::endl;
  single_dihedralXYZ << "" << std::endl;
  single_dihedralXYZ << "H          0.00000        0.00000        0.00000"
                     << std::endl;
  single_dihedralXYZ << "C          1.00000        0.00000        0.00000"
                     << std::endl;
  single_dihedralXYZ << "C          1.00000        1.00000        0.00000"
                     << std::endl;
  single_dihedralXYZ << "H          1.00000        1.00000        1.00000"
                     << std::endl;
  single_dihedralXYZ.close();
  Orbitals single_dihedral;
  single_dihedral.LoadFromXYZ("single_dihedral.xyz");
  InternalCoords ic(single_dihedral, true);

  BOOST_CHECK_EQUAL(ic.getNumBonds(), 3);
  BOOST_CHECK_EQUAL(ic.getNumAngles(), 2);
  BOOST_CHECK_EQUAL(ic.getNumDihedrals(), 1);
  std::cout << ic << std::endl;
}

BOOST_AUTO_TEST_CASE(linear_molecule) {
  std::cout << "linear_molecule1 co2" << std::endl;
  std::ofstream co2XYZ("co2.xyz");

  co2XYZ << "3" << std::endl;
  co2XYZ << " " << std::endl;
  co2XYZ << "O         0.00000        0.00000        0.00000" << std::endl;
  co2XYZ << "C         1.42000        0.00000        0.00000" << std::endl;
  co2XYZ << "O         2.84000        0.00000        0.00000" << std::endl;

  Orbitals co2;
  co2.LoadFromXYZ("co2.xyz");

  InternalCoords ic(co2);

  BOOST_CHECK_EQUAL(ic.getNumBonds(), 2);
  BOOST_CHECK_EQUAL(ic.getNumAngles(), 1);
  BOOST_CHECK_EQUAL(ic.getNumDihedrals(), 0);
  std::cout << ic << std::endl;
}

BOOST_AUTO_TEST_CASE(no_dihedrals) {
  std::cout << "no dihedrals - ethyne " << std::endl;
  std::ofstream ethyneXYZ("ethyne.xyz");

  ethyneXYZ << "4" << std::endl;
  ethyneXYZ << "" << std::endl;
  ethyneXYZ << "H         0.000000        0.00000        0.00000" << std::endl;
  ethyneXYZ << "C         1.000000        0.00000        0.00000" << std::endl;
  ethyneXYZ << "C         2.000000        0.00000        0.00000" << std::endl;
  ethyneXYZ << "H         3.000000        0.00000        0.00000" << std::endl;

  Orbitals ethyne;
  ethyne.LoadFromXYZ("ethyne.xyz");

  InternalCoords ic(ethyne);

  BOOST_CHECK_EQUAL(ic.getNumDihedrals(), 0);
  BOOST_CHECK_EQUAL(ic.getNumAngles(), 2);
  BOOST_CHECK_EQUAL(ic.getNumBonds(), 3);

  std::cout << ic << std::endl;
}

BOOST_AUTO_TEST_CASE(ammonia_internal_coords) {
  std::cout << "ammonia" << std::endl;

  std::ofstream ammoniaXYZ("ammonia.xyz");

  ammoniaXYZ << " 4" << std::endl;
  ammoniaXYZ << " ammonia" << std::endl;

  ammoniaXYZ << " N       0.000059      0.450035      0.427708" << std::endl;
  ammoniaXYZ << " H      -0.813712     -0.019797      0.024059" << std::endl;
  ammoniaXYZ << " H       0.813749     -0.019930      0.024104" << std::endl;
  ammoniaXYZ << " H      -0.000096      1.389693      0.024103" << std::endl;

  ammoniaXYZ.close();

  Orbitals ammonia;

  ammonia.LoadFromXYZ("ammonia.xyz");

  InternalCoords ammoniaIC(ammonia);

  BOOST_CHECK_EQUAL(ammoniaIC.getNumBonds(), 3);
  BOOST_CHECK_EQUAL(ammoniaIC.getNumAngles(), 3);
  BOOST_CHECK_EQUAL(ammoniaIC.getNumDihedrals(), 6);
  std::cout << ammoniaIC << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
