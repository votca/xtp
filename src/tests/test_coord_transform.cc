#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE coord_transform
#include <cassert>
#include <iostream>
#include <vector>
#include <fstream>
#include <boost/test/unit_test.hpp>
#include <votca/xtp/qmatom.h>
#include <votca/tools/constants.h>
#include <votca/xtp/internalcoords.h>
#include <votca/xtp/coordtransform.h>
#include <votca/xtp/orbitals.h>
#include <cmath>


BOOST_AUTO_TEST_SUITE(test_coord_transform)

using namespace votca::xtp;
using namespace votca::tools;

BOOST_AUTO_TEST_CASE(four_bonds){
    std::ofstream ammoniaXYZ("ammonia.xyz");

    ammoniaXYZ << " 4" << std::endl;
    ammoniaXYZ << " ammonia" << std::endl;

    ammoniaXYZ << " H      -0.813712     -0.019797      0.024059" << std::endl;
    ammoniaXYZ << " N       0.000059      0.450035      0.427708" << std::endl;
    ammoniaXYZ << " H       0.813749     -0.019930      0.024104" << std::endl;
    ammoniaXYZ << " H      -0.000096      1.389693      0.024103" << std::endl;

    ammoniaXYZ.close();

    Orbitals ammonia;

    ammonia.LoadFromXYZ("ammonia.xyz");

    InternalCoords ammoniaIC(ammonia);

    BOOST_CHECK_EQUAL(ammoniaIC.getNumBonds(), 3);
}
BOOST_AUTO_TEST_CASE(single_angle){
    std::ofstream single_angleXYZ("single_angle.xyz");

    single_angleXYZ << "3" << std::endl;
    single_angleXYZ << " " << std::endl;
    single_angleXYZ << "O          0.02272        2.28696        0.00000" << std::endl;
    single_angleXYZ << "H          0.99272        2.28696        0.00000" << std::endl;
    single_angleXYZ << "H         -0.30061        2.60542       -0.85729" << std::endl;

    Orbitals single_angle;
    single_angle.LoadFromXYZ("single_angle.xyz");
    InternalCoords ic(single_angle);

    BOOST_CHECK_EQUAL(ic.getNumAngles(), 1);
    double tol = 1e-6;

    BOOST_CHECK(std::abs(ic.Vector()[4] - acos(cos(109.5)) < tol));
}

BOOST_AUTO_TEST_CASE(single_dihedral_two_angles_ic_test){
    std::ofstream single_dihedralXYZ ("single_dihedral.xyz");

    single_dihedralXYZ << "4" << std::endl;
    single_dihedralXYZ << "" << std::endl;
    single_dihedralXYZ << "C         -6.56591        1.43600        0.00000" << std::endl;
    single_dihedralXYZ << "C         -5.49591        1.43600        0.00000" << std::endl;
    single_dihedralXYZ << "H         -6.92258        0.44220       -0.17333" << std::endl;
    single_dihedralXYZ << "H         -5.13925        2.42980        0.17333" << std::endl;
    single_dihedralXYZ.close();
    Orbitals single_dihedral;
    single_dihedral.LoadFromXYZ("single_dihedral.xyz");
    InternalCoords ic(single_dihedral, true);

    // std::cout << ic.getNumBonds() << std::endl;
    // std::cout << ic.getNumAuxBonds() << std::endl;

    BOOST_CHECK_EQUAL(ic.getNumDihedrals(), 1);
    BOOST_CHECK_EQUAL(ic.getNumAngles(), 2);

}

BOOST_AUTO_TEST_CASE(one_angles){
    std::ofstream co2XYZ("co2.xyz");

    co2XYZ << "3" << std::endl;
    co2XYZ << " " << std::endl;
    co2XYZ << "O         0.00000        0.00000        0.00000" << std::endl;
    co2XYZ << "C         1.42000        0.00000        0.00000" << std::endl;
    co2XYZ << "O         2.01000        0.00000        0.00000" << std::endl;


    Orbitals co2;
    co2.LoadFromXYZ("co2.xyz");

    InternalCoords co2IC(co2);

    BOOST_CHECK_EQUAL(co2IC.getNumBonds(), 2);
    BOOST_CHECK_EQUAL(co2IC.getNumAngles(), 1);
}

BOOST_AUTO_TEST_CASE(no_dihedrals){
    std::ofstream ethyneXYZ("ethyne.xyz");

    ethyneXYZ << "" << std::endl;
    ethyneXYZ << "" << std::endl;
    ethyneXYZ << "C         -3.93046        0.00000        0.00000" << std::endl;
    ethyneXYZ << "C         -2.14320        0.00000        0.00000" << std::endl;
    ethyneXYZ << "H         -4.92442        0.00000        0.00000" << std::endl;
    ethyneXYZ << "H         -1.14923        0.00000        0.00000" << std::endl;

    Orbitals ethyne; ethyne.LoadFromXYZ("ethyne.xyz");


    InternalCoords ethIC(ethyne, false);

    BOOST_CHECK_EQUAL(ethIC.getNumDihedrals(), 0);
    BOOST_CHECK_EQUAL(ethIC.getNumAngles(), 2);
    BOOST_CHECK_EQUAL(ethIC.getNumBonds(), 3);

}

BOOST_AUTO_TEST_SUITE_END()
