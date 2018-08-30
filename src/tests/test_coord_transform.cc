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


BOOST_AUTO_TEST_SUITE(test_coord_transform)

using namespace votca::xtp;
using namespace votca::tools;

BOOST_AUTO_TEST_CASE(incorrect_transform_selection){

    CoordinateTransform C2I(CARTESIAN, INTERNAL);

    // xyz files of different systems

    std::ofstream ammoniaXYZ("ammonia.xyz");

    ammoniaXYZ << " 4" << std::endl;
    ammoniaXYZ << " ammonia" << std::endl;

    ammoniaXYZ << " H      -0.813712     -0.019797      0.024059" << std::endl;
    ammoniaXYZ << " N       0.000059      0.450035      0.427708" << std::endl;
    ammoniaXYZ << " H       0.813749     -0.019930      0.024104" << std::endl;
    ammoniaXYZ << " H      -0.000096      1.389693      0.024103" << std::endl;

    ammoniaXYZ.close();

    // Read in ammonia

    Orbitals ammonia;

    ammonia.LoadFromXYZ("ammonia.xyz");

    CartesianCoords CC(ammonia.QMAtoms());

    InternalCoords ammoniaIC(ammonia.QMAtoms());

    BOOST_CHECK_EQUAL(ammoniaIC.getPossibleNumMols(), 1);


    std::ofstream waterAmmoniaXYZ("water_ammonia.xyz");

    waterAmmoniaXYZ << " 7" << std::endl;
    waterAmmoniaXYZ << " water and ammonia" << std::endl;

    waterAmmoniaXYZ << " N         -7.63720        1.47740       -0.42514"
                    << std::endl;
    waterAmmoniaXYZ << " H         -6.60208        1.56651       -0.42220"
                    << std::endl;
    waterAmmoniaXYZ << " H         -7.89991        0.83024        0.28238"
                    << std::endl;
    waterAmmoniaXYZ << " H         -7.89740        1.24524       -1.37010"
                    << std::endl;
    waterAmmoniaXYZ << " O         -6.00031        3.14833       -0.46360"
                    << std::endl;
    waterAmmoniaXYZ << " H         -5.06057        3.14460       -0.49868"
                    << std::endl;
    waterAmmoniaXYZ << " H         -6.39546        3.97859       -0.79116"
                    << std::endl;

    Orbitals waterAmmonia;
    waterAmmonia.LoadFromXYZ("water_ammonia.xyz");

    InternalCoords WA(waterAmmonia);

    BOOST_CHECK_EQUAL(WA.getPossibleNumMols(), 2);

BOOST_AUTO_TEST_SUITE_END()}
