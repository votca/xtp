/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE gwbse_exact_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/bse.h>
#include <votca/xtp/sigma.h>
#include <votca/xtp/convergenceacc.h>
#include <votca/xtp/gwbse_exact.h>
#include <votca/xtp/sigma_spectral.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(gwbse_exact_test)

BOOST_AUTO_TEST_CASE(bse_hamiltonian) {
    
    // ****** Write Test Files ******

    ofstream moleculefile;
    ofstream dftbasisfile;
    ofstream auxbasisfile;
    
    // <editor-fold defaultstate="collapsed" desc="methane">
    moleculefile = ofstream("methane.xyz");
    moleculefile << " 5"                                            << endl;
    moleculefile << " methane"                                      << endl;
    moleculefile << " C            .000000     .000000     .000000" << endl;
    moleculefile << " H            .629118     .629118     .629118" << endl;
    moleculefile << " H           -.629118    -.629118     .629118" << endl;
    moleculefile << " H            .629118    -.629118    -.629118" << endl;
    moleculefile << " H           -.629118     .629118    -.629118" << endl;
    moleculefile.close();
    // </editor-fold>

    // <editor-fold defaultstate="collapsed" desc="co">
    moleculefile = ofstream("co.xyz");
    moleculefile << " 2"                                            << endl;
    moleculefile << " CO"                                           << endl;
    moleculefile << " C            0.000000    .000000     .000000" << endl;
    moleculefile << " O            1.000000    .000000     .000000" << endl;
    moleculefile.close();
    // </editor-fold>

    // <editor-fold defaultstate="collapsed" desc="3-21G">
    dftbasisfile = ofstream("3-21G.xml");
    dftbasisfile << "<basis name=\"3-21G\">" << endl;
    dftbasisfile << "  <element name=\"H\">" << endl;
    dftbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
    dftbasisfile << "      <constant decay=\"5.447178e+00\">" << endl;
    dftbasisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>" << endl;
    dftbasisfile << "      </constant>" << endl;
    dftbasisfile << "      <constant decay=\"8.245470e-01\">" << endl;
    dftbasisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>" << endl;
    dftbasisfile << "      </constant>" << endl;
    dftbasisfile << "    </shell>" << endl;
    dftbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
    dftbasisfile << "      <constant decay=\"1.831920e-01\">" << endl;
    dftbasisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>" << endl;
    dftbasisfile << "      </constant>" << endl;
    dftbasisfile << "    </shell>" << endl;
    dftbasisfile << "  </element>" << endl;
    dftbasisfile << "  <element name=\"C\">" << endl;
    dftbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
    dftbasisfile << "      <constant decay=\"1.722560e+02\">" << endl;
    dftbasisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>" << endl;
    dftbasisfile << "      </constant>" << endl;
    dftbasisfile << "      <constant decay=\"2.591090e+01\">" << endl;
    dftbasisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>" << endl;
    dftbasisfile << "      </constant>" << endl;
    dftbasisfile << "      <constant decay=\"5.533350e+00\">" << endl;
    dftbasisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>" << endl;
    dftbasisfile << "      </constant>" << endl;
    dftbasisfile << "    </shell>" << endl;
    dftbasisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
    dftbasisfile << "      <constant decay=\"3.664980e+00\">" << endl;
    dftbasisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>" << endl;
    dftbasisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>" << endl;
    dftbasisfile << "      </constant>" << endl;
    dftbasisfile << "      <constant decay=\"7.705450e-01\">" << endl;
    dftbasisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>" << endl;
    dftbasisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>" << endl;
    dftbasisfile << "      </constant>" << endl;
    dftbasisfile << "    </shell>" << endl;
    dftbasisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
    dftbasisfile << "      <constant decay=\"1.958570e-01\">" << endl;
    dftbasisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>" << endl;
    dftbasisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>" << endl;
    dftbasisfile << "      </constant>" << endl;
    dftbasisfile << "    </shell>" << endl;
    dftbasisfile << "  </element>" << endl;
    dftbasisfile << "</basis>" << endl;
    dftbasisfile.close();
    // </editor-fold>

    // <editor-fold defaultstate="collapsed" desc="def2-SVP">
    dftbasisfile = ofstream("def2-SVP.xml");
    dftbasisfile << "<basis name=\"def2-SVP\">" << std::endl;
    dftbasisfile << "  <!--Basis set created by xtp_basisset from def2-svp.nw at Mon Aug  6 15:17:57 2018-->" << std::endl;
    dftbasisfile << "  <element name=\"H\">" << std::endl;
    dftbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
    dftbasisfile << "      <constant decay=\"1.301070e+01\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"1.968216e-02\" type=\"S\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"1.962257e+00\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"1.379652e-01\" type=\"S\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"4.445380e-01\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"4.783193e-01\" type=\"S\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"1.219496e-01\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"0.000000e+00\" type=\"S\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "    </shell>" << std::endl;
    dftbasisfile << "    <shell scale=\"1.0\" type=\"P\">" << std::endl;
    dftbasisfile << "      <constant decay=\"8.000000e-01\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "    </shell>" << std::endl;
    dftbasisfile << "  </element>" << std::endl;
    dftbasisfile << "  <element name=\"C\">" << std::endl;
    dftbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
    dftbasisfile << "      <constant decay=\"1.238402e+03\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"5.456883e-03\" type=\"S\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"1.862900e+02\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"4.063841e-02\" type=\"S\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"4.225118e+01\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"1.802559e-01\" type=\"S\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"1.167656e+01\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"4.631512e-01\" type=\"S\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"3.593051e+00\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"4.408717e-01\" type=\"S\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"4.024515e-01\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"0.000000e+00\" type=\"S\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"1.309018e-01\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"0.000000e+00\" type=\"S\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "    </shell>" << std::endl;
    dftbasisfile << "    <shell scale=\"1.0\" type=\"P\">" << std::endl;
    dftbasisfile << "      <constant decay=\"9.468097e+00\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"3.838787e-02\" type=\"P\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"2.010355e+00\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"2.111703e-01\" type=\"P\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"5.477100e-01\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"5.132817e-01\" type=\"P\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "      <constant decay=\"1.526861e-01\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"0.000000e+00\" type=\"P\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "    </shell>" << std::endl;
    dftbasisfile << "    <shell scale=\"1.0\" type=\"D\">" << std::endl;
    dftbasisfile << "      <constant decay=\"8.000000e-01\">" << std::endl;
    dftbasisfile << "        <contractions factor=\"1.000000e+00\" type=\"D\"/>" << std::endl;
    dftbasisfile << "      </constant>" << std::endl;
    dftbasisfile << "    </shell>" << std::endl;
    dftbasisfile << "  </element>" << std::endl;
    dftbasisfile << "</basis>" << std::endl;
    dftbasisfile.close();
    // </editor-fold>

    // <editor-fold defaultstate="collapsed" desc="aux-def2-SVP">
    auxbasisfile = ofstream("aux-def2-SVP.xml");
    auxbasisfile << "<basis name=\"aux-def2-SVP\">" << std::endl;
    auxbasisfile << "  <!--Basis set created by xtp_basisset from . at Tue Sep  6 17:03:05 2016-->" << std::endl;
    auxbasisfile << "  <element name=\"H\">" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
    auxbasisfile << "      <constant decay=\"1.567529e+01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"1.868860e-02\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"3.606358e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"6.316700e-02\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"1.208002e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"1.204609e-01\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
    auxbasisfile << "      <constant decay=\"4.726794e-01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"5.924850e-02\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
    auxbasisfile << "      <constant decay=\"2.018100e-01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"5.127200e-03\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"P\">" << std::endl;
    auxbasisfile << "      <constant decay=\"2.028137e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"5.358730e-01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"D\">" << std::endl;
    auxbasisfile << "      <constant decay=\"2.216512e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"3.311600e-03\" type=\"D\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "  </element>" << std::endl;
    auxbasisfile << "  <element name=\"C\">" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
    auxbasisfile << "      <constant decay=\"1.861092e+03\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"7.441710e-02\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"6.429940e+02\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"1.653957e-01\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"2.351106e+02\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"5.576484e-01\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"9.070289e+01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"1.310830e+00\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"3.677946e+01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"2.169468e+00\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"1.560463e+01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"1.766885e+00\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"6.890729e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"2.930769e-01\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
    auxbasisfile << "      <constant decay=\"3.147885e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"-1.708702e-01\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
    auxbasisfile << "      <constant decay=\"1.477729e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"1.641553e-01\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
    auxbasisfile << "      <constant decay=\"7.076466e-01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"4.149941e-01\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
    auxbasisfile << "      <constant decay=\"3.430122e-01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"1.624366e-01\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
    auxbasisfile << "      <constant decay=\"1.669453e-01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"2.076750e-02\" type=\"S\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"P\">" << std::endl;
    auxbasisfile << "      <constant decay=\"1.354729e+01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"-2.064770e-02\" type=\"P\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"5.466942e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"-1.152820e-02\" type=\"P\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"P\">" << std::endl;
    auxbasisfile << "      <constant decay=\"2.175172e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"4.559140e-02\" type=\"P\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"P\">" << std::endl;
    auxbasisfile << "      <constant decay=\"8.582194e-01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"2.836000e-03\" type=\"P\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"P\">" << std::endl;
    auxbasisfile << "      <constant decay=\"3.376720e-01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"1.818750e-02\" type=\"P\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"D\">" << std::endl;
    auxbasisfile << "      <constant decay=\"5.928725e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"-2.259480e-02\" type=\"D\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"1.980921e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"-4.768270e-02\" type=\"D\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"D\">" << std::endl;
    auxbasisfile << "      <constant decay=\"8.055417e-01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"-3.653720e-02\" type=\"D\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"D\">" << std::endl;
    auxbasisfile << "      <constant decay=\"3.531244e-01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"-1.454170e-02\" type=\"D\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"F\">" << std::endl;
    auxbasisfile << "      <constant decay=\"1.675563e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"8.879800e-03\" type=\"F\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "      <constant decay=\"5.997536e-01\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"6.990300e-03\" type=\"F\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "    <shell scale=\"1.0\" type=\"G\">" << std::endl;
    auxbasisfile << "      <constant decay=\"1.002460e+00\">" << std::endl;
    auxbasisfile << "        <contractions factor=\"-2.219200e-03\" type=\"G\"/>" << std::endl;
    auxbasisfile << "      </constant>" << std::endl;
    auxbasisfile << "    </shell>" << std::endl;
    auxbasisfile << "  </element>" << std::endl;
    auxbasisfile << "</basis>" << std::endl;
    auxbasisfile.close();
    // </editor-fold>
    
    // ****** Read Input Files, Initiate objects ******
    
//    std::string moleculefilename = "methane.xyz";
//    std::string dftbasisfilename = "3-21G.xml";
//    std::string auxbasisfilename = "3-21G.xml";
    
    std::string moleculefilename = "methane.xyz";
    std::string dftbasisfilename = "def2-SVP.xml";
    std::string auxbasisfilename = "aux-def2-SVP.xml";

    Orbitals orbitals;
    orbitals.LoadFromXYZ(moleculefilename);
    orbitals.setDFTbasis(dftbasisfilename);

    BasisSet basis;
    basis.LoadBasisSet(auxbasisfilename);

    AOBasis aobasis;
    aobasis.AOBasisFill(basis, orbitals.QMAtoms());

    // ****** Load Orbitals ******

    // Instead of loading the orbitals from a .ctp file,
    // we load them from a .xyz file and set the MO coefficients
    // and MO energies by hand.
    orbitals.setBasisSetSize(17);
    orbitals.setNumberOfLevels(4, 13);
    
    Eigen::MatrixXd& MOs = orbitals.MOCoefficients();
    MOs = Eigen::MatrixXd::Zero(17, 17);
    MOs << -0.00761992, -4.69664e-13, 8.35009e-15, -1.15214e-14, -0.0156169, -2.23157e-12, 1.52916e-14, 2.10997e-15, 8.21478e-15, 3.18517e-15, 2.89043e-13, -0.00949189, 1.95787e-12, 1.22168e-14, -2.63092e-15, -0.22227, 1.00844,
            0.233602, -3.18103e-12, 4.05093e-14, -4.70943e-14, 0.1578, 4.75897e-11, -1.87447e-13, -1.02418e-14, 6.44484e-14, -2.6602e-14, 6.5906e-12, -0.281033, -6.67755e-12, 2.70339e-14, -9.78783e-14, -1.94373, -0.36629,
            -1.63678e-13, -0.22745, -0.054851, 0.30351, 3.78688e-11, -0.201627, -0.158318, -0.233561, -0.0509347, -0.650424, 0.452606, -5.88565e-11, 0.453936, -0.165715, -0.619056, 7.0149e-12, 2.395e-14,
            -4.51653e-14, -0.216509, 0.296975, -0.108582, 3.79159e-11, -0.199301, 0.283114, -0.0198557, 0.584622, 0.275311, 0.461431, -5.93732e-11, 0.453057, 0.619523, 0.166374, 7.13235e-12, 2.56811e-14,
            -9.0903e-14, -0.21966, -0.235919, -0.207249, 3.75979e-11, -0.199736, -0.122681, 0.255585, -0.534902, 0.362837, 0.461224, -5.91028e-11, 0.453245, -0.453298, 0.453695, 7.01644e-12, 2.60987e-14,
            0.480866, 1.8992e-11, -2.56795e-13, 4.14571e-13, 2.2709, 4.78615e-10, -2.39153e-12, -2.53852e-13, -2.15605e-13, -2.80359e-13, 7.00137e-12, 0.145171, -1.96136e-11, -2.24876e-13, -2.57294e-14, 4.04176, 0.193617,
            -1.64421e-12, -0.182159, -0.0439288, 0.243073, 1.80753e-10, -0.764779, -0.600505, -0.885907, 0.0862014, 1.10077, -0.765985, 6.65828e-11, -0.579266, 0.211468, 0.789976, -1.41532e-11, -1.29659e-13,
            -1.64105e-12, -0.173397, 0.23784, -0.0869607, 1.80537e-10, -0.755957, 1.07386, -0.0753135, -0.989408, -0.465933, -0.78092, 6.72256e-11, -0.578145, -0.790571, -0.212309, -1.42443e-11, -1.31306e-13,
            -1.63849e-12, -0.17592, -0.188941, -0.165981, 1.79403e-10, -0.757606, -0.465334, 0.969444, 0.905262, -0.61406, -0.78057, 6.69453e-11, -0.578385, 0.578453, -0.578959, -1.40917e-11, -1.31002e-13,
            0.129798, -0.274485, 0.00256652, -0.00509635, -0.0118465, 0.141392, -0.000497905, -0.000510338, -0.000526798, -0.00532572, 0.596595, 0.65313, -0.964582, -0.000361559, -0.000717866, -0.195084, 0.0246232,
            0.0541331, -0.255228, 0.00238646, -0.0047388, -0.88576, 1.68364, -0.00592888, -0.00607692, -9.5047e-05, -0.000960887, 0.10764, -0.362701, 1.53456, 0.000575205, 0.00114206, -0.793844, -0.035336,
            0.129798, 0.0863299, -0.0479412, 0.25617, -0.0118465, -0.0464689, 0.0750316, 0.110468, -0.0436647, -0.558989, -0.203909, 0.65313, 0.320785, 0.235387, 0.878697, -0.195084, 0.0246232,
            0.0541331, 0.0802732, -0.0445777, 0.238198, -0.88576, -0.553335, 0.893449, 1.31541, -0.00787816, -0.100855, -0.0367902, -0.362701, -0.510338, -0.374479, -1.39792, -0.793844, -0.035336,
            0.129798, 0.0927742, -0.197727, -0.166347, -0.0118465, -0.0473592, 0.0582544, -0.119815, -0.463559, 0.320126, -0.196433, 0.65313, 0.321765, 0.643254, -0.642737, -0.195084, 0.0246232,
            0.0541331, 0.0862654, -0.183855, -0.154677, -0.88576, -0.563936, 0.693672, -1.42672, -0.0836372, 0.0577585, -0.0354411, -0.362701, -0.511897, -1.02335, 1.02253, -0.793844, -0.035336,
            0.129798, 0.0953806, 0.243102, -0.0847266, -0.0118465, -0.0475639, -0.132788, 0.00985812, 0.507751, 0.244188, -0.196253, 0.65313, 0.322032, -0.87828, -0.235242, -0.195084, 0.0246232,
            0.0541331, 0.088689, 0.226046, -0.0787824, -0.88576, -0.566373, -1.58119, 0.117387, 0.0916104, 0.0440574, -0.0354087, -0.362701, -0.512321, 1.39726, 0.374248, -0.793844, -0.035336;

    Eigen::VectorXd& mo_energies = orbitals.MOEnergies();
    mo_energies = Eigen::VectorXd::Zero(17);
    mo_energies << -0.612601,
            -0.341755,
            -0.341755,
            -0.341755,
            0.137304,
            0.16678,
            0.16678,
            0.16678,
            0.671592,
            0.671592,
            0.671592,
            0.974255,
            1.01205,
            1.01205,
            1.01205,
            1.64823,
            19.4429;
    
    // ****** Set-up Exchange-Correlation Potential Matrix ******
    
    Eigen::MatrixXd vxc = Eigen::MatrixXd::Zero(17, 17);
    vxc << -0.431767, -0.131967, -1.18442e-13, -1.26466e-13, -1.02288e-13, -0.10626, -3.92543e-13, -3.95555e-13, -3.91314e-13, -0.0116413, -0.0478527, -0.0116413, -0.0478527, -0.0116413, -0.0478527, -0.0116413, -0.0478527,
            -0.131967, -0.647421, 2.51812e-13, 1.39542e-13, 1.8995e-13, -0.465937, 1.53843e-14, -9.48305e-15, -5.94885e-15, -0.119833, -0.241381, -0.119833, -0.241381, -0.119833, -0.241381, -0.119833, -0.241381,
            -1.18442e-13, 2.51812e-13, -0.637843, 1.33983e-13, 9.6584e-14, 5.89028e-14, -0.296161, -4.97511e-13, -5.21849e-13, -0.103175, -0.0760583, -0.103175, -0.0760583, 0.103175, 0.0760583, 0.103175, 0.0760583,
            -1.26466e-13, 1.39542e-13, 1.33983e-13, -0.637843, 2.54059e-13, 4.95922e-15, -4.97536e-13, -0.296161, -4.56739e-13, -0.103175, -0.0760583, 0.103175, 0.0760583, 0.103175, 0.0760583, -0.103175, -0.0760583,
            -1.02288e-13, 1.8995e-13, 9.6584e-14, 2.54059e-13, -0.637843, 2.5538e-14, -5.21859e-13, -4.56639e-13, -0.296161, -0.103175, -0.0760583, 0.103175, 0.0760583, -0.103175, -0.0760583, 0.103175, 0.0760583,
            -0.10626, -0.465937, 5.89028e-14, 4.95922e-15, 2.5538e-14, -0.492236, -6.90263e-14, -8.71169e-14, -9.02027e-14, -0.180782, -0.300264, -0.180782, -0.300264, -0.180782, -0.300264, -0.180782, -0.300264,
            -3.92543e-13, 1.53843e-14, -0.296161, -4.97536e-13, -5.21859e-13, -6.90263e-14, -0.375768, -4.87264e-14, -6.59106e-14, -0.147757, -0.122087, -0.147757, -0.122087, 0.147757, 0.122087, 0.147757, 0.122087,
            -3.95555e-13, -9.48305e-15, -4.97511e-13, -0.296161, -4.56639e-13, -8.71169e-14, -4.87264e-14, -0.375768, -2.38269e-14, -0.147757, -0.122087, 0.147757, 0.122087, 0.147757, 0.122087, -0.147757, -0.122087,
            -3.91314e-13, -5.94885e-15, -5.21849e-13, -4.56739e-13, -0.296161, -9.02027e-14, -6.59106e-14, -2.38269e-14, -0.375768, -0.147757, -0.122087, 0.147757, 0.122087, -0.147757, -0.122087, 0.147757, 0.122087,
            -0.0116413, -0.119833, -0.103175, -0.103175, -0.103175, -0.180782, -0.147757, -0.147757, -0.147757, -0.571548, -0.31776, -0.00435077, -0.061678, -0.00435077, -0.061678, -0.00435077, -0.061678,
            -0.0478527, -0.241381, -0.0760583, -0.0760583, -0.0760583, -0.300264, -0.122087, -0.122087, -0.122087, -0.31776, -0.353709, -0.061678, -0.149893, -0.061678, -0.149893, -0.061678, -0.149893,
            -0.0116413, -0.119833, -0.103175, 0.103175, 0.103175, -0.180782, -0.147757, 0.147757, 0.147757, -0.00435077, -0.061678, -0.571548, -0.31776, -0.00435077, -0.061678, -0.00435077, -0.061678,
            -0.0478527, -0.241381, -0.0760583, 0.0760583, 0.0760583, -0.300264, -0.122087, 0.122087, 0.122087, -0.061678, -0.149893, -0.31776, -0.353709, -0.061678, -0.149893, -0.061678, -0.149893,
            -0.0116413, -0.119833, 0.103175, 0.103175, -0.103175, -0.180782, 0.147757, 0.147757, -0.147757, -0.00435077, -0.061678, -0.00435077, -0.061678, -0.571548, -0.31776, -0.00435077, -0.061678,
            -0.0478527, -0.241381, 0.0760583, 0.0760583, -0.0760583, -0.300264, 0.122087, 0.122087, -0.122087, -0.061678, -0.149893, -0.061678, -0.149893, -0.31776, -0.353709, -0.061678, -0.149893,
            -0.0116413, -0.119833, 0.103175, -0.103175, 0.103175, -0.180782, 0.147757, -0.147757, 0.147757, -0.00435077, -0.061678, -0.00435077, -0.061678, -0.00435077, -0.061678, -0.571548, -0.31776,
            -0.0478527, -0.241381, 0.0760583, -0.0760583, 0.0760583, -0.300264, 0.122087, -0.122087, 0.122087, -0.061678, -0.149893, -0.061678, -0.149893, -0.061678, -0.149893, -0.31776, -0.353709;
    vxc = MOs.transpose() * vxc * MOs;
    
    // We'll need Vxc later for eq. 6
    
    // ****** Set-up Overlap Matrix ******

    AOOverlap ov;
    ov.Fill(aobasis);
    
    // ****** Set-up Coulomb Matrix ******

    AOCoulomb cou;
    cou.Fill(aobasis);
    
    // We'll need the overlap and coulomb matrices for the three-center stuff
    
    // ****** Set-up Three-center Matrix ******
    
    TCMatrix_gwbse Mmn;
    Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
    Mmn.Fill(aobasis, aobasis, MOs);
    Mmn.MultiplyRightWithAuxMatrix(cou.Pseudo_InvSqrt_GWBSE(ov, 1e-7));

    // ***** BSE *****

    RPA rpa;
    rpa.configure(4, 0, 16);
    PPM ppm;
    Eigen::VectorXd screen_r = Eigen::VectorXd::Zero(1);
    screen_r(0) = ppm.getScreening_r();
    Eigen::VectorXd screen_i = Eigen::VectorXd::Zero(1);
    screen_i(0) = ppm.getScreening_i();
    rpa.setScreening(screen_r, screen_i);
    rpa.calculate_epsilon(mo_energies, Mmn);
    ppm.PPM_construct_parameters(rpa);
    Mmn.MultiplyRightWithAuxMatrix(ppm.getPpm_phi());
    votca::ctp::Logger log;
    Sigma sigma = Sigma(&log);
    sigma.configure(4, 0, 16, 20, 1e-5);
    sigma.setDFTdata(0.0, &vxc, &mo_energies);
    sigma.setGWAEnergies(mo_energies);
    sigma.CalcdiagElements(Mmn, ppm);
    sigma.CalcOffDiagElements(Mmn, ppm);
    
    // ****** Using new objects ******
    
    // Prepare spectral decomposition
    RPA_Spectral rpa_spectral = RPA_Spectral();
    rpa_spectral.configure_bse(4, 0, 16, 1);
    rpa_spectral.configure_qp(4, 0, 16);
    rpa_spectral.setGWAEnergies(mo_energies);
    rpa_spectral.prepare_decomp(Mmn);

    // Refine GWA energies
    Sigma_Spectral sigma_spectral = Sigma_Spectral();
    sigma_spectral.configure_bse(4, 0, 16, 1);
    sigma_spectral.configure_qp(4, 0, 16);
    sigma_spectral.configure_g_iter(40, 1e-5);
    sigma_spectral.setGWAEnergies(mo_energies);
    //sigma_spectral.refine_energies(Mmn, rpa_spectral, orbitals.getScaHFX(), vxc, orbitals.MOEnergies());

    sigma_spectral.setHedin(false);
    sigma_spectral.compute_sigma(Mmn, rpa_spectral, orbitals.getScaHFX());
    
    std::cout
            << "Correlation:"    << std::endl
            << "sigma:"          << std::endl << sigma.get_sigma_c().diagonal() << std::endl
            << "sigma_spectral:" << std::endl << sigma_spectral.get_sigma_c().diagonal() << std::endl
            << std::endl;
    
    sigma_spectral.setHedin(true);
    sigma_spectral.compute_sigma(Mmn, rpa_spectral, orbitals.getScaHFX());
    
    std::cout
            << "Correlation:"    << std::endl
            << "sigma:"          << std::endl << sigma.get_sigma_c().diagonal() << std::endl
            << "sigma_spectral:" << std::endl << sigma_spectral.get_sigma_c().diagonal() << std::endl
            << std::endl;
    
    //std::cout << "d_sigma_c:" << std::endl << sigma_spectral.get_sigma_c() - sigma.get_sigma_c() << std::endl;

}

BOOST_AUTO_TEST_SUITE_END()