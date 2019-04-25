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

#define BOOST_TEST_MODULE gaussian_quadrature_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/ctp/logger.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/gaussian_quadrature.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/threecenter.h>

using namespace votca::xtp;
using namespace std;
BOOST_AUTO_TEST_SUITE(gaussian_quadrature_test)

BOOST_AUTO_TEST_CASE(gaussian_quadrature_full) {

  std::ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  std::ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << endl;
  basisfile << "  <element name=\"H\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"8.245470e-01\">" << endl;
  basisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.831920e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();

  BasisSet basis;
  basis.LoadBasisSet("3-21G.xml");
  AOBasis aobasis;
  Orbitals orbitals;
  orbitals.LoadFromXYZ("molecule.xyz");

  aobasis.AOBasisFill(basis, orbitals.QMAtoms());

  orbitals.setBasisSetSize(17);
  orbitals.setNumberOfOccupiedLevels(4);

  AOKinetic kinetic;
  kinetic.Fill(aobasis);

  AOESP esp;
  esp.Fillnucpotential(aobasis, orbitals.QMAtoms());

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(kinetic.Matrix() +
                                                    esp.Matrix());

  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, es.eigenvectors());
  votca::ctp::Logger log;
  RPA rpa(log,Mmn);
  rpa.setRPAInputEnergies(es.eigenvalues());
  rpa.configure(4, 0, 17 - 1);

  GaussianQuadrature gq = GaussianQuadrature(es.eigenvalues(), Mmn);
  GaussianQuadrature::options opt;

  opt.homo = 4;
  opt.qptotal = 17;
  opt.qpmin = 0;
  opt.order = 12;
  opt.rpamin = 0;

  gq.configure(opt);
  // gq.configure(17,0,4);
  Eigen::MatrixXd result = gq.SigmaGQ(es.eigenvalues(), rpa);
  Eigen::VectorXd result2 = gq.SigmaGQDiag(es.eigenvalues(), rpa);
  Eigen::MatrixXd exactresult = gq.ExactSigmaGQ(es.eigenvalues(), rpa);
 
 
  bool check_c_diag = result2.isApprox(exactresult.diagonal(), 1e-5);

  if (!check_c_diag) {
    std::cout << "Sigma C" << std::endl;
    std::cout << result2 << std::endl;
    std::cout << "Sigma C ref" << std::endl;
    std::cout << exactresult.diagonal() << std::endl;
  }
  BOOST_CHECK_EQUAL(check_c_diag, true);
  bool check_c = exactresult.isApprox(result, 1e-5);
  if (!check_c) {
    std::cout << "Sigma C" << std::endl;
    std::cout << result << std::endl;
    std::cout << "Sigma C ref" << std::endl;
    std::cout << exactresult << std::endl;
  }
  BOOST_CHECK_EQUAL(check_c, true);
}

BOOST_AUTO_TEST_SUITE_END()
