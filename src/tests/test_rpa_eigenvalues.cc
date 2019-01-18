/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE rpa_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/rpa.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/threecenter.h>

using namespace std;
using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(rpa_test)

BOOST_AUTO_TEST_CASE(rpa_eigenvalues) {

  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << endl;
  basisfile << "  <element name=\"H\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"8.245470e-01\">" << endl;
  basisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.831920e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>" << endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>" << endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();

  Orbitals orbitals;
  orbitals.LoadFromXYZ("molecule.xyz");
  BasisSet basis;
  basis.LoadBasisSet("3-21G.xml");

  AOBasis aobasis;
  aobasis.AOBasisFill(basis, orbitals.QMAtoms());

  Eigen::VectorXd eigenvals = Eigen::VectorXd::Zero(17);
  eigenvals << 0.0468207, 0.0907801, 0.0907801, 0.104563, 0.592491, 0.663355,
          0.663355, 0.768373, 1.69292, 1.97724, 1.97724, 2.50877, 2.98732, 3.4418, 3.4418, 4.81084,
          17.1838;

  Eigen::MatrixXd eigenvectors = Eigen::MatrixXd::Zero(17, 17);
  eigenvectors << 0.0185815, 2.9133e-17, 8.49354e-17, -0.00312916, 0.0420075, 1.11356e-16, 1.85886e-17, -0.0334732, 0.0485113, -8.71556e-18, -3.79994e-17, -0.0346485, -0.0248392, -3.32286e-22, -4.62643e-17, 0.0144472, 0.996183,
          0.166534, 2.80578e-16, 7.85515e-16, -0.0299557, 0.409156, 1.10045e-15, 1.67666e-16, -0.336895, 0.608646, -1.00343e-16, -5.04622e-16, -0.437568, -0.299751, -2.77682e-17, -6.06812e-16, 0.176594, -0.0866677,
          0.0010572, -0.0210402, -0.0345975, 0.035778, 0.0611836, 0.0374747, -0.154443, 0.0892921, 0.0842611, -0.35309, -0.0759572, 0.278374, -0.409082, -0.64367, 0.308248, -0.261525, -0.000315534,
          -0.0010572, -0.0404824, 0.000922593, -0.035778, -0.0611836, -0.115015, -0.109676, -0.0892921, -0.0842611, -0.242326, 0.267806, -0.278374, 0.409082, -0.0548848, 0.711558, 0.261525, 0.000315534,
          0.0010572, -0.0194422, 0.0355201, 0.035778, 0.0611836, -0.152489, 0.0447677, 0.0892921, 0.0842611, 0.110764, 0.343764, 0.278374, -0.409082, 0.588785, 0.403311, -0.261525, -0.000315534,
          -0.823783, -9.8891e-16, -3.34692e-15, 0.103497, -0.277613, -5.51463e-16, 4.95594e-17, 0.163544, 0.121215, -7.23985e-17, -2.63149e-16, -0.259891, -0.284396, -1.74149e-16, -6.65818e-16, 0.208987, 0.00782842,
          -0.0333718, 0.22696, 0.373203, -0.337332, -0.251625, -0.144131, 0.594004, -0.329076, -0.0456626, 0.18588, 0.0399869, -0.0631275, -0.0704844, -0.231899, 0.111054, -0.189161, 0.000129868,
          0.0333718, 0.436683, -0.009952, 0.337332, 0.251625, 0.442357, 0.421824, 0.329076, 0.0456626, 0.12757, -0.140984, 0.0631275, 0.0704844, -0.0197737, 0.256357, 0.189161, -0.000129868,
          -0.0333718, 0.209723, -0.383155, -0.337332, -0.251625, 0.586489, -0.172181, -0.329076, -0.0456626, -0.0583106, -0.180971, -0.0631275, -0.0704844, 0.212125, 0.145303, -0.189161, 0.000129868,
          -0.00177478, 0.0553645, -0.00126176, -0.0164247, 0.23154, -0.262519, -0.250334, -0.0135392, -0.429472, 0.45567, -0.503583, -0.223493, -0.211802, -0.020461, 0.265268, 0.0023362, -0.00241145,
          0.294363, -0.686239, 0.0156394, 0.204055, -0.360136, 0.267096, 0.254698, 0.074687, -0.0228668, 0.132236, -0.14614, -0.174986, -0.185046, -0.0109958, 0.142556, 0.0661743, 0.0022999,
          -0.00177478, -0.0265895, 0.0485779, -0.0164247, 0.23154, 0.348055, -0.102182, -0.0135392, -0.429472, 0.208281, 0.646413, -0.223493, -0.211802, -0.219498, -0.150354, 0.0023362, -0.00241145,
          0.294363, 0.329576, -0.60212, 0.204055, -0.360136, -0.354123, 0.103963, 0.074687, -0.0228668, 0.0604434, 0.18759, -0.174986, -0.185046, -0.117959, -0.0808008, 0.0661743, 0.0022999,
          -0.00177478, -0.028775, -0.0473162, -0.0164247, 0.23154, -0.0855356, 0.352515, -0.0135392, -0.429472, -0.663951, -0.14283, -0.223493, -0.211802, 0.239959, -0.114914, 0.0023362, -0.00241145,
          0.294363, 0.356664, 0.586481, 0.204055, -0.360136, 0.0870267, -0.358661, 0.074687, -0.0228668, -0.192679, -0.0414494, -0.174986, -0.185046, 0.128955, -0.0617554, 0.0661743, 0.0022999,
          0.00741062, -3.87173e-16, -4.31863e-16, 0.0468488, -0.0476991, 7.27357e-16, 1.23654e-15, -0.43422, -0.159247, -4.34945e-17, 1.2743e-16, 0.503528, -0.228856, -7.97629e-17, -2.53026e-16, 0.689669, -0.00301027,
          0.173046, 7.91486e-15, 8.39419e-15, -0.717804, -0.0195249, -1.10754e-15, -1.66789e-15, 0.551371, 0.0684292, 4.15572e-17, -1.84233e-16, 0.0105378, -0.148396, -1.63792e-16, -4.6499e-16, 0.351571, 0.00210309;

  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, eigenvectors);

  RPA rpa(Mmn);
  rpa.setRPAInputEnergies(eigenvals);
  rpa.configure(4, 0, 16);

  Eigen::VectorXd rpa_eigenvalues_ref = Eigen::VectorXd(60);
  rpa_eigenvalues_ref << 6.08285e-09, 9.83127e-09, 0.0257186, 0.0399894, 0.0962028, 0.168569,
          0.557901, 0.557907, 0.568082, 0.569721, 0.569732, 0.572575,
          0.615673, 0.615676, 0.660257, 0.677119, 0.677121, 0.720781,
          1.02886, 1.09032, 1.3111, 1.36841, 1.36882, 1.58789,
          1.60137, 1.60138, 1.6413, 1.71693, 1.87207, 1.87207,
          1.8834, 1.88341, 1.88646, 1.88687, 1.9294, 1.9294,
          2.22202, 2.40138, 2.41599, 2.41599, 2.46099, 2.83006,
          2.83007, 2.87933, 2.89433, 2.89433, 2.93895, 3.33655,
          3.33655, 3.34817, 3.34882, 3.34882, 3.35102, 3.39442,
          3.39442, 4.21336, 4.70138, 4.71894, 4.71894, 4.76331;

  rpa_eigensolution sol = rpa.calculate_eigenvalues();
  Eigen::VectorXd rpa_eigenvalues = sol._Omega;

  bool check = rpa_eigenvalues_ref.isApprox(rpa_eigenvalues, 0.0001);

  if (!check) {
    cout << "rpa_eigenvalues" << endl;
    cout << rpa_eigenvalues << endl;
    cout << "rpa_eigenvalues_ref" << endl;
    cout << rpa_eigenvalues_ref << endl;
  }

  BOOST_CHECK_EQUAL(check, 1);

}

BOOST_AUTO_TEST_SUITE_END()
