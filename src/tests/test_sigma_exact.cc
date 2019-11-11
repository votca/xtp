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

#define BOOST_TEST_MODULE sigma_test
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/sigma_exact.h>
#include <votca/xtp/threecenter.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(sigma_test)

BOOST_AUTO_TEST_CASE(sigma_full) {

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

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile("molecule.xyz");
  BasisSet basis;
  basis.Load("3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::VectorXd mo_energy = Eigen::VectorXd::Zero(17);
  mo_energy << 0.0468207, 0.0907801, 0.0907801, 0.104563, 0.592491, 0.663355,
      0.663355, 0.768373, 1.69292, 1.97724, 1.97724, 2.50877, 2.98732, 3.4418,
      3.4418, 4.81084, 17.1838;

  Eigen::MatrixXd MOs = Eigen::MatrixXd::Zero(17, 17);
  MOs << 0.0185815, 2.9133e-17, 8.49354e-17, -0.00312916, 0.0420075,
      1.11356e-16, 1.85886e-17, -0.0334732, 0.0485113, -8.71556e-18,
      -3.79994e-17, -0.0346485, -0.0248392, -3.32286e-22, -4.62643e-17,
      0.0144472, 0.996183, 0.166534, 2.80578e-16, 7.85515e-16, -0.0299557,
      0.409156, 1.10045e-15, 1.67666e-16, -0.336895, 0.608646, -1.00343e-16,
      -5.04622e-16, -0.437568, -0.299751, -2.77682e-17, -6.06812e-16, 0.176594,
      -0.0866677, 0.0010572, -0.0210402, -0.0345975, 0.035778, 0.0611836,
      0.0374747, -0.154443, 0.0892921, 0.0842611, -0.35309, -0.0759572,
      0.278374, -0.409082, -0.64367, 0.308248, -0.261525, -0.000315534,
      -0.0010572, -0.0404824, 0.000922593, -0.035778, -0.0611836, -0.115015,
      -0.109676, -0.0892921, -0.0842611, -0.242326, 0.267806, -0.278374,
      0.409082, -0.0548848, 0.711558, 0.261525, 0.000315534, 0.0010572,
      -0.0194422, 0.0355201, 0.035778, 0.0611836, -0.152489, 0.0447677,
      0.0892921, 0.0842611, 0.110764, 0.343764, 0.278374, -0.409082, 0.588785,
      0.403311, -0.261525, -0.000315534, -0.823783, -9.8891e-16, -3.34692e-15,
      0.103497, -0.277613, -5.51463e-16, 4.95594e-17, 0.163544, 0.121215,
      -7.23985e-17, -2.63149e-16, -0.259891, -0.284396, -1.74149e-16,
      -6.65818e-16, 0.208987, 0.00782842, -0.0333718, 0.22696, 0.373203,
      -0.337332, -0.251625, -0.144131, 0.594004, -0.329076, -0.0456626, 0.18588,
      0.0399869, -0.0631275, -0.0704844, -0.231899, 0.111054, -0.189161,
      0.000129868, 0.0333718, 0.436683, -0.009952, 0.337332, 0.251625, 0.442357,
      0.421824, 0.329076, 0.0456626, 0.12757, -0.140984, 0.0631275, 0.0704844,
      -0.0197737, 0.256357, 0.189161, -0.000129868, -0.0333718, 0.209723,
      -0.383155, -0.337332, -0.251625, 0.586489, -0.172181, -0.329076,
      -0.0456626, -0.0583106, -0.180971, -0.0631275, -0.0704844, 0.212125,
      0.145303, -0.189161, 0.000129868, -0.00177478, 0.0553645, -0.00126176,
      -0.0164247, 0.23154, -0.262519, -0.250334, -0.0135392, -0.429472, 0.45567,
      -0.503583, -0.223493, -0.211802, -0.020461, 0.265268, 0.0023362,
      -0.00241145, 0.294363, -0.686239, 0.0156394, 0.204055, -0.360136,
      0.267096, 0.254698, 0.074687, -0.0228668, 0.132236, -0.14614, -0.174986,
      -0.185046, -0.0109958, 0.142556, 0.0661743, 0.0022999, -0.00177478,
      -0.0265895, 0.0485779, -0.0164247, 0.23154, 0.348055, -0.102182,
      -0.0135392, -0.429472, 0.208281, 0.646413, -0.223493, -0.211802,
      -0.219498, -0.150354, 0.0023362, -0.00241145, 0.294363, 0.329576,
      -0.60212, 0.204055, -0.360136, -0.354123, 0.103963, 0.074687, -0.0228668,
      0.0604434, 0.18759, -0.174986, -0.185046, -0.117959, -0.0808008,
      0.0661743, 0.0022999, -0.00177478, -0.028775, -0.0473162, -0.0164247,
      0.23154, -0.0855356, 0.352515, -0.0135392, -0.429472, -0.663951, -0.14283,
      -0.223493, -0.211802, 0.239959, -0.114914, 0.0023362, -0.00241145,
      0.294363, 0.356664, 0.586481, 0.204055, -0.360136, 0.0870267, -0.358661,
      0.074687, -0.0228668, -0.192679, -0.0414494, -0.174986, -0.185046,
      0.128955, -0.0617554, 0.0661743, 0.0022999, 0.00741062, -3.87173e-16,
      -4.31863e-16, 0.0468488, -0.0476991, 7.27357e-16, 1.23654e-15, -0.43422,
      -0.159247, -4.34945e-17, 1.2743e-16, 0.503528, -0.228856, -7.97629e-17,
      -2.53026e-16, 0.689669, -0.00301027, 0.173046, 7.91486e-15, 8.39419e-15,
      -0.717804, -0.0195249, -1.10754e-15, -1.66789e-15, 0.551371, 0.0684292,
      4.15572e-17, -1.84233e-16, 0.0105378, -0.148396, -1.63792e-16,
      -4.6499e-16, 0.351571, 0.00210309;

  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);

  Logger log;
  RPA rpa(log, Mmn);
  rpa.setRPAInputEnergies(mo_energy);
  rpa.configure(4, 0, 16);

  Sigma_Exact sigma = Sigma_Exact(Mmn, rpa);
  Sigma_Exact::options opt;
  opt.homo = 4;
  opt.qpmin = 0;
  opt.qpmax = 16;
  opt.rpamin = 0;
  opt.rpamax = 16;
  sigma.configure(opt);

  Eigen::MatrixXd x = sigma.CalcExchange();

  Eigen::MatrixXd x_ref = Eigen::MatrixXd::Zero(17, 17);
  x_ref << -0.00370412, -1.59566e-08, -6.75073e-09, -0.000684112, 0.0188911,
      -4.4624e-09, 1.22949e-08, -0.000795391, 0.00504307, 1.72119e-09,
      -9.19106e-09, 0.0140303, 0.01811, -2.89711e-09, 3.40918e-09, -0.00266758,
      2.00221e-05, -1.59566e-08, -0.00682224, -1.70145e-09, 2.72016e-09,
      1.6257e-07, 0.00561074, 0.00559988, 3.98481e-09, 2.85634e-08, 0.00746612,
      -0.00788248, 1.34435e-07, 1.84461e-07, -0.00064411, 0.0064346,
      -2.38641e-08, -1.82782e-09, -6.75073e-09, -1.70145e-09, -0.00682222,
      7.21988e-10, 6.40818e-08, -0.00559984, 0.00561072, 1.36378e-09,
      1.16052e-08, 0.00788246, 0.00746608, 5.28253e-08, 7.22461e-08,
      -0.00643455, -0.000644091, -9.27453e-09, -7.38479e-10, -0.000684112,
      2.72016e-09, 7.21988e-10, -0.00595334, -0.00819783, 6.55678e-09,
      -4.38787e-09, -0.0026928, 0.000460416, -8.03933e-10, 4.27309e-09,
      -0.00228116, -0.0103568, -2.40859e-09, 1.68126e-09, -0.00102159,
      -0.000111773, 0.0188911, 1.6257e-07, 6.40818e-08, -0.00819783, -0.209862,
      8.00987e-08, -1.22485e-07, -0.0159641, -0.0364151, -1.40586e-08,
      1.05127e-07, -0.154293, -0.21813, 5.43914e-09, -5.48907e-09, 0.00986884,
      5.90199e-05, -4.4624e-09, 0.00561074, -0.00559984, 6.55678e-09,
      8.00987e-08, -0.0291401, -5.4019e-09, 8.80569e-09, 1.14844e-08,
      0.000642242, 0.024554, 5.87742e-08, 8.70793e-08, -0.0162291, -0.0198789,
      -2.58402e-09, 4.9345e-10, 1.22949e-08, 0.00559988, 0.00561072,
      -4.38787e-09, -1.22485e-07, -5.4019e-09, -0.0291401, -1.60476e-08,
      -2.48907e-08, -0.024554, 0.000642261, -9.04449e-08, -1.27343e-07,
      0.0198789, -0.0162291, -2.8592e-09, 4.71519e-10, -0.000795391,
      3.98481e-09, 1.36378e-09, -0.0026928, -0.0159641, 8.80569e-09,
      -1.60476e-08, -0.0115818, -0.00171152, -2.83121e-09, 7.39872e-09,
      -0.00921298, -0.0142798, -2.05408e-09, -4.97901e-09, -0.0100267,
      -9.13292e-05, 0.00504307, 2.85634e-08, 1.16052e-08, 0.000460416,
      -0.0364151, 1.14844e-08, -2.48907e-08, -0.00171152, -0.0090557,
      -4.06853e-09, 1.74842e-08, -0.0276487, -0.0363408, 4.8481e-09,
      -6.02765e-09, 0.0021409, -2.62996e-05, 1.72119e-09, 0.00746612,
      0.00788246, -8.03933e-10, -1.40586e-08, 0.000642242, -0.024554,
      -2.83121e-09, -4.06853e-09, -0.0261176, 9.98743e-09, -2.24351e-08,
      -2.98729e-08, 0.0176435, -0.0136511, 5.75818e-09, 1.88005e-09,
      -9.19106e-09, -0.00788248, 0.00746608, 4.27309e-09, 1.05127e-07, 0.024554,
      0.000642261, 7.39872e-09, 1.74842e-08, 9.98743e-09, -0.0261176,
      8.24848e-08, 1.17485e-07, 0.0136511, 0.0176435, -8.51658e-09,
      -8.53563e-10, 0.0140303, 1.34435e-07, 5.28253e-08, -0.00228116, -0.154293,
      5.87742e-08, -9.04449e-08, -0.00921298, -0.0276487, -2.24351e-08,
      8.24848e-08, -0.117991, -0.162036, 6.73442e-09, -6.66614e-10, 0.0100442,
      0.000214376, 0.01811, 1.84461e-07, 7.22461e-08, -0.0103568, -0.21813,
      8.70793e-08, -1.27343e-07, -0.0142798, -0.0363408, -2.98729e-08,
      1.17485e-07, -0.162036, -0.234311, 6.81854e-09, -4.59511e-11, 0.0131863,
      -9.13286e-05, -2.89711e-09, -0.00064411, -0.00643455, -2.40859e-09,
      5.43914e-09, -0.0162291, 0.0198789, -2.05408e-09, 4.8481e-09, 0.0176435,
      0.0136511, 6.73442e-09, 6.81854e-09, -0.0242219, 6.69066e-09, -1.5142e-09,
      -5.43482e-10, 3.40918e-09, 0.0064346, -0.000644091, 1.68126e-09,
      -5.48907e-09, -0.0198789, -0.0162291, -4.97901e-09, -6.02765e-09,
      -0.0136511, 0.0176435, -6.66614e-10, -4.59511e-11, 6.69066e-09,
      -0.0242219, -9.38585e-09, 1.1836e-09, -0.00266758, -2.38641e-08,
      -9.27453e-09, -0.00102159, 0.00986884, -2.58402e-09, -2.8592e-09,
      -0.0100267, 0.0021409, 5.75818e-09, -8.51658e-09, 0.0100442, 0.0131863,
      -1.5142e-09, -9.38585e-09, -0.0129496, -0.000241474, 2.00221e-05,
      -1.82782e-09, -7.38479e-10, -0.000111773, 5.90199e-05, 4.9345e-10,
      4.71519e-10, -9.13292e-05, -2.62996e-05, 1.88005e-09, -8.53563e-10,
      0.000214376, -9.13286e-05, -5.43482e-10, 1.1836e-09, -0.000241474,
      -0.000593496;

  bool check_x = x_ref.isApprox(x, 1e-5);
  if (!check_x) {
    cout << "Sigma X" << endl;
    cout << x << endl;
    cout << "Sigma X ref" << endl;
    cout << x_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_x, true);

  sigma.PrepareScreening();
  sigma.set_COHSEX(false);

  Eigen::VectorXd c_diag = sigma.CalcCorrelationDiag(mo_energy);
  Eigen::MatrixXd c_off = sigma.CalcCorrelationOffDiag(mo_energy);

  c_off.diagonal() = c_diag;
  Eigen::MatrixXd c_ref = Eigen::MatrixXd::Zero(17, 17);
  c_ref << 0.000495273, 5.25897e-10, 1.14634e-09, -0.00023772, -0.00170784, 2.17043e-09,
      2.21014e-09, 0.0732799, 0.00236963, 7.93752e-10, -3.48238e-09, 0.00389167,
      0.00385436, -1.30855e-08, 7.63152e-09, -0.00925048, 0.000467863, 5.25897e-10,
      -0.000254943, -5.93472e-10, -1.8009e-09, 2.50202e-08, 0.00173095, 0.0017276,
      -1.16534e-06, 7.89587e-09, 0.00836747, -0.00883409, 1.54847e-07, 2.25869e-07,
      -0.00225727, 0.0225496, -9.56068e-08, -4.6909e-09, 1.14634e-09, -5.93472e-10,
      -0.000254935, -9.16956e-11, 6.17787e-09, -0.00172757, 0.00173093, -5.99279e-07,
      4.59628e-09, 0.00883405, 0.0083674, 5.84156e-08, 8.17389e-08, -0.0225493,
      -0.00225714, -9.86177e-09, -1.80381e-09, -0.00023772, -1.8009e-09, -9.16956e-11,
      -1.94885e-05, 0.00250376, 3.84451e-09, 2.74976e-09, 0.0774164, 0.00162667,
      -1.60759e-09, -1.79567e-09, 0.00369254, -0.000977076, 1.07776e-08, 2.14288e-08,
      -0.00635686, -3.46275e-05, -0.00170784, 2.50202e-08, 6.17787e-09, 0.00250376,
      -0.00773778, 1.36324e-08, -3.98219e-08, -0.450133, -0.00688453, -5.0166e-09,
      5.81583e-08, -0.0987464, -0.154698, 2.07321e-07, -1.48808e-08, 0.116362,
      0.00188979, 2.17043e-09, 0.00173095, -0.00172757, 3.84451e-09, 1.36324e-08,
      -0.0325109, -9.03448e-09, 7.78354e-08, 2.07287e-08, 0.0014161, 0.0541403,
      7.7947e-08, 2.0469e-07, -0.176466, -0.216152, -6.16114e-08, -1.735e-10,
      2.21014e-09, 0.0017276, 0.00173093, 2.74976e-09, -3.98219e-08, -9.03448e-09,
      -0.0325109, -7.63099e-07, -5.74164e-08, -0.0541403, 0.00141615, -2.29803e-07,
      -2.7685e-07, 0.216152, -0.176466, -7.62921e-08, 1.01264e-09, 0.0732799,
      -1.16534e-06, -5.99279e-07, 0.0774164, -0.450133, 7.78354e-08, -7.63099e-07,
      0.577854, -0.0434323, 9.27481e-07, -4.34245e-07, -0.43091, -0.509184,
      9.55096e-07, -4.23331e-07, 0.348581, -0.0106311, 0.00236963, 7.89587e-09,
      4.59628e-09, 0.00162667, -0.00688453, 2.07287e-08, -5.74164e-08, -0.0434323,
      -0.0509681, -4.97831e-08, -6.40814e-09, -0.0242454, -0.020095, 2.30051e-07,
      -6.03896e-08, 0.056842, -0.00858586, 7.93752e-10, 0.00836747, 0.00883405,
      -1.60759e-09, -5.0166e-09, 0.0014161, -0.0541403, 9.27481e-07, -4.97831e-08,
      -0.161677, -5.67937e-09, -6.79682e-08, -4.25427e-08, 0.21982, -0.170079,
      -1.04038e-06, -2.00888e-09, -3.48238e-09, -0.00883409, 0.0083674, -1.79567e-09,
      5.81583e-08, 0.0541403, 0.00141615, -4.34245e-07, -6.40814e-09, -5.67937e-09,
      -0.161676, 1.31124e-07, 2.36986e-07, 0.170079, 0.21982, 1.53735e-07,
      -6.59809e-09, 0.00389167, 1.54847e-07, 5.84156e-08, 0.00369254, -0.0987464,
      7.7947e-08, -2.29803e-07, -0.43091, -0.0242454, -6.79682e-08, 1.31124e-07,
      -0.342342, -0.47093, 4.83466e-07, -1.98528e-07, 0.333323, 0.016858,
      0.00385436, 2.25869e-07, 8.17389e-08, -0.000977076, -0.154698, 2.0469e-07,
      -2.7685e-07, -0.509184, -0.020095, -4.25427e-08, 2.36986e-07, -0.47093,
      -0.750009, 6.0018e-07, 1.39042e-07, 0.386863, 0.0253594, -1.30855e-08,
      -0.00225727, -0.0225493, 1.07776e-08, 2.07321e-07, -0.176466, 0.216152,
      9.55096e-07, 2.30051e-07, 0.21982, 0.170079, 4.83466e-07, 6.0018e-07,
      -0.653859, 1.58317e-07, 4.20836e-07, 7.24923e-09, 7.63152e-09, 0.0225496,
      -0.00225714, 2.14288e-08, -1.48808e-08, -0.216152, -0.176466, -4.23331e-07,
      -6.03896e-08, -0.170079, 0.21982, -1.98528e-07, 1.39042e-07, 1.58317e-07,
      -0.653859, -7.96593e-08, -7.12763e-11, -0.00925048, -9.56068e-08, -9.86177e-09,
      -0.00635686, 0.116362, -6.16114e-08, -7.62921e-08, 0.348581, 0.056842,
      -1.04038e-06, 1.53735e-07, 0.333323, 0.386863, 4.20836e-07, -7.96593e-08,
      -0.389375, -0.0138152, 0.000467863, -4.6909e-09, -1.80381e-09, -3.46275e-05,
      0.00188979, -1.735e-10, 1.01264e-09, -0.0106311, -0.00858586, -2.00888e-09,
      -6.59809e-09, 0.016858, 0.0253594, 7.24923e-09, -7.12763e-11, -0.0138152,
      -0.0743819;

  bool check_c_diag = c_diag.isApprox(c_ref.diagonal(), 1e-5);
  if (!check_c_diag) {
    cout << "Sigma C" << endl;
    cout << c_diag << endl;
    cout << "Sigma C ref" << endl;
    cout << c_ref.diagonal() << endl;
  }
  BOOST_CHECK_EQUAL(check_c_diag, true);

  bool check_c = c_ref.isApprox(c_off, 1e-5);
  if (!check_c) {
    cout << "Sigma C" << endl;
    cout << c_off << endl;
    cout << "Sigma C ref" << endl;
    cout << c_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_c, true);
}

BOOST_AUTO_TEST_SUITE_END()
