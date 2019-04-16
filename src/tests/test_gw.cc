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

#define BOOST_TEST_MODULE gw_test
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <votca/xtp/gw.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(gw_test)

BOOST_AUTO_TEST_CASE(gw_full) {

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
  orbitals.LoadFromXYZ("molecule.xyz");
  BasisSet basis;
  basis.LoadBasisSet("3-21G.xml");
  orbitals.setDFTbasisName("3-21G.xml");
  AOBasis aobasis;
  aobasis.AOBasisFill(basis, orbitals.QMAtoms());
  orbitals.setBasisSetSize(17);
  orbitals.setNumberOfOccupiedLevels(4);
  Eigen::MatrixXd& MOs = orbitals.MOCoefficients();
  MOs = Eigen::MatrixXd::Zero(17, 17);
  MOs << -0.00761992, -4.69664e-13, 8.35009e-15, -1.15214e-14, -0.0156169,
      -2.23157e-12, 1.52916e-14, 2.10997e-15, 8.21478e-15, 3.18517e-15,
      2.89043e-13, -0.00949189, 1.95787e-12, 1.22168e-14, -2.63092e-15,
      -0.22227, 1.00844, 0.233602, -3.18103e-12, 4.05093e-14, -4.70943e-14,
      0.1578, 4.75897e-11, -1.87447e-13, -1.02418e-14, 6.44484e-14, -2.6602e-14,
      6.5906e-12, -0.281033, -6.67755e-12, 2.70339e-14, -9.78783e-14, -1.94373,
      -0.36629, -1.63678e-13, -0.22745, -0.054851, 0.30351, 3.78688e-11,
      -0.201627, -0.158318, -0.233561, -0.0509347, -0.650424, 0.452606,
      -5.88565e-11, 0.453936, -0.165715, -0.619056, 7.0149e-12, 2.395e-14,
      -4.51653e-14, -0.216509, 0.296975, -0.108582, 3.79159e-11, -0.199301,
      0.283114, -0.0198557, 0.584622, 0.275311, 0.461431, -5.93732e-11,
      0.453057, 0.619523, 0.166374, 7.13235e-12, 2.56811e-14, -9.0903e-14,
      -0.21966, -0.235919, -0.207249, 3.75979e-11, -0.199736, -0.122681,
      0.255585, -0.534902, 0.362837, 0.461224, -5.91028e-11, 0.453245,
      -0.453298, 0.453695, 7.01644e-12, 2.60987e-14, 0.480866, 1.8992e-11,
      -2.56795e-13, 4.14571e-13, 2.2709, 4.78615e-10, -2.39153e-12,
      -2.53852e-13, -2.15605e-13, -2.80359e-13, 7.00137e-12, 0.145171,
      -1.96136e-11, -2.24876e-13, -2.57294e-14, 4.04176, 0.193617, -1.64421e-12,
      -0.182159, -0.0439288, 0.243073, 1.80753e-10, -0.764779, -0.600505,
      -0.885907, 0.0862014, 1.10077, -0.765985, 6.65828e-11, -0.579266,
      0.211468, 0.789976, -1.41532e-11, -1.29659e-13, -1.64105e-12, -0.173397,
      0.23784, -0.0869607, 1.80537e-10, -0.755957, 1.07386, -0.0753135,
      -0.989408, -0.465933, -0.78092, 6.72256e-11, -0.578145, -0.790571,
      -0.212309, -1.42443e-11, -1.31306e-13, -1.63849e-12, -0.17592, -0.188941,
      -0.165981, 1.79403e-10, -0.757606, -0.465334, 0.969444, 0.905262,
      -0.61406, -0.78057, 6.69453e-11, -0.578385, 0.578453, -0.578959,
      -1.40917e-11, -1.31002e-13, 0.129798, -0.274485, 0.00256652, -0.00509635,
      -0.0118465, 0.141392, -0.000497905, -0.000510338, -0.000526798,
      -0.00532572, 0.596595, 0.65313, -0.964582, -0.000361559, -0.000717866,
      -0.195084, 0.0246232, 0.0541331, -0.255228, 0.00238646, -0.0047388,
      -0.88576, 1.68364, -0.00592888, -0.00607692, -9.5047e-05, -0.000960887,
      0.10764, -0.362701, 1.53456, 0.000575205, 0.00114206, -0.793844,
      -0.035336, 0.129798, 0.0863299, -0.0479412, 0.25617, -0.0118465,
      -0.0464689, 0.0750316, 0.110468, -0.0436647, -0.558989, -0.203909,
      0.65313, 0.320785, 0.235387, 0.878697, -0.195084, 0.0246232, 0.0541331,
      0.0802732, -0.0445777, 0.238198, -0.88576, -0.553335, 0.893449, 1.31541,
      -0.00787816, -0.100855, -0.0367902, -0.362701, -0.510338, -0.374479,
      -1.39792, -0.793844, -0.035336, 0.129798, 0.0927742, -0.197727, -0.166347,
      -0.0118465, -0.0473592, 0.0582544, -0.119815, -0.463559, 0.320126,
      -0.196433, 0.65313, 0.321765, 0.643254, -0.642737, -0.195084, 0.0246232,
      0.0541331, 0.0862654, -0.183855, -0.154677, -0.88576, -0.563936, 0.693672,
      -1.42672, -0.0836372, 0.0577585, -0.0354411, -0.362701, -0.511897,
      -1.02335, 1.02253, -0.793844, -0.035336, 0.129798, 0.0953806, 0.243102,
      -0.0847266, -0.0118465, -0.0475639, -0.132788, 0.00985812, 0.507751,
      0.244188, -0.196253, 0.65313, 0.322032, -0.87828, -0.235242, -0.195084,
      0.0246232, 0.0541331, 0.088689, 0.226046, -0.0787824, -0.88576, -0.566373,
      -1.58119, 0.117387, 0.0916104, 0.0440574, -0.0354087, -0.362701,
      -0.512321, 1.39726, 0.374248, -0.793844, -0.035336;
  Eigen::MatrixXd vxc = Eigen::MatrixXd::Zero(17, 17);
  vxc << -0.431767, -0.131967, -1.18442e-13, -1.26466e-13, -1.02288e-13,
      -0.10626, -3.92543e-13, -3.95555e-13, -3.91314e-13, -0.0116413,
      -0.0478527, -0.0116413, -0.0478527, -0.0116413, -0.0478527, -0.0116413,
      -0.0478527, -0.131967, -0.647421, 2.51812e-13, 1.39542e-13, 1.8995e-13,
      -0.465937, 1.53843e-14, -9.48305e-15, -5.94885e-15, -0.119833, -0.241381,
      -0.119833, -0.241381, -0.119833, -0.241381, -0.119833, -0.241381,
      -1.18442e-13, 2.51812e-13, -0.637843, 1.33983e-13, 9.6584e-14,
      5.89028e-14, -0.296161, -4.97511e-13, -5.21849e-13, -0.103175, -0.0760583,
      -0.103175, -0.0760583, 0.103175, 0.0760583, 0.103175, 0.0760583,
      -1.26466e-13, 1.39542e-13, 1.33983e-13, -0.637843, 2.54059e-13,
      4.95922e-15, -4.97536e-13, -0.296161, -4.56739e-13, -0.103175, -0.0760583,
      0.103175, 0.0760583, 0.103175, 0.0760583, -0.103175, -0.0760583,
      -1.02288e-13, 1.8995e-13, 9.6584e-14, 2.54059e-13, -0.637843, 2.5538e-14,
      -5.21859e-13, -4.56639e-13, -0.296161, -0.103175, -0.0760583, 0.103175,
      0.0760583, -0.103175, -0.0760583, 0.103175, 0.0760583, -0.10626,
      -0.465937, 5.89028e-14, 4.95922e-15, 2.5538e-14, -0.492236, -6.90263e-14,
      -8.71169e-14, -9.02027e-14, -0.180782, -0.300264, -0.180782, -0.300264,
      -0.180782, -0.300264, -0.180782, -0.300264, -3.92543e-13, 1.53843e-14,
      -0.296161, -4.97536e-13, -5.21859e-13, -6.90263e-14, -0.375768,
      -4.87264e-14, -6.59106e-14, -0.147757, -0.122087, -0.147757, -0.122087,
      0.147757, 0.122087, 0.147757, 0.122087, -3.95555e-13, -9.48305e-15,
      -4.97511e-13, -0.296161, -4.56639e-13, -8.71169e-14, -4.87264e-14,
      -0.375768, -2.38269e-14, -0.147757, -0.122087, 0.147757, 0.122087,
      0.147757, 0.122087, -0.147757, -0.122087, -3.91314e-13, -5.94885e-15,
      -5.21849e-13, -4.56739e-13, -0.296161, -9.02027e-14, -6.59106e-14,
      -2.38269e-14, -0.375768, -0.147757, -0.122087, 0.147757, 0.122087,
      -0.147757, -0.122087, 0.147757, 0.122087, -0.0116413, -0.119833,
      -0.103175, -0.103175, -0.103175, -0.180782, -0.147757, -0.147757,
      -0.147757, -0.571548, -0.31776, -0.00435077, -0.061678, -0.00435077,
      -0.061678, -0.00435077, -0.061678, -0.0478527, -0.241381, -0.0760583,
      -0.0760583, -0.0760583, -0.300264, -0.122087, -0.122087, -0.122087,
      -0.31776, -0.353709, -0.061678, -0.149893, -0.061678, -0.149893,
      -0.061678, -0.149893, -0.0116413, -0.119833, -0.103175, 0.103175,
      0.103175, -0.180782, -0.147757, 0.147757, 0.147757, -0.00435077,
      -0.061678, -0.571548, -0.31776, -0.00435077, -0.061678, -0.00435077,
      -0.061678, -0.0478527, -0.241381, -0.0760583, 0.0760583, 0.0760583,
      -0.300264, -0.122087, 0.122087, 0.122087, -0.061678, -0.149893, -0.31776,
      -0.353709, -0.061678, -0.149893, -0.061678, -0.149893, -0.0116413,
      -0.119833, 0.103175, 0.103175, -0.103175, -0.180782, 0.147757, 0.147757,
      -0.147757, -0.00435077, -0.061678, -0.00435077, -0.061678, -0.571548,
      -0.31776, -0.00435077, -0.061678, -0.0478527, -0.241381, 0.0760583,
      0.0760583, -0.0760583, -0.300264, 0.122087, 0.122087, -0.122087,
      -0.061678, -0.149893, -0.061678, -0.149893, -0.31776, -0.353709,
      -0.061678, -0.149893, -0.0116413, -0.119833, 0.103175, -0.103175,
      0.103175, -0.180782, 0.147757, -0.147757, 0.147757, -0.00435077,
      -0.061678, -0.00435077, -0.061678, -0.00435077, -0.061678, -0.571548,
      -0.31776, -0.0478527, -0.241381, 0.0760583, -0.0760583, 0.0760583,
      -0.300264, 0.122087, -0.122087, 0.122087, -0.061678, -0.149893, -0.061678,
      -0.149893, -0.061678, -0.149893, -0.31776, -0.353709;
  vxc = MOs.transpose() * vxc * MOs;
  Eigen::VectorXd mo_energy = Eigen::VectorXd::Zero(17);
  mo_energy << -0.612601, -0.341755, -0.341755, -0.341755, 0.137304, 0.16678,
      0.16678, 0.16678, 0.671592, 0.671592, 0.671592, 0.974255, 1.01205,
      1.01205, 1.01205, 1.64823, 19.4429;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);
  votca::ctp::Logger log;

  GW::options opt;
  opt.ScaHFX = 0;
  opt.homo = 4;
  opt.qpmax = 16;
  opt.qpmin = 0;
  opt.rpamax = 16;
  opt.rpamin = 0;
  opt.gw_sc_max_iterations = 1;
  GW gw(log, Mmn, vxc, mo_energy);
  gw.configure(opt);
  gw.CalculateGWPerturbation();

  Eigen::MatrixXd diag = gw.getGWAResults();
  Eigen::MatrixXd diag_ref = Eigen::MatrixXd::Zero(17, 5);
  diag_ref << -0.612601, -0.898354, 0.0770146, -0.535716, -0.898224, -0.341755,
      -0.690315, 0.0816318, -0.507503, -0.442936, -0.341755, -0.690316,
      0.0816323, -0.507503, -0.442935, -0.341755, -0.690315, 0.0816323,
      -0.507503, -0.442935, 0.137304, -0.419561, -0.0243086, -0.307146,
      0.000580756, 0.16678, -0.162521, -0.0243167, -0.340866, 0.320808, 0.16678,
      -0.162521, -0.0243172, -0.340865, 0.320807, 0.16678, -0.162521,
      -0.0243177, -0.340866, 0.320807, 0.671592, -0.112664, -0.0642448,
      -0.391959, 0.886642, 0.671592, -0.112664, -0.064245, -0.391959, 0.886642,
      0.671592, -0.112664, -0.0642449, -0.391959, 0.886642, 0.974255, -0.17259,
      -0.0675618, -0.523599, 1.2577, 1.01205, -0.131712, 0.0133325, -0.492266,
      1.38594, 1.01205, -0.131712, 0.0133343, -0.492266, 1.38594, 1.01205,
      -0.131712, 0.0133358, -0.492267, 1.38594, 1.64823, -0.10267, -0.0299524,
      -0.414866, 1.93047, 19.4429, -0.0285864, -0.403971, -0.419726, 19.4301;
  std::cout << "okay2.5" << std::endl;
  bool check_diag = diag_ref.isApprox(diag, 1e-5);
  if (!check_diag) {
    cout << "GW energies" << endl;
    cout << diag << endl;
    cout << "GW energies ref" << endl;
    cout << diag_ref << endl;
  }

  gw.CalculateHQP();
  Eigen::MatrixXd offdiag = gw.getHQP();

  Eigen::MatrixXd offdiag_ref = Eigen::MatrixXd::Zero(17, 17);
  offdiag_ref << -0.898224, 4.23678e-07, -1.29796e-07, 1.34937e-07, 0.0277142,
      1.18961e-06, -6.70857e-07, -1.26851e-07, -1.74807e-07, 2.29427e-08,
      -1.71338e-07, -0.0173535, -7.27197e-07, -1.16343e-06, 5.83176e-08,
      0.00771153, -0.008621, 4.23678e-07, -0.442936, 1.10777e-07, -1.45968e-07,
      -1.30323e-07, 0.0444199, 0.000554857, 0.000428135, 9.17985e-05,
      0.000316496, -0.0111513, 1.08681e-07, 0.0169894, -0.000176982,
      -0.000289165, 6.29926e-07, 2.41805e-09, -1.29796e-07, 1.10777e-07,
      -0.442935, 1.73406e-07, 2.76257e-08, -0.000336507, 0.0406831, -0.0178448,
      0.0110944, 0.00116727, 0.000124475, 3.87122e-08, -0.000153526, -0.0169374,
      0.00133593, -4.78396e-07, 4.99241e-09, 1.34937e-07, -1.45968e-07,
      1.73406e-07, -0.442935, -4.34706e-08, 0.000615049, -0.0178393, -0.0406827,
      0.00117032, -0.0110905, -0.000305187, 1.21357e-07, 0.000302317,
      0.00133415, 0.0169342, 1.25256e-07, 1.45205e-09, 0.0277142, -1.30323e-07,
      2.76257e-08, -4.34706e-08, 0.000580756, -1.757e-07, 1.16925e-07,
      -1.14996e-10, -4.40827e-08, 4.31279e-09, -2.11565e-08, -0.0136937,
      9.50177e-08, 1.13961e-07, -4.064e-09, -0.00208028, -0.0204426,
      1.18961e-06, 0.0444199, -0.000336507, 0.000615049, -1.757e-07, 0.320808,
      1.66058e-07, -1.59393e-08, -2.84745e-05, -0.000182575, 0.0132206,
      2.98342e-07, -0.0331867, 5.83871e-05, 0.000126789, -1.83644e-07,
      3.06443e-08, -6.70857e-07, 0.000554857, 0.0406831, -0.0178393,
      1.16925e-07, 1.66058e-07, 0.320807, 2.50411e-09, -0.0114849, -0.0065493,
      -0.000114984, -1.69854e-07, 9.71481e-05, 0.0313476, 0.0109008,
      1.24985e-07, -1.92004e-08, -1.26851e-07, 0.000428135, -0.0178448,
      -0.0406827, -1.14996e-10, -1.59393e-08, 2.50411e-09, 0.320807, 0.00655019,
      -0.011484, -0.000144427, 1.58299e-07, 0.000100547, -0.0109009, 0.0313489,
      -2.78769e-08, 2.5969e-09, -1.74807e-07, 9.17985e-05, 0.0110944,
      0.00117032, -4.40827e-08, -2.84745e-05, -0.0114849, 0.00655019, 0.886642,
      5.29529e-08, -9.26177e-08, 1.79143e-07, -5.20185e-06, 0.00459359,
      -0.000853556, -3.75826e-07, -9.06493e-09, 2.29427e-08, 0.000316496,
      0.00116727, -0.0110905, 4.31279e-09, -0.000182575, -0.0065493, -0.011484,
      5.29529e-08, 0.886642, -1.34264e-08, -1.19884e-07, -4.54307e-05,
      0.0008528, 0.00459404, -1.58829e-07, -6.78976e-09, -1.71338e-07,
      -0.0111513, 0.000124475, -0.000305187, -2.11565e-08, 0.0132206,
      -0.000114984, -0.000144427, -9.26177e-08, -1.34264e-08, 0.886642,
      1.32967e-08, 0.00467098, 1.35001e-05, 4.37768e-05, -6.10398e-07,
      -1.98658e-08, -0.0173535, 1.08681e-07, 3.87122e-08, 1.21357e-07,
      -0.0136937, 2.98342e-07, -1.69854e-07, 1.58299e-07, 1.79143e-07,
      -1.19884e-07, 1.32967e-08, 1.25766, 4.95172e-07, -2.47872e-07,
      3.82491e-07, -0.0371831, 0.0233571, -7.27197e-07, 0.0169894, -0.000153526,
      0.000302317, 9.50177e-08, -0.0331867, 9.71481e-05, 0.000100547,
      -5.20185e-06, -4.54307e-05, 0.00467098, 4.95172e-07, 1.34081, 2.88483e-07,
      2.04738e-07, 7.89636e-07, -6.29457e-08, -1.16343e-06, -0.000176982,
      -0.0169374, 0.00133415, 1.13961e-07, 5.83871e-05, 0.0313476, -0.0109009,
      0.00459359, 0.0008528, 1.35001e-05, -2.47872e-07, 2.88483e-07, 1.34083,
      -1.13626e-07, 6.27971e-07, -6.31283e-08, 5.83176e-08, -0.000289165,
      0.00133593, 0.0169342, -4.064e-09, 0.000126789, 0.0109008, 0.0313489,
      -0.000853556, 0.00459404, 4.37768e-05, 3.82491e-07, 2.04738e-07,
      -1.13626e-07, 1.34083, 1.29113e-07, -2.59622e-09, 0.00771153, 6.29926e-07,
      -4.78396e-07, 1.25256e-07, -0.00208028, -1.83644e-07, 1.24985e-07,
      -2.78769e-08, -3.75826e-07, -1.58829e-07, -6.10398e-07, -0.0371831,
      7.89636e-07, 6.27971e-07, 1.29113e-07, 1.93047, 0.0322435, -0.008621,
      2.41805e-09, 4.99241e-09, 1.45205e-09, -0.0204426, 3.06443e-08,
      -1.92004e-08, 2.5969e-09, -9.06493e-09, -6.78976e-09, -1.98658e-08,
      0.0233571, -6.29457e-08, -6.31283e-08, -2.59622e-09, 0.0322435, 19.4301;

  bool check_offdiag = offdiag_ref.isApprox(offdiag, 1e-5);
  if (!check_offdiag) {
    cout << "GW energies" << endl;
    cout << offdiag << endl;
    cout << "GW energies ref" << endl;
    cout << offdiag_ref << endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
