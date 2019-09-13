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
#include <votca/xtp/ppm.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/sigma_ppm.h>
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

  Eigen::MatrixXd MOs = Eigen::MatrixXd::Zero(17, 17);
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

  Eigen::VectorXd mo_energy = Eigen::VectorXd::Zero(17);
  mo_energy << -0.612601, -0.341755, -0.341755, -0.341755, 0.137304, 0.16678,
      0.16678, 0.16678, 0.671592, 0.671592, 0.671592, 0.974255, 1.01205,
      1.01205, 1.01205, 1.64823, 19.4429;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);

  votca::ctp::Logger log;
  RPA rpa(log, Mmn);
  rpa.configure(4, 0, 16);
  rpa.setRPAInputEnergies(mo_energy);

  Sigma_PPM sigma = Sigma_PPM(Mmn, rpa);

  Sigma_PPM::options opt;
  opt.homo = 4;
  opt.qpmin = 0;
  opt.qpmax = 16;
  opt.rpamin = 0;
  opt.rpamax = 16;
  sigma.configure(opt);

  Eigen::MatrixXd x = sigma.CalcExchange();

  Eigen::MatrixXd x_ref = Eigen::MatrixXd::Zero(17, 17);
  x_ref << -0.898354, 5.71768e-07, -1.37292e-08, 1.98142e-07, -0.149915,
      2.16108e-06, -1.16525e-06, -3.17829e-07, -1.88839e-07, 2.21419e-07,
      -2.98344e-08, -0.0269881, -2.13239e-06, -3.09447e-06, -4.98945e-08,
      0.11249, -0.00337645, 5.71768e-07, -0.690315, 4.39806e-07, -3.6005e-07,
      -5.82737e-08, -0.0266566, -0.000334105, -0.000256846, -0.00109173,
      -0.00375978, 0.132508, 1.1112e-07, -0.0402288, 0.000418149, 0.000685308,
      -1.62313e-08, 4.06964e-09, -1.37292e-08, 4.39806e-07, -0.690316,
      1.49062e-07, 5.35669e-09, 0.000201816, -0.0244124, 0.0107079, -0.13183,
      -0.0138696, -0.00147968, -1.15888e-07, 0.000363237, 0.0401121,
      -0.00316551, 7.90554e-09, 3.08003e-09, 1.98142e-07, -3.6005e-07,
      1.49062e-07, -0.690315, -1.48218e-08, -0.000368948, 0.0107048, 0.0244113,
      -0.0139054, 0.131785, 0.00362461, 3.26075e-08, -0.000716146, -0.00315591,
      -0.0401082, -7.83667e-09, 1.75034e-09, -0.149915, -5.82737e-08,
      5.35669e-09, -1.48218e-08, -0.419561, 2.94693e-08, 1.4477e-08,
      2.02462e-10, 1.3573e-08, -1.31684e-08, 7.55841e-10, 0.111589,
      -1.43489e-08, 1.23778e-08, -9.2956e-09, 0.101203, -0.00966821,
      2.16108e-06, -0.0266566, 0.000201816, -0.000368948, 2.94693e-08,
      -0.162521, -2.50885e-08, 2.37541e-08, -6.18131e-05, -0.000396018,
      0.0286842, 7.29432e-08, 0.0853922, -0.000150722, -0.000326281,
      -1.86634e-07, -1.08635e-09, -1.16525e-06, -0.000334105, -0.0244124,
      0.0107048, 1.4477e-08, -2.50885e-08, -0.162521, -6.437e-08, -0.0249187,
      -0.01421, -0.000249755, -3.68182e-08, -0.000249614, -0.0806557,
      -0.0280463, 8.95564e-08, 1.72272e-09, -3.17829e-07, -0.000256846,
      0.0107079, 0.0244113, 2.02462e-10, 2.37541e-08, -6.437e-08, -0.162521,
      0.0142121, -0.0249168, -0.000313469, -4.22882e-08, -0.00025868, 0.0280471,
      -0.0806557, 1.90873e-08, 8.04122e-10, -1.88839e-07, -0.00109173, -0.13183,
      -0.0139054, 1.3573e-08, -6.18131e-05, -0.0249187, 0.0142121, -0.112664,
      5.15694e-09, 2.67503e-08, -6.82025e-08, -3.58992e-06, 0.00311123,
      -0.000577989, 1.57751e-08, 2.46073e-09, 2.21419e-07, -0.00375978,
      -0.0138696, 0.131785, -1.31684e-08, -0.000396018, -0.01421, -0.0249168,
      5.15694e-09, -0.112664, -3.26454e-09, 2.92286e-08, -3.07662e-05,
      0.000577593, 0.00311108, -4.33364e-09, -2.83729e-09, -2.98344e-08,
      0.132508, -0.00147968, 0.00362461, 7.55841e-10, 0.0286842, -0.000249755,
      -0.000313469, 2.67503e-08, -3.26454e-09, -0.112664, -2.55154e-08,
      0.00316414, 9.12634e-06, 2.94706e-05, 6.29731e-09, -3.74354e-09,
      -0.0269881, 1.1112e-07, -1.15888e-07, 3.26075e-08, 0.111589, 7.29432e-08,
      -3.68182e-08, -4.22882e-08, -6.82025e-08, 2.92286e-08, -2.55154e-08,
      -0.17259, -5.78498e-08, -9.45104e-08, 2.34917e-08, -0.0427377, 0.00896338,
      -2.13239e-06, -0.0402288, 0.000363237, -0.000716146, -1.43489e-08,
      0.0853922, -0.000249614, -0.00025868, -3.58992e-06, -3.07662e-05,
      0.00316414, -5.78498e-08, -0.131712, -1.06678e-08, 2.79713e-08,
      1.92221e-07, -7.78508e-10, -3.09447e-06, 0.000418149, 0.0401121,
      -0.00315591, 1.23778e-08, -0.000150722, -0.0806557, 0.0280471, 0.00311123,
      0.000577593, 9.12634e-06, -9.45104e-08, -1.06678e-08, -0.131712,
      2.65503e-07, 2.62621e-07, 2.69838e-10, -4.98945e-08, 0.000685308,
      -0.00316551, -0.0401082, -9.2956e-09, -0.000326281, -0.0280463,
      -0.0806557, -0.000577989, 0.00311108, 2.94706e-05, 2.34917e-08,
      2.79713e-08, 2.65503e-07, -0.131712, 5.71293e-09, -9.80113e-10, 0.11249,
      -1.62313e-08, 7.90554e-09, -7.83667e-09, 0.101203, -1.86634e-07,
      8.95564e-08, 1.90873e-08, 1.57751e-08, -4.33364e-09, 6.29731e-09,
      -0.0427377, 1.92221e-07, 2.62621e-07, 5.71293e-09, -0.10267, 0.0171175,
      -0.00337645, 4.06964e-09, 3.08003e-09, 1.75034e-09, -0.00966821,
      -1.08635e-09, 1.72272e-09, 8.04122e-10, 2.46073e-09, -2.83729e-09,
      -3.74354e-09, 0.00896338, -7.78508e-10, 2.69838e-10, -9.80113e-10,
      0.0171175, -0.0285864;

  bool check_x = x_ref.isApprox(x, 1e-5);
  if (!check_x) {
    cout << "Sigma X" << endl;
    cout << x << endl;
    cout << "Sigma X ref" << endl;
    cout << x_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_x, true);

  sigma.PrepareScreening();
  Eigen::VectorXd c_diag = sigma.CalcCorrelationDiag(mo_energy);
  Eigen::MatrixXd c_off = sigma.CalcCorrelationOffDiag(mo_energy);

  c_off.diagonal() = c_diag;
  Eigen::MatrixXd c_ref = Eigen::MatrixXd::Zero(17, 17);
  c_ref << 0.120676, 2.58689e-07, -2.52037e-07, 3.99968e-08, 0.0405292,
      -1.25428e-07, 5.37756e-08, 2.99233e-08, -8.10766e-08, -5.95507e-08,
      -1.4014e-07, -0.0233041, -3.41069e-07, -2.17655e-07, -2.87835e-08,
      -0.0147014, -0.00144565, 2.58689e-07, 0.0628008, -1.2626e-07, 3.62623e-08,
      -2.5785e-08, -7.87579e-05, -1.00701e-06, -7.48821e-07, -6.18615e-05,
      -0.000213051, 0.00750794, -1.27887e-07, 0.0193896, -0.000201918,
      -0.000330184, 3.4198e-07, -1.1817e-09, -2.52037e-07, -1.2626e-07,
      0.0628013, 9.54326e-08, -7.24038e-09, 5.49519e-07, -7.20056e-05,
      3.15477e-05, -0.00746948, -0.000785852, -8.38336e-05, 1.39395e-07,
      -0.000175035, -0.0193322, 0.00152515, -2.68216e-07, 2.67092e-09,
      3.99968e-08, 3.62623e-08, 9.54326e-08, 0.0628013, -5.97396e-09,
      -1.06723e-06, 3.1592e-05, 7.19606e-05, -0.000787889, 0.00746695,
      0.000205411, 1.02681e-08, 0.000345011, 0.00152207, 0.01933, 7.70523e-08,
      3.24027e-10, 0.0405292, -2.5785e-08, -7.24038e-09, -5.97396e-09,
      -0.0451823, -2.18631e-07, 1.28638e-07, 2.19836e-08, -5.17734e-08,
      -2.52021e-09, -4.23752e-08, 0.00254096, 7.35127e-08, 6.70113e-08,
      -2.26926e-09, 0.00843914, -0.00771886, -1.25428e-07, -7.87579e-05,
      5.49519e-07, -1.06723e-06, -2.18631e-07, -0.010061, 9.8763e-08,
      1.40982e-08, -1.52749e-05, -9.76693e-05, 0.0070754, 2.23156e-07,
      0.0520061, -9.16759e-05, -0.000198702, -7.91971e-08, 3.73239e-08,
      5.37756e-08, -1.00701e-06, -7.20056e-05, 3.1592e-05, 1.28638e-07,
      9.8763e-08, -0.0100613, -6.99128e-08, -0.00614656, -0.00350514,
      -6.16153e-05, -1.81029e-07, -0.000152159, -0.0491215, -0.017081,
      8.98441e-08, -2.4424e-08, 2.99233e-08, -7.48821e-07, 3.15477e-05,
      7.19606e-05, 2.19836e-08, 1.40982e-08, -6.99128e-08, -0.0100616,
      0.00350559, -0.00614613, -7.73261e-05, 5.52026e-08, -0.000157503,
      0.0170815, -0.049122, -4.06344e-08, 2.22078e-09, -8.10766e-08,
      -6.18615e-05, -0.00746948, -0.000787889, -5.17734e-08, -1.52749e-05,
      -0.00614656, 0.00350559, -0.0453151, 3.12997e-08, -3.5253e-09,
      2.21917e-07, 2.0425e-05, -0.0179271, 0.00332992, -1.20113e-07,
      -1.13638e-08, -5.95507e-08, -0.000213051, -0.000785852, 0.00746695,
      -2.52021e-09, -9.76693e-05, -0.00350514, -0.00614613, 3.12997e-08,
      -0.0453153, -4.16704e-08, 8.99497e-08, 0.000177348, -0.00332983,
      -0.0179262, -5.24173e-08, -4.27551e-09, -1.4014e-07, 0.00750794,
      -8.38336e-05, 0.000205411, -4.23752e-08, 0.0070754, -6.16153e-05,
      -7.73261e-05, -3.5253e-09, -4.16704e-08, -0.0453152, 3.19659e-07,
      -0.0182333, -5.22712e-05, -0.000170735, -2.15984e-07, -1.51435e-08,
      -0.0233041, -1.27887e-07, 1.39395e-07, 1.02681e-08, 0.00254096,
      2.23156e-07, -1.81029e-07, 5.52026e-08, 2.21917e-07, 8.99497e-08,
      3.19659e-07, -0.0937988, -5.91266e-07, -3.95018e-07, -4.00936e-09,
      -0.0251495, 0.0103057, -3.41069e-07, 0.0193896, -0.000175035, 0.000345011,
      7.35127e-08, 0.0520061, -0.000152159, -0.000157503, 2.0425e-05,
      0.000177348, -0.0182333, -5.91266e-07, -0.159212, -3.1902e-07,
      2.58648e-10, 2.82618e-07, -8.13748e-08, -2.17655e-07, -0.000201918,
      -0.0193322, 0.00152207, 6.70113e-08, -9.16759e-05, -0.0491215, 0.0170815,
      -0.0179271, -0.00332983, -5.22712e-05, -3.95018e-07, -3.1902e-07,
      -0.159212, 2.89916e-07, 3.1149e-07, -8.13703e-08, -2.87835e-08,
      -0.000330184, 0.00152515, 0.01933, -2.26926e-09, -0.000198702, -0.017081,
      -0.049122, 0.00332992, -0.0179262, -0.000170735, -4.00936e-09,
      2.58648e-10, 2.89916e-07, -0.159213, 1.48831e-08, -3.54064e-09,
      -0.0147014, 3.4198e-07, -2.68216e-07, 7.70523e-08, 0.00843914,
      -7.91971e-08, 8.98441e-08, -4.06344e-08, -1.20113e-07, -5.24173e-08,
      -2.15984e-07, -0.0251495, 2.82618e-07, 3.1149e-07, 1.48831e-08,
      -0.0229993, 0.012047, -0.00144565, -1.1817e-09, 2.67092e-09, 3.24027e-10,
      -0.00771886, 3.73239e-08, -2.4424e-08, 2.22078e-09, -1.13638e-08,
      -4.27551e-09, -1.51435e-08, 0.0103057, -8.13748e-08, -8.13703e-08,
      -3.54064e-09, 0.012047, -0.40848;

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
