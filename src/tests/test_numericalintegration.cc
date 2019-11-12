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

#define BOOST_TEST_MODULE aomatrix_test
#include "votca/xtp/orbitals.h"
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <votca/xtp/numerical_integrations.h>
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(numerical_integration_test)

AOBasis CreateBasis(const QMMolecule& mol) {
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

  BasisSet basis;
  basis.Load("3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  return aobasis;
}

Eigen::MatrixXd DMat() {
  Eigen::MatrixXd dmat = Eigen::MatrixXd::Zero(17, 17);
  dmat << 0.00157507, 0.0337454, 4.48905e-16, -5.93152e-16, 7.87133e-17,
      0.030876, 2.51254e-16, -1.49094e-16, 5.77899e-17, 0.00415998, -0.00445632,
      0.00415998, -0.00445632, 0.00415998, -0.00445632, 0.00415998, -0.00445632,
      0.0337454, 0.722983, 2.66427e-15, -4.44783e-15, 3.45846e-16, 0.661507,
      4.39854e-15, -2.02475e-15, 1.04832e-15, 0.0891262, -0.095475, 0.0891262,
      -0.095475, 0.0891262, -0.095475, 0.0891262, -0.095475, 4.48905e-16,
      2.66427e-15, 1.52199, 2.88658e-15, 2.09034e-15, -7.94212e-15, 0.215492,
      2.8727e-15, -1.40513e-15, 0.141933, -0.0402359, 0.141933, -0.0402359,
      -0.141933, 0.0402359, -0.141933, 0.0402359, -5.93152e-16, -4.44783e-15,
      2.88658e-15, 1.52199, -2.31759e-15, 9.21105e-15, -2.22045e-15, 0.215492,
      1.6263e-15, 0.141933, -0.0402359, -0.141933, 0.0402359, -0.141933,
      0.0402359, 0.141933, -0.0402359, 7.87133e-17, 3.45846e-16, 2.09034e-15,
      -2.31759e-15, 1.52199, 2.98902e-15, -2.04958e-15, 4.79738e-15, 0.215492,
      0.141933, -0.0402359, -0.141933, 0.0402359, 0.141933, -0.0402359,
      -0.141933, 0.0402359, 0.030876, 0.661507, -7.94212e-15, 9.21105e-15,
      2.98902e-15, 0.605259, 2.55488e-15, 2.7779e-17, 1.33759e-15, 0.0815477,
      -0.0873567, 0.0815477, -0.0873567, 0.0815477, -0.0873567, 0.0815477,
      -0.0873567, 2.51254e-16, 4.39854e-15, 0.215492, -2.22045e-15,
      -2.04958e-15, 2.55488e-15, 0.0305108, 3.29597e-17, -5.29036e-16,
      0.0200958, -0.00569686, 0.0200958, -0.00569686, -0.0200958, 0.00569686,
      -0.0200958, 0.00569686, -1.49094e-16, -2.02475e-15, 2.8727e-15, 0.215492,
      4.79738e-15, 2.7779e-17, 3.29597e-17, 0.0305108, 9.55941e-16, 0.0200958,
      -0.00569686, -0.0200958, 0.00569686, -0.0200958, 0.00569686, 0.0200958,
      -0.00569686, 5.77899e-17, 1.04832e-15, -1.40513e-15, 1.6263e-15, 0.215492,
      1.33759e-15, -5.29036e-16, 9.55941e-16, 0.0305108, 0.0200958, -0.00569686,
      -0.0200958, 0.00569686, 0.0200958, -0.00569686, -0.0200958, 0.00569686,
      0.00415998, 0.0891262, 0.141933, 0.141933, 0.141933, 0.0815477, 0.0200958,
      0.0200958, 0.0200958, 0.0506951, -0.0230264, -0.00224894, -0.00801753,
      -0.00224894, -0.00801753, -0.00224894, -0.00801753, -0.00445632,
      -0.095475, -0.0402359, -0.0402359, -0.0402359, -0.0873567, -0.00569686,
      -0.00569686, -0.00569686, -0.0230264, 0.0157992, -0.00801753, 0.0115445,
      -0.00801753, 0.0115445, -0.00801753, 0.0115445, 0.00415998, 0.0891262,
      0.141933, -0.141933, -0.141933, 0.0815477, 0.0200958, -0.0200958,
      -0.0200958, -0.00224894, -0.00801753, 0.0506951, -0.0230264, -0.00224894,
      -0.00801753, -0.00224894, -0.00801753, -0.00445632, -0.095475, -0.0402359,
      0.0402359, 0.0402359, -0.0873567, -0.00569686, 0.00569686, 0.00569686,
      -0.00801753, 0.0115445, -0.0230264, 0.0157992, -0.00801753, 0.0115445,
      -0.00801753, 0.0115445, 0.00415998, 0.0891262, -0.141933, -0.141933,
      0.141933, 0.0815477, -0.0200958, -0.0200958, 0.0200958, -0.00224894,
      -0.00801753, -0.00224894, -0.00801753, 0.0506951, -0.0230264, -0.00224894,
      -0.00801753, -0.00445632, -0.095475, 0.0402359, 0.0402359, -0.0402359,
      -0.0873567, 0.00569686, 0.00569686, -0.00569686, -0.00801753, 0.0115445,
      -0.00801753, 0.0115445, -0.0230264, 0.0157992, -0.00801753, 0.0115445,
      0.00415998, 0.0891262, -0.141933, 0.141933, -0.141933, 0.0815477,
      -0.0200958, 0.0200958, -0.0200958, -0.00224894, -0.00801753, -0.00224894,
      -0.00801753, -0.00224894, -0.00801753, 0.0506951, -0.0230264, -0.00445632,
      -0.095475, 0.0402359, -0.0402359, 0.0402359, -0.0873567, 0.00569686,
      -0.00569686, 0.00569686, -0.00801753, 0.0115445, -0.00801753, 0.0115445,
      -0.00801753, 0.0115445, -0.0230264, 0.0157992;
  return dmat;
}

BOOST_AUTO_TEST_CASE(vxc_test) {

  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  QMMolecule mol("none", 0);

  mol.LoadFromFile("molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);

  Eigen::MatrixXd dmat = DMat();

  NumericalIntegration num;
  num.GridSetup("medium", mol, aobasis);
  num.setXCfunctional("XC_GGA_X_PBE XC_GGA_C_PBE");

  BOOST_CHECK_EQUAL(num.getGridSize(), num.getGridpoints().size());
  BOOST_CHECK_EQUAL(num.getGridSize(), 53404);
  BOOST_CHECK_EQUAL(num.getBoxesSize(), 51);

  BOOST_CHECK_CLOSE(num.getExactExchange("XC_GGA_X_PBE XC_GGA_C_PBE"), 0.0,
                    1e-5);
  BOOST_CHECK_CLOSE(num.getExactExchange("XC_HYB_GGA_XC_PBEH"), 0.25, 1e-5);
  Mat_p_Energy e_vxc = num.IntegrateVXC(dmat);
  Eigen::MatrixXd vxc_ref = Eigen::MatrixXd::Zero(17, 17);
  vxc_ref << -0.604846, -0.193724, 1.4208e-12, 1.24779e-12, 1.33915e-12,
      -0.158347, 2.77358e-12, 2.74891e-12, 2.76197e-12, -0.0171771, -0.0712393,
      -0.0171771, -0.0712393, -0.0171771, -0.0712393, -0.0171771, -0.0712393,
      -0.193724, -0.830219, 3.50498e-13, 4.0751e-16, 1.87309e-13, -0.566479,
      1.38973e-13, 3.9268e-14, 9.4247e-14, -0.131012, -0.287372, -0.131012,
      -0.287372, -0.131012, -0.287372, -0.131012, -0.287372, 1.4208e-12,
      3.50498e-13, -0.814676, 5.25758e-14, 6.17368e-14, 2.75977e-13, -0.328319,
      6.87396e-14, 7.04428e-14, -0.10473, -0.0814687, -0.10473, -0.0814687,
      0.10473, 0.0814687, 0.10473, 0.0814687, 1.24779e-12, 4.0751e-16,
      5.25758e-14, -0.814676, 1.44859e-13, 5.42212e-14, 6.87022e-14, -0.328319,
      9.0481e-14, -0.10473, -0.0814687, 0.10473, 0.0814687, 0.10473, 0.0814687,
      -0.10473, -0.0814687, 1.33915e-12, 1.87309e-13, 6.17368e-14, 1.44859e-13,
      -0.814676, 1.76012e-13, 7.03983e-14, 9.05174e-14, -0.328319, -0.10473,
      -0.0814687, 0.10473, 0.0814687, -0.10473, -0.0814687, 0.10473, 0.0814687,
      -0.158347, -0.566479, 2.75977e-13, 5.42212e-14, 1.76012e-13, -0.508333,
      8.56063e-15, -6.62991e-14, -2.11309e-14, -0.157267, -0.288594, -0.157267,
      -0.288594, -0.157267, -0.288594, -0.157267, -0.288594, 2.77358e-12,
      1.38973e-13, -0.328319, 6.87022e-14, 7.03983e-14, 8.56063e-15, -0.308367,
      -4.39422e-14, -4.05135e-14, -0.113309, -0.0917206, -0.113309, -0.0917206,
      0.113309, 0.0917206, 0.113309, 0.0917206, 2.74891e-12, 3.9268e-14,
      6.87396e-14, -0.328319, 9.05174e-14, -6.62991e-14, -4.39422e-14,
      -0.308367, -3.9348e-14, -0.113309, -0.0917206, 0.113309, 0.0917206,
      0.113309, 0.0917206, -0.113309, -0.0917206, 2.76197e-12, 9.4247e-14,
      7.04428e-14, 9.0481e-14, -0.328319, -2.11309e-14, -4.05135e-14,
      -3.9348e-14, -0.308367, -0.113309, -0.0917206, 0.113309, 0.0917206,
      -0.113309, -0.0917206, 0.113309, 0.0917206, -0.0171771, -0.131012,
      -0.10473, -0.10473, -0.10473, -0.157267, -0.113309, -0.113309, -0.113309,
      -0.416577, -0.237903, -0.00491594, -0.0554951, -0.00491594, -0.0554951,
      -0.00491594, -0.0554951, -0.0712393, -0.287372, -0.0814687, -0.0814687,
      -0.0814687, -0.288594, -0.0917206, -0.0917206, -0.0917206, -0.237903,
      -0.2758, -0.0554951, -0.140797, -0.0554951, -0.140797, -0.0554951,
      -0.140797, -0.0171771, -0.131012, -0.10473, 0.10473, 0.10473, -0.157267,
      -0.113309, 0.113309, 0.113309, -0.00491594, -0.0554951, -0.416577,
      -0.237903, -0.00491594, -0.0554951, -0.00491594, -0.0554951, -0.0712393,
      -0.287372, -0.0814687, 0.0814687, 0.0814687, -0.288594, -0.0917206,
      0.0917206, 0.0917206, -0.0554951, -0.140797, -0.237903, -0.2758,
      -0.0554951, -0.140797, -0.0554951, -0.140797, -0.0171771, -0.131012,
      0.10473, 0.10473, -0.10473, -0.157267, 0.113309, 0.113309, -0.113309,
      -0.00491594, -0.0554951, -0.00491594, -0.0554951, -0.416577, -0.237903,
      -0.00491594, -0.0554951, -0.0712393, -0.287372, 0.0814687, 0.0814687,
      -0.0814687, -0.288594, 0.0917206, 0.0917206, -0.0917206, -0.0554951,
      -0.140797, -0.0554951, -0.140797, -0.237903, -0.2758, -0.0554951,
      -0.140797, -0.0171771, -0.131012, 0.10473, -0.10473, 0.10473, -0.157267,
      0.113309, -0.113309, 0.113309, -0.00491594, -0.0554951, -0.00491594,
      -0.0554951, -0.00491594, -0.0554951, -0.416577, -0.237903, -0.0712393,
      -0.287372, 0.0814687, -0.0814687, 0.0814687, -0.288594, 0.0917206,
      -0.0917206, 0.0917206, -0.0554951, -0.140797, -0.0554951, -0.140797,
      -0.0554951, -0.140797, -0.237903, -0.2758;
  bool check_vxc = e_vxc.matrix().isApprox(vxc_ref, 0.0001);

  BOOST_CHECK_CLOSE(e_vxc.energy(), -4.6303432151572643, 1e-5);
  if (!check_vxc) {
    std::cout << "ref" << std::endl;
    std::cout << vxc_ref << std::endl;
    std::cout << "calc" << std::endl;
    std::cout << e_vxc.matrix() << std::endl;
  }
  BOOST_CHECK_EQUAL(check_vxc, 1);
}

BOOST_AUTO_TEST_CASE(density_test) {

  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  QMMolecule mol("none", 0);

  mol.LoadFromFile("molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);

  Eigen::MatrixXd dmat = DMat();

  NumericalIntegration num;
  num.GridSetup("medium", mol, aobasis);

  double ntot = num.IntegrateDensity(dmat);
  BOOST_CHECK_CLOSE(ntot, 8.000000, 1e-5);

  Eigen::Vector3d pos = {3, 3, 3};

  BOOST_CHECK_CLOSE(num.IntegratePotential(pos), -1.543242, 1e-4);

  Eigen::Vector3d field = num.IntegrateField(pos);
  Eigen::Vector3d field_ref = {0.172802, 0.172802, 0.172802};
  bool field_check = field.isApprox(field_ref, 1e-5);
  if (!field_check) {
    std::cout << "field" << std::endl;
    std::cout << field.transpose() << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << field_ref.transpose() << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(gyration_test) {

  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            1.000000     1.000000     1.000000" << endl;
  xyzfile << " H            1.629118     1.629118     1.629118" << endl;
  xyzfile << " H           0.370882    0.370882     1.629118" << endl;
  xyzfile << " H            1.629118    0.370882    0.370882" << endl;
  xyzfile << " H           0.370882     1.629118   0.370882" << endl;
  xyzfile.close();

  QMMolecule mol("none", 0);

  mol.LoadFromFile("molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);

  Eigen::MatrixXd dmat = DMat();

  NumericalIntegration num;
  num.GridSetup("medium", mol, aobasis);

  Gyrationtensor tensor = num.IntegrateGyrationTensor(dmat);
  BOOST_CHECK_CLOSE(tensor.mass, 8.0000005, 1e-5);

  Eigen::Vector3d dip_ref = Eigen::Vector3d::Zero();
  dip_ref << 1.88973, 1.88973, 1.88973;
  bool centroid_check = dip_ref.isApprox(tensor.centroid, 1e-5);
  BOOST_CHECK_EQUAL(centroid_check, true);
  if (!centroid_check) {
    std::cout << "centroid" << std::endl;
    std::cout << tensor.centroid.transpose() << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << dip_ref.transpose() << std::endl;
  }
  Eigen::Matrix3d gyro_ref = Eigen::Matrix3d::Zero();
  gyro_ref << 0.596158, 2.85288e-12, 2.86873e-12, 2.85289e-12, 0.596158,
      2.87163e-12, 2.86874e-12, 2.87161e-12, 0.596158;
  bool gyro_check = gyro_ref.isApprox(tensor.gyration, 1e-5);
  BOOST_CHECK_EQUAL(gyro_check, true);
  if (!gyro_check) {
    std::cout << "gyro" << std::endl;
    std::cout << tensor.gyration << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << gyro_ref << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
