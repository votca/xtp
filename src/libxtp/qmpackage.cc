/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include "votca/xtp/qmpackage.h"
#include <boost/algorithm/string.hpp>
#include <votca/xtp/ecpaobasis.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {
using std::flush;

void QMPackage::ParseCommonOptions(tools::Property& options) {

  std::string key = "package";
  std::string name = options.get(key + ".name").as<std::string>();

  if (name != getPackageName()) {
    throw std::runtime_error("Tried to use " + name +
                             " package. Wrong options file");
  }

  if (getPackageName() != "xtp") {
    _executable = options.ifExistsReturnElseThrowRuntimeError<std::string>(
        key + ".executable");
    _memory = options.ifExistsReturnElseThrowRuntimeError<std::string>(
        key + ".memory");
    _options = options.ifExistsReturnElseThrowRuntimeError<std::string>(
        key + ".options");
    _scratch_dir = options.ifExistsReturnElseThrowRuntimeError<std::string>(
        key + ".scratch");
  }

  _charge = options.ifExistsReturnElseThrowRuntimeError<Index>(key + ".charge");
  _spin = options.ifExistsReturnElseThrowRuntimeError<Index>(key + ".spin");
  _cleanup =
      options.ifExistsReturnElseReturnDefault(key + ".cleanup", _cleanup);
  _dpl_spacing = options.ifExistsReturnElseReturnDefault(
      key + ".dipole_spacing", _dpl_spacing);

  _write_guess = options.ifExistsReturnElseReturnDefault<bool>(
      key + ".read_guess", _write_guess);

  if (options.exists(key + ".basisset")) {
    _basisset_name = options.get(key + ".basisset").as<std::string>();
    _write_basis_set = true;
  }
  if (options.exists(key + ".auxbasisset")) {
    _auxbasisset_name = options.get(key + ".auxbasisset").as<std::string>();
    _write_auxbasis_set = true;
  }

  if (options.exists(key + ".ecp")) {
    _write_pseudopotentials = true;
    _ecp_name = options.get(key + ".ecp").as<std::string>();
  }
}

void QMPackage::ReorderOutput(Orbitals& orbitals) const {
  if (!orbitals.hasQMAtoms()) {
    throw std::runtime_error("Orbitals object has no QMAtoms");
  }

  AOBasis dftbasis = orbitals.SetupDftBasis();
  // necessary to update nuclear charges on qmatoms
  if (orbitals.hasECPName()) {
    ECPBasisSet ecps;
    ecps.Load(orbitals.getECPName());
    ECPAOBasis ecpbasis;
    ecpbasis.Fill(ecps, orbitals.QMAtoms());
  }

  if (orbitals.hasMOs()) {
    dftbasis.ReorderMOs(orbitals.MOs().eigenvectors(), getPackageName(), "xtp");
    XTP_LOG(logDEBUG, *_pLog) << "Reordered MOs" << flush;
  }

  return;
}

Eigen::MatrixXd QMPackage::ReorderMOsBack(const Orbitals& orbitals) const {
  BasisSet dftbasisset;
  dftbasisset.Load(_basisset_name);
  if (!orbitals.hasQMAtoms()) {
    throw std::runtime_error("Orbitals object has no QMAtoms");
  }
  AOBasis dftbasis;
  dftbasis.Fill(dftbasisset, orbitals.QMAtoms());
  Eigen::MatrixXd result = orbitals.MOs().eigenvectors();
  dftbasis.ReorderMOs(result, "xtp", getPackageName());
  return result;
}

std::vector<QMPackage::MinimalMMCharge> QMPackage::SplitMultipoles(
    const StaticSite& aps) const {

  std::vector<QMPackage::MinimalMMCharge> multipoles_split;
  // Calculate virtual charge positions
  double a = _dpl_spacing;                // this is in a0
  double mag_d = aps.getDipole().norm();  // this is in e * a0
  const Eigen::Vector3d dir_d = aps.getDipole().normalized();
  const Eigen::Vector3d A = aps.getPos() + 0.5 * a * dir_d;
  const Eigen::Vector3d B = aps.getPos() - 0.5 * a * dir_d;
  double qA = mag_d / a;
  double qB = -qA;
  if (std::abs(qA) > 1e-12) {
    multipoles_split.push_back(MinimalMMCharge(A, qA));
    multipoles_split.push_back(MinimalMMCharge(B, qB));
  }

  if (aps.getRank() > 1) {
    const Eigen::Matrix3d components = aps.CalculateCartesianMultipole();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.computeDirect(components);
    double a2 = 2 * _dpl_spacing;
    for (Index i = 0; i < 3; i++) {
      double q = es.eigenvalues()[i] / (a2 * a2);
      const Eigen::Vector3d vec1 =
          aps.getPos() + 0.5 * a2 * es.eigenvectors().col(i);
      const Eigen::Vector3d vec2 =
          aps.getPos() - 0.5 * a2 * es.eigenvectors().col(i);
      multipoles_split.push_back(MinimalMMCharge(vec1, q));
      multipoles_split.push_back(MinimalMMCharge(vec2, q));
    }
  }
  return multipoles_split;
}

std::vector<std::string> QMPackage::GetLineAndSplit(
    std::ifstream& input_file, const std::string separators) const {
  std::string line;
  getline(input_file, line);
  boost::trim(line);
  tools::Tokenizer tok(line, separators.c_str());
  return tok.ToVector();
}

}  // namespace xtp
}  // namespace votca
