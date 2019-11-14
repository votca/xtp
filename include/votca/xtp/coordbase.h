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

#pragma once
#ifndef _VOTCA_XTP_COORDBASE
#define _VOTCA_XTP_COORDBASE
#include <string>
#include <votca/xtp/eigen.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmmolecule.h>

namespace votca {
namespace xtp {

enum CoordType { CARTESIAN, INTERNAL };

class CoordBase {
 public:
  CoordBase();

  const Eigen::VectorXd& Vector() const;
  const Eigen::VectorXd& operator()() const;

  void Increment(const Eigen::VectorXd& dx);
  Index getNumAtoms() const;
  bool isApprox(CoordBase& other, const double& tol) const;

 protected:
  // CoordBase cannot be instantiated
  CoordBase(const CoordType& type, const QMMolecule& system);

  const CoordType _type;
  const QMMolecule& _qmMolecule;
  const Index _numAtoms;
  std::vector<double> _vector;
  Eigen::VectorXd _coords;
};
using CoordSystem = CoordBase;
}  // namespace xtp
}  // namespace votca

#endif
