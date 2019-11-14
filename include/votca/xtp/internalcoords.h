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
#ifndef _VOTCA_XTP_INTERNAL_COORDS_H
#define _VOTCA_XTP_INTERNAL_COORDS_H
#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <vector>
#include <votca/tools/elements.h>
#include <votca/xtp/CoordContainer.h>
#include <votca/xtp/cartesiancoords.h>
#include <votca/xtp/coordbase.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmatom.h>

namespace votca {
namespace xtp {

class InternalCoords : public CoordBase {
 public:
  InternalCoords(const QMMolecule& mol, bool withAuxiliary);
  InternalCoords(const QMMolecule& mol);

  Index getNumBonds() const;
  Index getNumHBonds() const;
  Index getNumAngles() const;
  Index getNumAuxBonds() const;
  Index getNumDihedrals() const;
  const Eigen::MatrixXd& getWilsonBMatrix() const;

  Eigen::MatrixXd CalculatePseudoInverseB() const {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    // svd.setThreshold(1e-9);
    svd.compute(_wilsonBMatrix, Eigen::DecompositionOptions::ComputeThinU |
                                    Eigen::DecompositionOptions::ComputeThinV);
    return svd.matrixV() *
           svd.singularValues().head(svd.nonzeroSingularValues()) *
           svd.matrixU().transpose();
  }

  Eigen::VectorXd CalcInternalForces(
      const Eigen::VectorXd& cartesianForces) const {
    return CalculateP() * CalculatePseudoInverseB().transpose() *
           cartesianForces;
  }

  Eigen::MatrixXd CalculateP() const {
    return _wilsonBMatrix * CalculatePseudoInverseB();
  }

  friend std::ostream& operator<<(std::ostream& s, const InternalCoords& ic);

 private:
  // Bonds will be stored in an undirected graph
  using BglGraph =
      boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;

  Index _numBonds = 0;
  Index _numInterMolBonds = 0;  // interFragmentBonds
  Index _numHBonds = 0;
  Index _numAngles = 0;
  Index _numAuxBonds = 0;
  Index _numDihedrals = 0;

  const bool _withAuxiliary;

  Eigen::MatrixXd _wilsonBMatrix;
  BglGraph _bondGraph;

  CoordContainer<BondIdx, double> _bonds;
  CoordContainer<AngleIdx, double> _angles;
  CoordContainer<DihedralIdx, double> _dihedrals;

  CoordContainer<BondIdx, double> _auxBonds;

  CartesianCoords _cartCoords;

  void ConnectBonds();
  void ConnectMolecules();
  void ConnectHBonds();
  void CalculateAnglesDihedrals();
  void PopulateWilsonMatrix();
};

}  // namespace xtp
}  // namespace votca

#endif  // _VOTCA_XTP_INTERNAL_COORDS_H
