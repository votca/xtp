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
#ifndef VOTCA_XTP_DIIS_H
#define VOTCA_XTP_DIIS_H

#include <vector>
#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {

class DIIS {
 public:
  void Update(Index maxerrorindex, const Eigen::MatrixXd& errormatrix);
  Eigen::VectorXd CalcCoeff();

  void setHistLength(Index length) { _histlength = length; }

  bool Info() { return success; }

 private:
  bool success = true;
  Index _histlength;
  std::vector<std::vector<double> > _Diis_Bs;
  std::vector<Eigen::MatrixXd> _errormatrixhist;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DIIS_H
