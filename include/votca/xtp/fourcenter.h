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
#ifndef VOTCA_XTP_FOURCENTER_H
#define VOTCA_XTP_FOURCENTER_H
#include <votca/xtp/eigen.h>
#include <votca/xtp/multiarray.h>

/**
 * \brief Calculates four center electron overlap integrals for DFT.
 *
 *
 *
 */

namespace votca {
namespace xtp {

class AOBasis;
class AOShell;
class FCMatrix {

 public:
  void Fill_4c_small_molecule(const AOBasis& dftbasis);

  const Eigen::VectorXd& get_4c_vector() const { return _4c_vector; }

  bool FillFourCenterRepBlock(tensor4d& block, const AOShell& shell_1,
                              const AOShell& shell_2, const AOShell& shell_3,
                              const AOShell& shell_4) const;

 private:
  Eigen::VectorXd _4c_vector;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_FOURCENTER_H
