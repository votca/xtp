/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Local VOTCA includes
#include "votca/xtp/orbreorder.h"

namespace votca {
namespace xtp {

std::vector<std::array<Index, 2>> OrbReorder::getTranspositions(
    Index nrOfFunctions) {
  switch (nrOfFunctions) {
    case 1:
      return _transpositions.sOrb;
    case 3:
      return _transpositions.pOrb;
    case 5:
      return _transpositions.dOrb;
    case 7:
      return _transpositions.fOrb;
    case 9:
      return _transpositions.gOrb;
    default:
      throw std::runtime_error("Impossible number of functions in shell");
  }
}

void OrbReorder::reorderOrbitals(Eigen::MatrixXd& moCoefficients,
                                 AOBasis& basis) {
  std::vector<Index> multiplier;
  multiplier.reserve(basis.AOBasisSize());

  // reorder and get multiplier vector
  Index currentFunction = 0;
  Index nrOfFunctions = 0;
  for (const AOShell& shell : basis) {
    // Make multiplier vector
    std::vector<Index> shellmultiplier{
        _multipliers.begin() + shell.getOffset(),
        _multipliers.begin() + shell.getOffset() + shell.getNumFunc()};
    multiplier.insert(multiplier.end(), shellmultiplier.begin(),
                      shellmultiplier.end());
    // reorder
    nrOfFunctions = shell.getNumFunc();
    for (auto& transposition : getTranspositions(nrOfFunctions)) {
      moCoefficients.row(currentFunction + transposition[0])
          .swap(moCoefficients.row(currentFunction + transposition[1]));
    }
    currentFunction += nrOfFunctions;
  }

  // Multiply by multiplier
  for (Index i = 0; i < moCoefficients.cols(); i++) {
    moCoefficients.row(i) = multiplier[i] * moCoefficients.row(i);
  }
}

}  // namespace xtp
}  // namespace votca
