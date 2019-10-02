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

#include <string>
#include <votca/xtp/cartesiancoords.h>

namespace votca {
namespace xtp {
CartesianCoords::CartesianCoords(const QMMolecule& system)
    : CoordBase(CARTESIAN, system) {

  for (const QMAtom& qma : _qmMolecule) {
    const Eigen::Vector3d& pos = qma.getPos();
    _vector.push_back(pos.x());
    _vector.push_back(pos.y());
    _vector.push_back(pos.z());
  }

  _coords = Eigen::Map<Eigen::VectorXd>(_vector.data(), _vector.size());
}

}  // namespace xtp
}  // namespace votca
