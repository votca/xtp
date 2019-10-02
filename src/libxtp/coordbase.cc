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

#include <votca/xtp/coordbase.h>

namespace votca {
namespace xtp {

CoordBase::CoordBase(const CoordType& type, const QMMolecule& system)
    : _type(type), _qmMolecule(system), _numAtoms(_qmMolecule.size()) {}

const Eigen::VectorXd& CoordBase::Vector() const { return _coords; }

const Eigen::VectorXd& CoordBase::operator()() const { return Vector(); }

void CoordBase::Increment(const Eigen::VectorXd& dx) {
  if (dx.size() != _coords.size()) {
    std::ostringstream stream;
    stream << "Dimensions do not match." << std::endl
           << "I am a " << _coords.size() << " dimensional vector." << std::endl
           << "You gave me a " << dx.size() << " dimensional vector."
           << std::endl;
    throw std::runtime_error(stream.str());
  }

  _coords += dx;
}

int CoordBase::getNumAtoms() const { return _numAtoms; }

bool CoordBase::isApprox(CoordBase& other, const double& tol) const {
  return _coords.isApprox(other.Vector(), tol);
}

}  // namespace xtp
}  // namespace votca
