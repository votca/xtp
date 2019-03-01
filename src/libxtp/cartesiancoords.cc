#include <string>
#include <votca/xtp/cartesiancoords.h>
namespace votca {
namespace xtp {
CartesianCoords::CartesianCoords(const Orbitals& system)
    : CoordBase(CARTESIAN, system) {

  for (auto& qma : _qmMolecule) {
    auto pos = qma->getPos();
    _vector.push_back(pos.x());
    _vector.push_back(pos.y());
    _vector.push_back(pos.z());
  }

  _coords = Eigen::Map<Eigen::VectorXd>(_vector.data(), _vector.size());
}

}  // namespace xtp
}  // namespace votca
