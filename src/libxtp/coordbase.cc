#include <votca/xtp/coordbase.h>

namespace votca {
namespace xtp {

CoordBase::CoordBase(const CoordType& type, const Orbitals& system)
    : _type(type),
      _qmMolecule(system.QMAtoms()),
      _numAtoms(_qmMolecule.size()),
      _system(system) {}

Eigen::VectorXd CoordBase::Vector() { return _coords; }

Eigen::VectorXd CoordBase::operator()() { return Vector(); }

const Orbitals& CoordBase::System() { return _system; }

void CoordBase::Increment(Eigen::VectorXd dx) {
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

int CoordBase::getNumAtoms() { return _numAtoms; }

bool CoordBase::isApprox(CoordBase& other, const double& tol) {
  return _coords.isApprox(other.Vector(), tol);
}

}  // namespace xtp
}  // namespace votca
