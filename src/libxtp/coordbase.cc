#include <votca/xtp/coordbase.h>

namespace votca {
namespace xtp {

CoordBase::CoordBase(const CoordType& type, const Orbitals& system)
    : _type(type),
      _qmMolecule(system.QMAtoms()),
      _numAtoms(_qmMolecule.size()),
      _system(system) {}

const Eigen::VectorXd& CoordBase::Vector() const { return _coords; }

const Eigen::VectorXd& CoordBase::operator()() const { return Vector(); }

const Orbitals& CoordBase::System() const { return _system; }

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
