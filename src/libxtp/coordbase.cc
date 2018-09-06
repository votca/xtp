#include<votca/xtp/coordbase.h>

namespace votca { namespace xtp{

CoordBase::CoordBase(const CoordType& type, const Orbitals& system):
    _type(type),
    _qmMolecule(system.QMAtoms()),
    _numAtoms(_qmMolecule.size()),
    _system(system){}



Eigen::VectorXd CoordBase::Vector(){
    Eigen::VectorXd temp = Eigen::Map<Eigen::VectorXd>(_vector.data(), _vector.size());
    return temp;
}

Eigen::VectorXd CoordBase::operator()(){
    return Vector();
}

const Orbitals& CoordBase::System(){
    return _system;
}

}//xtp
}//votca
