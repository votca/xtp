#include<votca/xtp/coordbase.h>

namespace votca { namespace xtp{

CoordBase::CoordBase(const CoordType& _type, const Orbitals& _orb):
    CoordBase(_type, _orb.QMAtoms()){};

CoordBase::CoordBase(const CoordType& _type, const std::vector<QMAtom*>& _qmm):
    type(_type), qmMolecule(_qmm), numAtoms(_qmm.size()){}

Eigen::VectorXd CoordBase::Vector(){
    Eigen::VectorXd temp = Eigen::Map<Eigen::VectorXd>(vector.data(), vector.size());
    return temp;
}

}//xtp
}//votca
