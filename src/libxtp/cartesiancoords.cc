#include<votca/xtp/cartesiancoords.h>

namespace votca { namespace xtp{
CartesianCoords::CartesianCoords(const std::vector<QMAtom*>& _qmm):
    CoordBase(CARTESIAN, _qmm){

    for (auto& qma:qmMolecule){
        auto pos = qma->getPos();
        vector.push_back(pos.x());
        vector.push_back(pos.y());
        vector.push_back(pos.z());
    }
}

CartesianCoords::CartesianCoords(const Orbitals& orb):
    CartesianCoords(orb.QMAtoms()){}

} // xtp
} // votca
