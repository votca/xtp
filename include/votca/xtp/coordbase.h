#ifndef _VOTCA_XTP_COORDBASE
#define _VOTCA_XTP_COORDBASE
#include <string>
#include <votca/xtp/qmatom.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/eigen.h>

namespace votca { namespace xtp{

enum CoordType{CARTESIAN, INTERNAL};

class CoordBase{
public:
    Eigen::VectorXd Vector();
protected:
    // CoordBase cannot be instantiated
    CoordBase(const CoordType& _type, const std::vector<QMAtom*>& _qmm);
    CoordBase(const CoordType& , const Orbitals& orb);
    const CoordType type;
    const std::vector<QMAtom*>& qmMolecule;
    const int numAtoms;
    std::vector<double> vector;
};
} //xtp
} //votca
#endif
