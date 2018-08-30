#ifndef _VOTCA_XTP_INTERNAL_COORDS_H
#define _VOTCA_XTP_INTERNAL_COORDS_H
#include<vector>
#include<votca/xtp/qmatom.h>
#include<votca/xtp/coordbase.h>
#include<votca/xtp/cartesiancoords.h>
#include<votca/tools/elements.h>
#include<votca/xtp/eigen.h>
#include<votca/xtp/orbitals.h>

namespace votca { namespace xtp {

class InternalCoords: public CoordBase {
public:
    InternalCoords(const Orbitals& orb, const bool _withAuxiliary);
    InternalCoords(const std::vector<QMAtom*>& _qmm, const bool _withAuxiliary);

    InternalCoords(const Orbitals& orb);
    InternalCoords(const std::vector<QMAtom*>& _qmm);

    int getPossibleNumMols();

private:
    int numBonds;
    int numAngles;
    int numDihedrals;
    int possibleNumMols;

    const bool withAuxiliary;

    Eigen::MatrixXd bondMatrix;
};
} // xtp
} // votca

#endif // _VOTCA_XTP_INTERNAL_COORDS_H
