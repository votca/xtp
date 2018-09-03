#ifndef _VOTCA_XTP_INTERNAL_COORDS_H
#define _VOTCA_XTP_INTERNAL_COORDS_H
#include<vector>
#include<votca/xtp/qmatom.h>
#include<votca/xtp/coordbase.h>
#include<votca/xtp/cartesiancoords.h>
#include<votca/tools/elements.h>
#include<votca/xtp/eigen.h>
#include<votca/xtp/orbitals.h>
#include<boost/graph/adjacency_list.hpp>

namespace votca { namespace xtp {

//Bonds will be stored in an undirected graph
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> BglGraph;

class InternalCoords: public CoordBase {
public:
    InternalCoords(const Orbitals& orb, const bool _withAuxiliary);
    InternalCoords(const std::vector<QMAtom*>& _qmm, const bool _withAuxiliary);

    InternalCoords(const Orbitals& orb);
    InternalCoords(const std::vector<QMAtom*>& _qmm);

    int getPossibleNumMols();
    int getNumBonds();
    int getNumHBonds();

private:
    int numBonds;
    int numInterMolBonds; // interFragmentBonds
    int numHBonds;
    int numAngles;
    int numDihedrals;
    int possibleNumMols;

    const bool withAuxiliary;

    Eigen::MatrixXd bondMatrix;
    BglGraph bondGraph;


    void ConnectBonds();
    void ConnectMolecules();
    void ConnectHBonds();
};
} // xtp
} // votca

#endif // _VOTCA_XTP_INTERNAL_COORDS_H
