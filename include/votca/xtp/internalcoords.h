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
    InternalCoords(const Orbitals& orb);

    int getPossibleNumMols();
    int getNumBonds();
    int getNumHBonds();
    int getNumAngles();
    int getNumAuxBonds();
    int getNumDihedrals();

private:
    int _numBonds;
    int _numInterMolBonds; // interFragmentBonds
    int _numHBonds;
    int _numAngles;
    int _numAuxBonds;
    int _numDihedrals;
    int _possibleNumMols;

    const bool _withAuxiliary;

    Eigen::MatrixXd _bondMatrix;
    Eigen::MatrixXd _wilsonBMatrix;
    BglGraph _bondGraph;


    std::map< std::tuple<int, int>, double> _bonds;
    std::map< std::tuple<int, int, int> , double> _angles;
    std::map< std::tuple<int, int, int, int>, double> _dihedrals;

    std::vector< std::pair<int, int> > _auxBonds;

    void ConnectBonds();
    void ConnectMolecules();
    void ConnectHBonds();
    void CalculateAnglesDihedrals();
    void IsAuxBond();
};
} // xtp
} // votca

#endif // _VOTCA_XTP_INTERNAL_COORDS_H
