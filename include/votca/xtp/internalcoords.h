#ifndef _VOTCA_XTP_INTERNAL_COORDS_H
#define _VOTCA_XTP_INTERNAL_COORDS_H
#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <vector>
#include <votca/tools/elements.h>
#include <votca/xtp/CoordContainer.h>
#include <votca/xtp/cartesiancoords.h>
#include <votca/xtp/coordbase.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmatom.h>

namespace votca {
namespace xtp {

// Bonds will be stored in an undirected graph
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
    BglGraph;

class InternalCoords : public CoordBase {
 public:
  InternalCoords(const Orbitals& orb, const bool _withAuxiliary);
  InternalCoords(const Orbitals& orb);

  int getPossibleNumMols();
  int getNumBonds();
  int getNumHBonds();
  int getNumAngles();
  int getNumAuxBonds();
  int getNumDihedrals();
  using CoordBase::Increment;
  void Increment(Eigen::VectorXd dx);
  Eigen::MatrixXd getWilsonBMatrix();

  friend std::ostream& operator<<(std::ostream& s, const InternalCoords& ic);

 private:
  int _numBonds;
  int _numInterMolBonds;  // interFragmentBonds
  int _numHBonds;
  int _numAngles;
  int _numAuxBonds;
  int _numDihedrals;
  int _possibleNumMols;

  const bool _withAuxiliary;

  Eigen::MatrixXd _wilsonBMatrix;
  BglGraph _bondGraph;

  /* std::map< std::tuple<int, int>, double> _bonds; */
  /* std::map< std::tuple<int, int, int> , double> _angles; */
  /* std::map< std::tuple<int, int, int, int>, double> _dihedrals; */

  CoordContainer<BondIdx, double> _bonds;
  CoordContainer<AngleIdx, double> _angles;
  CoordContainer<DihedralIdx, double> _dihedrals;

  /* std::vector<std::tuple<int, int>> _bondInsOrder; */
  /* std::vector<std::tuple<int, int, int>> _angleInsOrder; */
  /* std::vector<std::tuple<int, int, int,int >> _dihedralInsOrder; */

  CoordContainer<BondIdx, double> _auxBonds;

  CartesianCoords _cartCoords;

  void ConnectBonds();
  void ConnectMolecules();
  void ConnectHBonds();
  void CalculateAnglesDihedrals();
  void PopulateWilsonMatrix();
};

}  // namespace xtp
}  // namespace votca

#endif  // _VOTCA_XTP_INTERNAL_COORDS_H
