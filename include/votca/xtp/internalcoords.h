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
  InternalCoords(const Orbitals& orb, bool _withAuxiliary);
  InternalCoords(const Orbitals& orb);

  int getNumBonds() const;
  int getNumHBonds() const;
  int getNumAngles() const;
  int getNumAuxBonds() const;
  int getNumDihedrals() const;
  const Eigen::MatrixXd& getWilsonBMatrix() const;

  friend std::ostream& operator<<(std::ostream& s, const InternalCoords& ic);

 private:
  int _numBonds = 0;
  int _numInterMolBonds = 0;  // interFragmentBonds
  int _numHBonds = 0;
  int _numAngles = 0;
  int _numAuxBonds = 0;
  int _numDihedrals = 0;

  const bool _withAuxiliary;

  Eigen::MatrixXd _wilsonBMatrix;
  BglGraph _bondGraph;

  CoordContainer<BondIdx, double> _bonds;
  CoordContainer<AngleIdx, double> _angles;
  CoordContainer<DihedralIdx, double> _dihedrals;

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
