#ifndef _VOTCA_XTP_COORDBASE
#define _VOTCA_XTP_COORDBASE
#include <string>
#include <votca/xtp/eigen.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmmolecule.h>

namespace votca {
namespace xtp {

enum CoordType { CARTESIAN, INTERNAL };

class CoordBase {
 public:
  CoordBase();

  Eigen::VectorXd Vector();
  Eigen::VectorXd operator()();

  const Orbitals& System();

  void Increment(Eigen::VectorXd dx);
  int getNumAtoms();
  bool isApprox(CoordBase& other, const double& tol);

 protected:
  // CoordBase cannot be instantiated
  CoordBase(const CoordType& type, const Orbitals& system);

  const CoordType _type;
  const QMMolecule _qmMolecule;
  const int _numAtoms;
  const Orbitals& _system;
  std::vector<double> _vector;
  Eigen::VectorXd _coords;
};
typedef CoordBase CoordSystem;
}  // namespace xtp
}  // namespace votca

#endif
