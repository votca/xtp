#ifndef _VOTCA_XTP_CARTESIAN_COORD
#define _VOTCA_XTP_CARTESIAN_COORD

#include <votca/xtp/coordbase.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmatom.h>

namespace votca {
namespace xtp {

class CartesianCoords : public CoordBase {
 public:
  CartesianCoords(const Orbitals& orb);
};

}  // namespace xtp
}  // namespace votca
#endif  // _VOTCA_XTP_CARTESIAN_COORD
