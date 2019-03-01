#ifndef _VOTCA_XTP_COORDINATE_TRANSFORM_H
#define _VOTCA_XTP_COORDINATE_TRANSFORM_H

#include <boost/algorithm/string.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <votca/xtp/cartesiancoords.h>
#include <votca/xtp/coordbase.h>
#include <votca/xtp/internalcoords.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {

class CoordinateTransform {
 public:
  CoordinateTransform(CoordType from, CoordType to, const Orbitals& system);

  void operator()(const CoordSystem& inCoord, CoordSystem& outCoord);

 private:
  std::pair<CoordType, CoordType> _conversion;
  Eigen::MatrixXd _jac;
  const Orbitals& _system;
  std::unique_ptr<CoordSystem> _from;
  std::unique_ptr<CoordSystem> _to;
};
}  // namespace xtp
}  // namespace votca
#endif  // _VOTCA_XTP_COORDINATE_TRANSFORM_H
