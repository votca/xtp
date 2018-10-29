#ifndef _VOTCA_XTP_COORDINATE_TRANSFORM_H
#define _VOTCA_XTP_COORDINATE_TRANSFORM_H

#include <string>
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include <votca/xtp/coordbase.h>
#include <votca/xtp/cartesiancoords.h>
#include <votca/xtp/internalcoords.h>
#include <votca/xtp/orbitals.h>
#include <memory>

namespace votca { namespace xtp {

class CoordinateTransform{
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
} // xtp
} // votca
#endif // _VOTCA_XTP_COORDINATE_TRANSFORM_H
