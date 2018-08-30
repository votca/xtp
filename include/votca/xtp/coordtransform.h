#ifndef _VOTCA_XTP_COORDINATE_TRANSFORM_H
#define _VOTCA_XTP_COORDINATE_TRANSFORM_H

#include <string>
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include <votca/xtp/coordbase.h>
#include <votca/xtp/cartesiancoords.h>
#include <votca/xtp/internalcoords.h>

namespace votca { namespace xtp {

using namespace boost::algorithm;


class CoordinateTransform{
public:
CoordinateTransform(CoordType _from, CoordType _to):
    conversion{std::make_pair(_from, _to)}{

    }

    /* CoordBase& operator()(CoordBase& coord){ */
    /*     if (conversion == "cartesian_to_internal"){ */
    /*         return InternalCoords(coord); */
    /*     } */
    /* } */

private:
    std::pair<CoordType, CoordType> conversion;
};
} // xtp
} // votca
#endif // _VOTCA_XTP_COORDINATE_TRANSFORM_H
