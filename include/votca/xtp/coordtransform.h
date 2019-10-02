/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef _VOTCA_XTP_COORDINATE_TRANSFORM_H
#define _VOTCA_XTP_COORDINATE_TRANSFORM_H

#include <boost/algorithm/string.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <votca/xtp/cartesiancoords.h>
#include <votca/xtp/coordbase.h>
#include <votca/xtp/internalcoords.h>

namespace votca {
namespace xtp {

class CoordinateTransform {
 public:
  CoordinateTransform(CoordType from, CoordType to, const QMMolecule& system);

  void operator()(const CoordSystem& inCoord, CoordSystem& outCoord);

 private:
  std::pair<CoordType, CoordType> _conversion;
  Eigen::MatrixXd _jac;
  const QMMolecule& _system;
  std::unique_ptr<CoordSystem> _from;
  std::unique_ptr<CoordSystem> _to;
};
}  // namespace xtp
}  // namespace votca
#endif  // _VOTCA_XTP_COORDINATE_TRANSFORM_H
