/*
 *            Copyright 2009-2020 The VOTCA Development Team
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
#ifndef VOTCA_XTP_AITKEN_H
#define VOTCA_XTP_AITKEN_H

// VOTCA includes
#include <votca/tools/types.h>

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief Aitken Extrapolation Delta Squared Process
 *
 * https://en.wikipedia.org/wiki/Aitken%27s_delta-squared_process
 *
 */
template <class Func>
class Aitken {
 public:
  enum Errors { success, smalldenom, notconverged };
  Aitken(Index max_iterations, double tolerance)
      : _max_iterations(max_iterations), _tolerance(tolerance) {}

  std::vector<double> FindRoot(const Func& f, double x0) {
    _info = Errors::notconverged;
    double x = x0;
    double derivative =0.0;
    for (_iter = 0; _iter < _max_iterations; _iter++) {

      double x1 = f.value(x)+x;     // this is solving f(x)=x, but f.value is f(x)-x
      double x2 = f.value(x1)+x1;
      double denominator = (x2 - x1) - (x1 - x);

      if (std::abs(denominator) < 1e-12) {
        _info = Errors::smalldenom;
        break;
      }

      double aitkenX = x2 - ( std::pow(x2 - x1,2) )/denominator;

      if (std::abs(aitkenX - x2) < _tolerance) {
        _info = Errors::success;
        break;
      }

      x = aitkenX;

      std::cout << " Aitken " << x << std::endl;
    }

    std::vector<double> fixedpoint(2);
    fixedpoint[0] = x;
    fixedpoint[1] = -1.0/f.deriv(x);

    return fixedpoint;
  }

  Errors getInfo() const { return _info; }
  Index getIterations() const { return _iter; }

 private:
  Errors _info = Errors::notconverged;
  Index _max_iterations;
  Index _iter;
  double _tolerance;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_AITKEN_H
