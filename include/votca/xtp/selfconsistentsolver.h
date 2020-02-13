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
#ifndef XTP_VXC_SELFCONSISTENTSOLVER_H
#define XTP_VXC_SELFCONSISTENTSOLVER_H

#include <votca/tools/types.h>
#include <votca/xtp/eigen.h>
namespace votca {
namespace xtp {

/**
 * \brief Simple self-consistent fixed point solver for self-energy
 */
template <class Func>
class Selfconsistentsolver {
 public:
  enum Errors { success, notconverged };
  Selfconsistentsolver(Index max_iterations, double tolerance)
      : _max_iterations(max_iterations), _tolerance(tolerance) {}

  Selfconsistentsolver(Index max_iterations, double tolerance, double alpha)
      : _max_iterations(max_iterations), _tolerance(tolerance), _alpha(alpha) {}

  double FindRoot(const Func& f, double x0) {
    _info = Errors::notconverged;
    double x = x0;
    for (_iter = 0; _iter < _max_iterations; _iter++) {
      std::pair<double, double> res = f(x);
      double step = std::abs(res.first - x0);
      if (std::abs(step) < _tolerance) {
        _info = Errors::success;
        break;
      }
    }
    return x;
  }

  Errors getInfo() const { return _info; }
  Index getIterations() const { return _iter; }

 private:
  Errors _info = Errors::notconverged;
  Index _max_iterations;
  Index _iter;
  double _tolerance;
  double _alpha = 1.0;
};

}  // namespace xtp
}  // namespace votca
#endif  // XTP_VXC_SELFCONSISTENTSOLVER_H
