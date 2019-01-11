/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_RPA_EIGENSOLUTION_H
#define _VOTCA_XTP_RPA_EIGENSOLUTION_H

namespace votca
{
namespace xtp
{

struct rpa_eigensolution {
    Eigen::VectorXd _Omega; // Eigenvalues
    Eigen::MatrixXd _XpY; // Eigenvector components (X + Y)
};

}
}

#endif /* _VOTCA_XTP_RPA_EIGENSOLUTION_H */
