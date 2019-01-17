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

#ifndef _VOTCA_XTP_SIGMA_SPECTRAL_H
#define _VOTCA_XTP_SIGMA_SPECTRAL_H

#include <complex>
#include <cmath>
#include <votca/xtp/sigma_base.h>
#include <votca/xtp/vc2index.h>
#include <votca/xtp/rpa.h>

namespace votca {
namespace xtp {

class TCMatrix_gwbse;
class RPA;

class Sigma_Spectral : public Sigma_base {
    
public:
    
    Sigma_Spectral(TCMatrix_gwbse& Mmn, RPA& rpa)
        : Sigma_base(Mmn, rpa), _vc2index(0, 0, 0) {};

    // Sets up the screening parametrisation
    void PrepareScreening();
    // Calculates Sigma_c diag elements
    Eigen::VectorXd CalcCorrelationDiag(const Eigen::VectorXd& frequencies) const;
    // Calculates Sigma_c offdiag elements
    Eigen::MatrixXd CalcCorrelationOffDiag(const Eigen::VectorXd& frequencies) const;

private:
    
    // Internal Options
    bool _HedinApprox; // Hedin's Static Approximation
    
    vc2index _vc2index;
    
    // Eigenvalues, eigenvectors from RPA
    rpa_eigensolution _EigenSol;

    // Used to compute correlation part of sigma (eq. 47, 48)(*)
    Eigen::MatrixXd CalcResidues(int s) const;
    double Equation47(int m, int n, const Eigen::VectorXd& energies, double w, double omega, Eigen::MatrixXd& residues) const;
    double Equation48(int m, int n, double omega, Eigen::MatrixXd& residues) const;
    
    // (*) Bruneval, F. et al. molgw 1: Many-body perturbation theory software
    // for atoms, molecules, and clusters. Computer Physics Communications 208,
    // 149–161 (2016).
};
}
}

#endif /* _VOTCA_XTP_SIGMA_SPECTRAL_H */
