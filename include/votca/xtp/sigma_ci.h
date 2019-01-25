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

#ifndef _VOTCA_XTP_SIGMA_CI_H
#define _VOTCA_XTP_SIGMA_CI_H
#include <votca/xtp/eigen.h>
#include <votca/xtp/sigma_base.h>
#include <votca/ctp/logger.h>

#include <votca/xtp/rpa.h>
#include <votca/xtp/gaussian_quadrature.h>

namespace votca {
    namespace xtp {
        
        class TCMatrix_gwbse;
        class RPA;
        
        class Sigma_CI : public Sigma_base {
            
            public:
            
            Sigma_CI(TCMatrix_gwbse& Mmn, RPA& rpa ):Sigma_base(Mmn,rpa),_gq(rpa.getRPAInputEnergies(),Mmn){};
          
            //Sets up the screening parametrisation: empty for the while
            void PrepareScreening();
            
            //This function returns the diagonal of the self-energy correlation
            //expectation matrix for given frequencies, RPA (Kohn-Sham) energies
            //and a calculate Gaussian Quadrature contribution matrix
            Eigen::VectorXd CalcCorrelationDiag(const Eigen::VectorXd&
            frequencies)const;
            //This function returns the whole of the aforementioned matrix
            Eigen::MatrixXd CalcCorrelationOffDiag(const Eigen::VectorXd&
            frequencies)const;
            
            private:
            
            //We add the Gaussian Quadrature contribution to the one from the
            //contour integral, together making the full contribution to Sigma
            GaussianQuadrature _gq;
            
            //Here we use again an imaginary perturbation called eta
            double _eta;
            
            };
    
            }
    
        }

#endif /* _VOTCA_XTP_SIGMA_RESIDUAL_H */
