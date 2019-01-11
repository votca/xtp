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

namespace votca {
    namespace xtp {
   
        class TCMatrix_gwbse;
        class RPA;

        class Sigma_CI : public Sigma_base {
            
            public:
            
            
            Sigma_CI(TCMatrix_gwbse& Mmn):Sigma_base(Mmn){};
  
            //Sets up the screening parametrisation
            void PrepareScreening(const RPA& rpa);
            
            //This function returns the whole self-energy expectation matrix
            //for given Kohn-Sham energies and M coefficients (hence "RPA")
            Eigen::VectorXd CalcCorrelationDiag(const Eigen::VectorXd&
            frequencies, const Eigen::VectorXd& RPAEnergies, RPA& rpa,
            GaussianQuadrature& gq)const;
            //Calculates Sigma_c offdiag elements
            Eigen::MatrixXd CalcCorrelationOffDiag(const Eigen::VectorXd&
            frequencies, const Eigen::VectorXd& RPAEnergies, RPA& rpa,
            GaussianQuadrature& gq)const;
            
            //This function returns the chosen order (default=12)
            int Order()const{return _order;}
            
        
            private:
                
                
                
            
            //This function returns the coordinate transformation applied to the
            //quadrature points. These vector entries will serve as frequencies
            //for the dielectric matrix inverses; hence the name
            //Eigen::VectorXd CooTfFreq();
            
            //This function calculates the inverse of the microscopic dielectric
            //matrix for given complex frequency and Kohn-Sham energies
            //Eigen::MatrixXd CalcDielInv(double frequencyReal,
            //                double frequencyImag, RPA& rpa);
            
            
            //Here, we load in the constant vector containing the Kohn-Sham
            //energies, which appears in the constructor. The number of KS wave
            //functions, J, which is the length of that vector, is also defined
            //here
            const Eigen::VectorXd& _qpenergies;
            int noqplevels = _qpenergies.size();
            
            //Here, we load in the M coefficients, which are stored in a matrix
            //array (tensor), which appear in the constructor
            const TCMatrix_gwbse& _Mmn;
            
            
                                    };
    
            }
    
        }

#endif /* _VOTCA_XTP_SIGMA_RESIDUAL_H */
