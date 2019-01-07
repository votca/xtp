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

#ifndef __XTP_GAUSSIAN_QUADRATURE__H
#define __XTP_GAUSSIAN_QUADRATURE__H

#include <votca/xtp/eigen.h>
#include <iostream>
#include <votca/xtp/rpa.h>

//Computes the self-energy expectation matrix for given Kohn-Sham energies and
//M coefficients, using a Gaussian quadrature and expressions for the residuals
//arising from complex contour integration

namespace votca {
    namespace xtp {
       
        class GaussianQuadrature {
        
        public:
            
            //Default constructor
            GaussianQuadrature(const Eigen::VectorXd& qpenergies,const TCMatrix_gwbse& Mmn);
            
            //This function returns the whole self-energy expectation matrix
            //for given Kohn-Sham energies and M coefficients (hence "RPA")
            Eigen::MatrixXcd Sigma(RPA& rpa);
            
            //This function only returns the diagonal of aforementioned matrix
            //Sigma in vector form
            Eigen::VectorXcd SigmaDiag(RPA& rpa);
            
            //This function returns the chosen order (default=12)
            int Order()const{return _order;}
        
        private:
            
            //This function returns the coordinate transformation applied to the
            //quadrature points. These vector entries will serve as frequencies
            //for the dielectric matrix inverses; hence the name
            Eigen::VectorXd CooTfFreq();
            
            //This function calculates the inverse of the microscopic dielectric
            //matrix for given complex frequency and Kohn-Sham energies
            Eigen::MatrixXd CalcDielInv(double frequencyReal, double frequencyImag, RPA& rpa);
            
            //This function returns the B matrix, which contains the sum of
            //the dielectric matrices evaluated at the translated frequency vector 
            //entries from aforementioned vector CooTfFreq
            Eigen::MatrixXd SumDielInvMinId(RPA& rpa);
            
            //This function returns the residual contribution matrix
            Eigen::MatrixXd SigmaRes(RPA& rpa);
            
            //This function returns the diagonal of aforementioned matrix 
            //SigmaRes in vector form
            Eigen::VectorXd SigmaResDiag(RPA& rpa);
            
            //Here, we pick the value of int Order() (default=12)            
            int _order=12;
            
            //Here, we load in the constant vector containing the Kohn-Sham
            //energies, which appears in the constructor. The number of KS wave
            //functions, J, which is the length of that vector, is also defined
            //here
            const Eigen::VectorXd& _qpenergies;
            int noqplevels = _qpenergies.size();
            
            //Here, we load in the M coefficients, which are stored in a matrix
            //array (tensor), which appear in the constructor
            const TCMatrix_gwbse& _Mmn;
            
            //Here, the vectors containing the quadrature evaluation points, 
            //resp. the quadrature weights, are stored
            Eigen::VectorXd _quadpoints;            
            Eigen::VectorXd _quadweights;
                
            };
        }
    }
#endif /* GAUSSIAN_QUADRATURE_H */
