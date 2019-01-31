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

//Computes the contribution from the Gauss-Hermite quadrature to the 
//self-energy expectation matrix for given RPA and frequencies
namespace votca {
    namespace xtp {
    
        class TCMatrix_gwbse;
        class RPA;
       
        class GaussianQuadrature {
        
        public:
            
            GaussianQuadrature(const Eigen::VectorXd& energies,const TCMatrix_gwbse& Mmn);
            
            void configure(int qptotal, int qpmin, int homo){
                _qptotal = qptotal;
                _qpmin = qpmin;
                _homo = homo;
            }
            
            Eigen::MatrixXd SigmaGQ(const Eigen::VectorXd& frequencies,
                const RPA& rpa)const;
            
            Eigen::VectorXd SigmaGQDiag(const Eigen::VectorXd& frequencies,
                const RPA& rpa)const;
            
            int Order()const{return _order;}

        private:
            
            
            Eigen::VectorXd AdaptedWeights() const ;
            
            //This function calculates and stores inverses of the microscopic dielectric
            //matrix in a matrix vector
            std::vector<Eigen::MatrixXd> CalcDielInvVector(const RPA& rpa)const;
            
            //This function returns the sum of the matrix vector minus the identity
            Eigen::MatrixXd SumDielInvMinId(const RPA& rpa)const;

            int _order=12;
            
            const Eigen::VectorXd& _energies;
            int _qptotal;
            int _homo;
            int _qpmin;
            double _eta;
            const TCMatrix_gwbse& _Mmn;
            Eigen::VectorXd _quadpoints;            
            Eigen::VectorXd _quadweights;
            };
        }
    }
#endif /* GAUSSIAN_QUADRATURE_H */
