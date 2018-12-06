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

// Computes the Gaussian quadrature used for numerical integration of the complex integrals used for the self-energy

namespace votca {
    namespace xtp {
       
        class GaussianQuadrature {
        public:
            
            GaussianQuadrature(const Eigen::VectorXd& qpenergies,const TCMatrix_gwbse& Mmn);
            
            Eigen::MatrixXcd Integrate(RPA& rpa);
            
            int Order()const{return _order;}
        private:   
            
            Eigen::MatrixXd CalcInverse(double frequency, RPA& rpa);
            
            Eigen::VectorXd TranslateFrequencies();
            
            Eigen::MatrixXd SumInversesMinusOne(RPA& rpa);
                        
            int _order=12;
            
            const Eigen::VectorXd& _qpenergies;
            int numofqplevels = _qpenergies.size();
            
            const TCMatrix_gwbse& _Mmn;
            
            Eigen::VectorXd _integrationpoints;
            
            Eigen::VectorXd _integrationweights;
            
          
            
            
                        

            
        
        };

    }
}
#endif /* GAUSSIAN_QUADRATURE_H */
