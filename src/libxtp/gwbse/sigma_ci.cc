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



#include <votca/xtp/sigma_ci.h>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <votca/tools/constants.h>
#include <votca/xtp/threecenter.h>
#include <votca/xtp/eigen.h>
#include <iostream>
#include <vector>
#include <complex>

namespace votca {
  namespace xtp {
    
    //This function initialises the important indices QPtotal, QPmin, and HOMO
    void Sigma_CI::PrepareScreening(){
         _gq.configure(_qptotal,_qpmin,_homo);
        }

    //This function returns the whole self-energy expectation matrix
    //for given Kohn-Sham energies and M coefficients (hence "RPA")
    Eigen::MatrixXd Sigma_CI::CalcCorrelationOffDiag(const Eigen::VectorXd&
            frequencies)const {
        const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
        int DFTsize = energies.size();
        //First, we initialise the result, which is a sum, at zero
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
        //Then, using a double for loop, we fill the resulting matrix. We use 
        //a distinction between occupied and unoccupied states, resulting in 
        //0 < k <= HOMO, resp. HOMO < k <= QPtotal
        for (int k = 0; k < _homo; ++k) {
            //making sure the M tensor is double-valued
            #if (GWBSE_DOUBLE)
            const Eigen::MatrixXd& MMatrix = _Mmn[k];
            #else
            const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
            #endif
            for (int m = 0; m < _qptotal; ++m) {
                if (frequencies(m) < energies(k)) {
                    //We only have to make the dielectric inverse when it is 
                    //needed; hence, we do not use matrix vectors here
                    Eigen::MatrixXd DielInvMinId =
                            _rpa.calculate_epsilon_r(energies(k)
                            - frequencies(m)).inverse();
                    //To subtract the identity, we must know the right dimension
                    int nobasisfcs = DielInvMinId.rows();
                    DielInvMinId -=
                            Eigen::MatrixXd::Identity(nobasisfcs,
                            nobasisfcs);
                    Eigen::MatrixXd ResPart =
                            MMx * DielInvMinId * MMx.transpose();
                    //the resulting matrix entry is the sum of all these residual 
                    //partial contribution, which we fill in in the second for loop
                    for (int n = 0; n < _qptotal; ++n) {
                        result(m, n) += ResPart(m, n);
                    }
                }
            }
            //Analogously, we do the same with m and n interchanged
            for (int n = 0; n < _qptotal; ++n) {
                if (frequencies(n) < energies(k)) {
                    Eigen::MatrixXd DielInvMinId =
                            _rpa.calculate_epsilon_r(energies(k)
                            - frequencies(n)).inverse();
                    int nobasisfcs = DielInvMinId.rows();
                    DielInvMinId -=
                            Eigen::MatrixXd::Identity(nobasisfcs,
                            nobasisfcs);
                    Eigen::MatrixXd ResPart =
                            MMx * DielInvMinId * MMx.transpose();
                    for (int m = 0; m < _qptotal; ++m) {
                        result(m, n) += ResPart(m, n);
                    }
                }
            }
        }
        //Now, we deal with the unoccupied states in a fully similar way
        for (int k = _homo; k < DFTsize; k++) {
            #if (GWBSE_DOUBLE)
            const Eigen::MatrixXd& MMatrix = _Mmn[k];
            #else
            const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
            #endif
            for (int m = 0; m < _qptotal; ++m) {
                if (frequencies(m) > energies(k)) {
                    Eigen::MatrixXd DielInvMinId =
                            _rpa.calculate_epsilon_r(energies(k)
                            - frequencies(m)).inverse();
                    int nobasisfcs = DielInvMinId.rows();
                    DielInvMinId -=
                            Eigen::MatrixXd::Identity(nobasisfcs, nobasisfcs);
                    Eigen::MatrixXd ResPart =
                            MMx * DielInvMinId * MMx.transpose();
                        for (int n = 0; n < _qptotal; ++n) {
                        //the resulting matrix entry is the sum of all these
                        result(m, n) += ResPart(m, n);
                        }
                }
            }
            for (int n = 0; n < _qptotal; ++n) {
                if (frequencies(n) > energies(k)) {
                    Eigen::MatrixXd DielInvMinId =
                            _rpa.calculate_epsilon_r(energies(k)
                            - frequencies(n)).inverse();
                    int nobasisfcs = DielInvMinId.rows();
                    DielInvMinId -=
                            Eigen::MatrixXd::Identity(nobasisfcs, nobasisfcs);
                    Eigen::MatrixXd ResPart =
                            MMx * DielInvMinId * MMx.transpose();
                    for (int m = 0; m < _qptotal; ++m) {
                        //the resulting matrix entry is the sum of all these
                        result(m, n) += ResPart(m, n);
                    }
                }
            }
        }
        //Now, we add a minus, and the Gaussian quadrature contribution
        return _gq.SigmaGQ(frequencies, _rpa)-result ;
    }
    
    //This function returns the diagonal of aforementioned matrix SigmaRes
    //in vector form, which happens fully similar, albeit a factor 2 is added
    Eigen::VectorXd Sigma_CI::CalcCorrelationDiag(const Eigen::VectorXd& 
    frequencies)const{
        const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
        int DFTsize=energies.size();
        //First, we initialise the result, which is a sum, at zero
        Eigen::VectorXd result=Eigen::VectorXd::Zero(_qptotal);
            //Now, everything is the same, but now we only have one sum to split
            //since m = n on the diagonal
                    for ( int k = 0 ; k < _homo ; ++k ){
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMatrix = _Mmn[k];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
                    #endif
                        for (int m = 0 ; m < _qptotal ; ++m){
                            if ( frequencies(m) < energies(k) ){
                            Eigen::MatrixXd DielInvMinId = 
                                _rpa.calculate_epsilon_r(energies(k)
                                    - frequencies(m)).inverse();
                            int nobasisfcs = DielInvMinId.rows();
                            DielInvMinId -= 
                                Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                            Eigen::MatrixXd ResPart = 
                                MMx * DielInvMinId * MMx.transpose();
                            result(m) += ResPart(m,m);
                            }
                        }
                    }
                    for ( int k = _homo ; k < DFTsize ; ++k ){
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMatrix = _Mmn[k];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
                    #endif
                        for (int m = 0 ; m < _qptotal ; ++m){
                            
                            if ( frequencies(m) > energies(k) ){
                            Eigen::MatrixXd DielInvMinId = 
                                _rpa.calculate_epsilon_r(energies(k) -
                                    frequencies(m)).inverse();
                            int nobasisfcs = DielInvMinId.rows();
                            DielInvMinId -= 
                                Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                            Eigen::MatrixXd ResPart = 
                                MMx * DielInvMinId * MMx.transpose();
                            result(m) += ResPart(m,m);
                            }
                        }
                    }
        std::cout << "GQ contribution" << std::endl;
        std::cout << _gq.SigmaGQDiag(frequencies,_rpa) << std::endl;
        std::cout << "residual contribution" << std::endl;
        std::cout << (-2)*result << std::endl;
        
        //Because m = n on the diagonal, we also get a doubling factor in front
        return (-2)*result+_gq.SigmaGQDiag(frequencies,_rpa);
        }
  }
}

