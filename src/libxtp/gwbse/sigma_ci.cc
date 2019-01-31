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

namespace votca {
  namespace xtp {
    
    void Sigma_CI::PrepareScreening(){
         _gq.configure(_qptotal,_qpmin,_homo);
        }

    Eigen::MatrixXd Sigma_CI::CalcCorrelationOffDiag(const Eigen::VectorXd&
            frequencies)const {
        const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
        int DFTsize = energies.size();
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
        //occupied states
        for (int k = 0; k < _homo; ++k) {
            //making sure the M tensor is double-valued
            #if (GWBSE_DOUBLE)
            const Eigen::MatrixXd& MMatrix = _Mmn[k];
            #else
            const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
            #endif
            for (int m = 0; m < _qptotal; ++m) {
                if (frequencies(m) < energies(k)) {
                    //We do not store dielectric inverses now, since we do not 
                    //need all of them
                    Eigen::MatrixXd DielInvMinId = _rpa.calculate_epsilon_r(
                        energies(k)- frequencies(m)).inverse();
                    int nobasisfcs = DielInvMinId.rows();
                    DielInvMinId -= Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                    Eigen::MatrixXd ResPart = MMx * DielInvMinId * MMx.transpose();
                    for (int n = 0; n < _qptotal; ++n) {
                        result(m,n) += ResPart(m,n);
                    }
                }
            }
            for (int n = 0; n < _qptotal; ++n) {
                if (frequencies(n) < energies(k)) {
                    Eigen::MatrixXd DielInvMinId = _rpa.calculate_epsilon_r(
                        energies(k) - frequencies(n)).inverse();
                    int nobasisfcs = DielInvMinId.rows();
                    DielInvMinId -= Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                    Eigen::MatrixXd ResPart = MMx * DielInvMinId * MMx.transpose();
                    for (int m = 0; m < _qptotal; ++m) {
                        result(m,n) += ResPart(m,n);
                    }
                }
            }
        }
        //unoccupied states
        for (int k = _homo; k < DFTsize; k++) {
            #if (GWBSE_DOUBLE)
            const Eigen::MatrixXd& MMatrix = _Mmn[k];
            #else
            const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
            #endif
            for (int m = 0; m < _qptotal; ++m) {
                if (frequencies(m) > energies(k)) {
                    Eigen::MatrixXd DielInvMinId = _rpa.calculate_epsilon_r(
                        energies(k) - frequencies(m)).inverse();
                    int nobasisfcs = DielInvMinId.rows();
                    DielInvMinId -= Eigen::MatrixXd::Identity(nobasisfcs, nobasisfcs);
                    Eigen::MatrixXd ResPart = MMx * DielInvMinId * MMx.transpose();
                        for (int n = 0; n < _qptotal; ++n) {
                            result(m,n) += ResPart(m,n);
                        }
                }
            }
            for (int n = 0; n < _qptotal; ++n) {
                if (frequencies(n) > energies(k)) {
                    Eigen::MatrixXd DielInvMinId = _rpa.calculate_epsilon_r(
                        energies(k) - frequencies(n)).inverse();
                    int nobasisfcs = DielInvMinId.rows();
                    DielInvMinId -= Eigen::MatrixXd::Identity(nobasisfcs, nobasisfcs);
                    Eigen::MatrixXd ResPart = MMx * DielInvMinId * MMx.transpose();
                    for (int m = 0; m < _qptotal; ++m) {
                        result(m,n) += ResPart(m,n);
                    }
                }
            }
        }
        //Now, we add a minus, and the Gauss-Hermite quadrature contribution
        return _gq.SigmaGQ(frequencies, _rpa)-result ;
    }
    
    Eigen::VectorXd Sigma_CI::CalcCorrelationDiag(const Eigen::VectorXd& 
    frequencies)const{
        const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
        int DFTsize=energies.size();
        Eigen::VectorXd result=Eigen::VectorXd::Zero(_qptotal);
        for ( int k = 0 ; k < _homo ; ++k ){
            #if (GWBSE_DOUBLE)
            const Eigen::MatrixXd& MMatrix = _Mmn[k];
            #else
            const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
            #endif
            for (int m = 0 ; m < _qptotal ; ++m){
                if ( frequencies(m) < energies(k) ){
                    Eigen::MatrixXd DielInvMinId = _rpa.calculate_epsilon_r(
                        energies(k) - frequencies(m)).inverse();
                    int nobasisfcs = DielInvMinId.rows();
                    DielInvMinId -= Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                    Eigen::MatrixXd ResPart = MMx * DielInvMinId * MMx.transpose();
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
                    Eigen::MatrixXd DielInvMinId = _rpa.calculate_epsilon_r(
                            energies(k) - frequencies(m)).inverse();
                            int nobasisfcs = DielInvMinId.rows();
                            DielInvMinId -= Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                            Eigen::MatrixXd ResPart = MMx * DielInvMinId * MMx.transpose();
                            result(m) += ResPart(m,m);
                }
            }
        }
        //Doubling factor because m = n
        return _gq.SigmaGQDiag(frequencies,_rpa)-2*result;
        }
  }
}

