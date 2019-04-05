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
        GaussianQuadrature::options opt;
        opt.homo =_opt.homo;
        opt.order = _opt.order;
        opt.qptotal = _qptotal;
        opt.qpmin = _opt.qpmin;
        opt.rpamin = _opt.rpamin;
         _gq.configure(opt);
        }
    
    Eigen::VectorXd Sigma_CI::CalcCorrelationDiag(const Eigen::VectorXd& 
    frequencies)const{
        Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);
        const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
        int rpatotal=energies.size();
        int auxsize = _Mmn.auxsize();
        for ( int m = 0; m < _qptotal; ++m ){
            Eigen::MatrixXd Rmx = Eigen::MatrixXd::Zero(rpatotal,auxsize);
            #if (GWBSE_DOUBLE)
            const Eigen::MatrixXd& Imx = _Mmn[m];
            #else
            const Eigen::MatrixXd Imx = _Mmn[m].cast<double>();       
            #endif
            for ( int i = 0; i < _opt.homo-_opt.rpamin+1; ++i ){
                double omega = energies(i)-frequencies(m)+_eta;
                //std::cd omega = 0;
                //omega.real() = energies(i)-frequencies(m);
                //omega.imag() = - _eta;
                Eigen::MatrixXd DielMxInv = Eigen::MatrixXd::Zero(auxsize,auxsize);
                DielMxInv = (_rpa.calculate_epsilon_r(omega).inverse()).real();
                Eigen::MatrixXd Jmx = Eigen::MatrixXd::Zero(rpatotal,auxsize);
                if ( frequencies(m) < energies(i) ){
                    for ( int mu = 0; mu < auxsize; mu++ ){
                        Jmx(i,mu) = Imx(i,mu);
                   }
                }
                for ( int mu = 0; mu < auxsize; mu++ ){
                    Rmx(i,mu) = DielMxInv(i,mu)-Jmx(i,mu);
                } 
                }
            for ( int i = _opt.homo-_opt.rpamin+1; i < rpatotal; ++i ){
                double omega = energies(i)-frequencies(m)+_eta;
                //std::cd omega = 0;
                //omega.real() = energies(i)-frequencies(m);
                //omega.imag() = _eta;
                Eigen::MatrixXd DielMxInv = Eigen::MatrixXd::Zero(auxsize,auxsize);
                DielMxInv = (_rpa.calculate_epsilon_r(omega).inverse()).real();
                Eigen::MatrixXd Jmx = Eigen::MatrixXd::Zero(rpatotal,auxsize);
                if ( frequencies(m) > energies(i) ){
                    for ( int mu = 0; mu < auxsize; mu++ ){
                        Jmx(i,mu) = Imx(i,mu);
                   }
                }
                for ( int mu = 0; mu < auxsize; mu++ ){
                    Rmx(i,mu) = DielMxInv(i,mu)-Jmx(i,mu);
                } 
                }
            result(m) = -(Rmx.cwiseProduct(Imx)).sum();
        }
        result += _gq.SigmaGQDiag(frequencies,_rpa);

        return result;
    }

    Eigen::MatrixXd Sigma_CI::CalcCorrelationOffDiag(const Eigen::VectorXd&
            frequencies)const {
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal,_qptotal);
        const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
        int rpatotal=energies.size();
        int auxsize = _Mmn.auxsize();
        for ( int m = 0; m < _qptotal; ++m ){
            Eigen::MatrixXd Rmxm = Eigen::MatrixXd::Zero(rpatotal,auxsize);
            #if (GWBSE_DOUBLE)
            const Eigen::MatrixXd& Imxm = _Mmn[m];
            #else
            const Eigen::MatrixXd Imxm = _Mmn[m].cast<double>();       
            #endif
            for ( int n = 0; n < m; ++n ){
                Eigen::MatrixXd Rmxn = Eigen::MatrixXd::Zero(rpatotal,auxsize);
                #if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& Imxn = _Mmn[n];
                #else
                const Eigen::MatrixXd Imxn = _Mmn[n].cast<double>();       
                #endif
            for ( int i = 0; i < _opt.homo-_opt.rpamin+1; ++i ){
                double omegam = energies(i)-frequencies(m)+_eta;
                //std::cd omegam = 0;
                //omegam.real() = energies(i)-frequencies(m);
                //omegam.imag() = - _eta;
                double omegan = energies(i)-frequencies(n)+_eta;
                //std::cd omegan = 0;
                //omegan.real() = energies(i)-frequencies(n);
                //omegan.imag() = - _eta;
                Eigen::MatrixXd DielMxInvm = Eigen::MatrixXd::Zero(auxsize,auxsize);
                Eigen::MatrixXd DielMxInvn = Eigen::MatrixXd::Zero(auxsize,auxsize);
                DielMxInvm = (_rpa.calculate_epsilon_r(omegam).inverse()).real();
                DielMxInvn = (_rpa.calculate_epsilon_r(omegan).inverse()).real();
                Eigen::MatrixXd Jmxm = Eigen::MatrixXd::Zero(rpatotal,auxsize);
                Eigen::MatrixXd Jmxn = Eigen::MatrixXd::Zero(rpatotal,auxsize);
                if ( frequencies(m) < energies(i) ){
                    for ( int mu = 0; mu < auxsize; mu++ ){
                        Jmxm(i,mu) = Imxm(i,mu);
                   }
                }
                if ( frequencies(n) < energies(i) ){
                    for ( int mu = 0; mu < auxsize; mu++ ){
                        Jmxn(i,mu) = Imxn(i,mu);
                   }
                }
                for ( int mu = 0; mu < auxsize; mu++ ){
                    Rmxm(i,mu) = DielMxInvm(i,mu)-Jmxm(i,mu);
                    Rmxn(i,mu) = DielMxInvn(i,mu)-Jmxn(i,mu);
                } 
                }
            for ( int i = _opt.homo-_opt.rpamin+1; i < rpatotal; ++i ){
                double omegam = energies(i)-frequencies(m)+_eta;
                //std::cd omegam = 0;
                //omegam.real() = energies(i)-frequencies(m);
                //omegam.imag() = _eta;
                double omegan = energies(i)-frequencies(n)+_eta;
                //std::cd omegan = 0;
                //omegan.real() = energies(i)-frequencies(n);
                //omegan.imag() = _eta;
                Eigen::MatrixXd DielMxInvm = Eigen::MatrixXd::Zero(auxsize,auxsize);
                Eigen::MatrixXd DielMxInvn = Eigen::MatrixXd::Zero(auxsize,auxsize);
                DielMxInvm = (_rpa.calculate_epsilon_r(omegam).inverse()).real();
                DielMxInvn = (_rpa.calculate_epsilon_r(omegan).inverse()).real();
                Eigen::MatrixXd Jmxm = Eigen::MatrixXd::Zero(rpatotal,auxsize);
                Eigen::MatrixXd Jmxn = Eigen::MatrixXd::Zero(rpatotal,auxsize);
                if ( frequencies(m) > energies(i) ){
                    for ( int mu = 0; mu < auxsize; mu++ ){
                        Jmxm(i,mu) = Imxm(i,mu);
                   }
                }
                if ( frequencies(n) > energies(i) ){
                    for ( int mu = 0; mu < auxsize; mu++ ){
                        Jmxn(i,mu) = Imxn(i,mu);
                   }
                }
                for ( int mu = 0; mu < auxsize; mu++ ){
                    Rmxm(i,mu) = DielMxInvm(i,mu)-Jmxm(i,mu);
                    Rmxn(i,mu) = DielMxInvn(i,mu)-Jmxn(i,mu);
                } 
                }
            result(m,n) = (Rmxm.cwiseProduct(Imxn)+Rmxn.cwiseProduct(Imxm)).sum()/(-2);
        }
        }
        result += result.transpose();
        result += _gq.SigmaGQ(frequencies,_rpa);
        result.diagonal()=CalcCorrelationDiag(frequencies);
        return result;
    }
    
  }
}

