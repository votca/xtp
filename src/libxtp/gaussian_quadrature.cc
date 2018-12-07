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

#include <votca/xtp/gaussian_quadrature.h>
#include <votca/xtp/eigen.h>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <votca/xtp/threecenter.h>
#include <unsupported/Eigen/CXX11/Tensor>

namespace votca {
  namespace xtp {
GaussianQuadrature::GaussianQuadrature(const Eigen::VectorXd& qpenergies,
        const TCMatrix_gwbse& Mmn):_qpenergies(qpenergies),_Mmn(Mmn){
    //initialise qudarture weights and points
    _quadpoints=Eigen::VectorXd::Zero(_order);
    _quadweights=Eigen::VectorXd::Zero(_order);
    //fill in quadrature weights and points
    _quadpoints << 0.9815606342467191, 0.9041172563704749, 0.7699026741943047,
  0.5873179542866174, 0.3678314989981801, 0.1252334085114689, -0.1252334085114689,
  0.3678314989981801, -0.5873179542866174, -0.7699026741943047, -0.9041172563704749,
  -0.9815606342467191;
    _quadweights << 0.04717533638651180, 0.1069393259953181, 0.1600783285433461,
  0.2031674267230658, 0.2334925365383548, 0.2491470458134029, 0.2491470458134029,
  0.2334925365383548, 0.2031674267230658, 0.1600783285433461, 0.1069393259953181, 
  0.04717533638651180;
}


Eigen::VectorXd GaussianQuadrature::CooTfFreq() {
    Eigen::ArrayXd Freq_p1  = _quadpoints.array() + 1;
    Eigen::ArrayXd Inv_Freq_m1  = (1-_quadpoints.array()).inverse();
    return (Freq_p1*Inv_Freq_m1).log().matrix();
}

Eigen::MatrixXd GaussianQuadrature::CalcDielInv(double freqReal,double freqImag, RPA& rpa){
    Eigen::VectorXd real=Eigen::VectorXd::Ones(1)*freqReal;
    Eigen::VectorXd imag=Eigen::VectorXd::Ones(1)*freqImag;
    rpa.setScreening(real,imag);
    rpa.calculate_epsilon(_qpenergies,_Mmn);
    const Eigen::MatrixXd& DielMx=rpa.GetEpsilon_i()[0];
    return DielMx.inverse();
}

Eigen::MatrixXd GaussianQuadrature::SumDielInvMinId(RPA& rpa){
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(noqplevels,noqplevels);
    Eigen::VectorXd CooTfFreqs = CooTfFreq();
    for (int k = 0 ; k < _order ; ++k){
        result+=CalcDielInv(0,CooTfFreqs(k),rpa);
    }
    int nobasisfcs=result.size();
    result-=Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
    return result;
}

Eigen::MatrixXd GaussianQuadrature::SigmaRes(RPA& rpa){
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(noqplevels,noqplevels);
    
    for (int m = 0 ; m < noqplevels ; ++m){
        for (int n = 0 ; n < noqplevels ; ++n){
            for (int k = 0 ; k < m - 1 ; ++k){
                if (_qpenergies(k) > 0){
                    #if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& MMatrix = _Mmn[k];
#else
                const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
#endif
                Eigen::MatrixXd DielInvMinId = CalcDielInv(_qpenergies(k)-_qpenergies(m),0,rpa);
                int nobasisfcs=DielInvMinId.size();        
                DielInvMinId-=Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                Eigen::MatrixXd ResPart=MMx*DielInvMinId*MMx.transpose();
                result(m,n)+=ResPart(m,n);
                }
            }
            for (int k = m ; k < noqplevels ; ++k){
                if (_qpenergies(k) < 0){
                    #if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& MMatrix = _Mmn[k];
#else
                const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
#endif
                Eigen::MatrixXd DielInvMinId = CalcDielInv(_qpenergies(k)-_qpenergies(m),0,rpa);
                int nobasisfcs=DielInvMinId.size();        
                DielInvMinId-=Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                Eigen::MatrixXd ResPart=MMx*DielInvMinId*MMx.transpose();
                result(m,n)+=ResPart(m,n);
                }
            }
                        for (int l = 0 ; l < n - 1 ; ++l){
                if (_qpenergies(l) > 0){
                    #if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& MMatrix = _Mmn[k];
#else
                const Eigen::MatrixXd MMx = _Mmn[l].cast<double>();       
#endif
                Eigen::MatrixXd DielInvMinId = CalcDielInv(_qpenergies(l)-_qpenergies(m),0,rpa);
                int nobasisfcs=DielInvMinId.size();        
                DielInvMinId-=Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                Eigen::MatrixXd ResPart=MMx*DielInvMinId*MMx.transpose();
                result(m,n)+=ResPart(m,n);
                }
            }
            for (int l = n ; l < noqplevels ; ++l){
                if (_qpenergies(l) < 0){
                    #if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& MMatrix = _Mmn[k];
#else
                const Eigen::MatrixXd MMx = _Mmn[l].cast<double>();       
#endif
                Eigen::MatrixXd DielInvMinId = CalcDielInv(_qpenergies(l)-_qpenergies(m),0,rpa);
                int nobasisfcs=DielInvMinId.size();        
                DielInvMinId-=Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                Eigen::MatrixXd ResPart=MMx*DielInvMinId*MMx.transpose();
                result(m,n)+=ResPart(m,n);
                }
            }
    }
    }
            return (-1)*result;
        }

Eigen::VectorXd GaussianQuadrature::SigmaResDiag(RPA& rpa){
    Eigen::VectorXd result=Eigen::VectorXd::Zero(noqplevels,noqplevels);
    
    for (int m = 0 ; m < noqplevels ; ++m){
            for (int k = 0 ; k < m - 1 ; ++k){
                if (_qpenergies(k) > 0){
                    #if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& MMatrix = _Mmn[k];
#else
                const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
#endif
                Eigen::MatrixXd DielInvMinId = CalcDielInv(_qpenergies(k)-_qpenergies(m),0,rpa);
                int nobasisfcs=DielInvMinId.size();        
                DielInvMinId-=Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                Eigen::MatrixXd ResPart=MMx*DielInvMinId*MMx.transpose();
                result(m)+=ResPart(m,m);
                }
            }
            for (int k = m ; k < noqplevels ; ++k){
                if (_qpenergies(k) < 0){
                    #if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& MMatrix = _Mmn[k];
#else
                const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
#endif
                Eigen::MatrixXd DielInvMinId = CalcDielInv(_qpenergies(k)-_qpenergies(m),0,rpa);
                int nobasisfcs=DielInvMinId.size();        
                DielInvMinId-=Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                Eigen::MatrixXd ResPart=MMx*DielInvMinId*MMx.transpose();
                result(m)+=ResPart(m,m);
                }
            }
        
    }
            return (-2)*result;
        }

        Eigen::VectorXcd GaussianQuadrature::SigmaDiag(RPA& rpa) {
            
            Eigen::VectorXcd result= Eigen::MatrixXcd::Zero(noqplevels,noqplevels);
            Eigen::VectorXcd SigmaGQDiag= Eigen::VectorXcd::Zero(noqplevels,noqplevels);
            Eigen::VectorXd CooTfFreqs = CooTfFreq();
            Eigen::MatrixXd SummedDielInvMinId = SumDielInvMinId(rpa);

            for (int k = 0; k < noqplevels; ++k) {
                Eigen::MatrixXd ResFreqsReal = Eigen::MatrixXd::Zero(_order, noqplevels);
                Eigen::MatrixXd ResFreqsImag = Eigen::MatrixXd::Zero(_order, noqplevels);
#if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& MMatrix = _Mmn[ k ];
#else
                const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
#endif
                Eigen::MatrixXd MMxXSumInvMinId = MMx*SummedDielInvMinId;
                for (int m = 0; m < noqplevels; ++m) {
                    for (int j = 0; j < _order; ++j) {
                        Eigen::VectorXd DeltaE = _qpenergies.array() - _qpenergies(k);
                        Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                        Eigen::VectorXd CooTfFreqsSq = CooTfFreqs.cwiseAbs2();
                        // no cwise stuff, since we deal with different dimensions for j and m
                        ResFreqsReal(j, m) = CooTfFreqs(j) / (DeltaESq(m) + CooTfFreqsSq(j));
                        ResFreqsImag(j, m) = DeltaE(m) / (DeltaESq(m) + CooTfFreqsSq(j));
                    }
                }
                SigmaGQDiag.real() += (_quadweights.transpose() * ResFreqsReal).asDiagonal()*(MMxXSumInvMinId * MMx.transpose());
                SigmaGQDiag.imag() += (_quadweights.transpose() * ResFreqsImag).asDiagonal()*(MMxXSumInvMinId * MMx.transpose());
                }
            return 2*SigmaGQDiag+SigmaResDiag(rpa);
        }


        Eigen::MatrixXcd GaussianQuadrature::Sigma(RPA& rpa) {
            
            Eigen::MatrixXcd result= Eigen::MatrixXcd::Zero(noqplevels,noqplevels);
            Eigen::MatrixXcd SigmaGQ= Eigen::MatrixXcd::Zero(noqplevels,noqplevels);
            Eigen::VectorXd CooTfFreqs = CooTfFreq();
            Eigen::MatrixXd SummedDielInvMinId = SumDielInvMinId(rpa);

            for (int k = 0; k < noqplevels; ++k) {
                Eigen::MatrixXd ResFreqsReal = Eigen::MatrixXd::Zero(_order, noqplevels);
                Eigen::MatrixXd ResFreqsImag = Eigen::MatrixXd::Zero(_order, noqplevels);
#if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& MMatrix = _Mmn[ k ];
#else
                const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
#endif
                Eigen::MatrixXd MMxXSumInvMinId = MMx*SummedDielInvMinId;
                for (int m = 0; m < noqplevels; ++m) {
                    for (int j = 0; j < _order; ++j) {
                        Eigen::VectorXd DeltaE = _qpenergies.array() - _qpenergies(k);
                        Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                        Eigen::VectorXd CooTfFreqsSq = CooTfFreqs.cwiseAbs2();
                        // no cwise stuff, since we deal with different dimensions for j and m
                        ResFreqsReal(j, m) = CooTfFreqs(j) / (DeltaESq(m) + CooTfFreqsSq(j));
                        ResFreqsImag(j, m) = DeltaE(m) / (DeltaESq(m) + CooTfFreqsSq(j));
                    }
                }
                SigmaGQ.real() += (_quadweights.transpose() * ResFreqsReal).asDiagonal()*(MMxXSumInvMinId * MMx.transpose());
                SigmaGQ.real() += (MMxXSumInvMinId * MMx.transpose())*(_quadweights.transpose() * ResFreqsReal).asDiagonal();
                SigmaGQ.imag() += (_quadweights.transpose() * ResFreqsImag).asDiagonal()*(MMxXSumInvMinId * MMx.transpose());
                SigmaGQ.imag() += (MMxXSumInvMinId * MMx.transpose())*(_quadweights.transpose() * ResFreqsImag).asDiagonal();
            }
            return SigmaGQ+SigmaRes(rpa);
        }

}
}
