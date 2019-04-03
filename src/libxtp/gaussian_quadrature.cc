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
#include <votca/xtp/rpa.h>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <votca/xtp/threecenter.h>
#include <votca/xtp/gauss_hermite_quadrature_constants.h>

namespace votca {
  namespace xtp {
      
    //Constructor  
    GaussianQuadrature::GaussianQuadrature(const Eigen::VectorXd& energies,
        const TCMatrix_gwbse& Mmn):_energies(energies),_Mmn(Mmn){}
    
    void GaussianQuadrature::configure(options opt){
                _opt=opt;
                Gauss_Hermite_Quadrature_Constants ghqc;
                _quadpoints = ghqc.getPoints(_opt.order);
                _quadweights = ghqc.getWeights(_opt.order);
            }
    
    Eigen::VectorXd GaussianQuadrature::AdaptedWeights() const{
        //We temporarily move to arrays to enable component-wise operations 
        Eigen::ArrayXd quadpoints_sq = (_quadpoints.array()).square();
        Eigen::ArrayXd quadpoints_sq_exp = (quadpoints_sq).exp();
        return (_quadweights.array() * (quadpoints_sq_exp)).matrix();
    }

    //This function calculates and stores inverses of the microscopic dielectric
    //matrix in a matrix vector
    std::vector<Eigen::MatrixXd> GaussianQuadrature::CalcDielInvVector(const RPA& rpa)const{
        std::vector<Eigen::MatrixXd> result;
        for ( int i = 0 ; i < _opt.order ; i++ ) {
            result.push_back(rpa.calculate_epsilon_i(_quadpoints(i)).inverse());
            }
        return result;
        }

    //This function returns the sum of the matrix vector minus the identity
    Eigen::MatrixXd GaussianQuadrature::SumDielInvMinId(const RPA& rpa)const{
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_Mmn.auxsize(),_Mmn.auxsize());
        std::vector<Eigen::MatrixXd> DielInvVector = CalcDielInvVector(rpa);
            for (int k = 0 ; k < _opt.order ; ++k){
                result += DielInvVector[k];
                }
        result -= Eigen::MatrixXd::Identity(_Mmn.auxsize(),_Mmn.auxsize());
        return result;
    }
    
    Eigen::MatrixXd GaussianQuadrature::SigmaGQ(const Eigen::VectorXd&
        frequencies, const RPA& rpa)const{
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_opt.qptotal,_opt.qptotal);
        Eigen::MatrixXd SummedDielInvMinId = SumDielInvMinId(rpa);
        Eigen::VectorXd AdapWeights = AdaptedWeights();
        const double eta = rpa.getEta();
        //For occupied states, where the RPA energies are negative, we take the 
        //indices that correspond to the DFT/RPA energies being negative (where
        //the indices are between 1 and HOMO): this results in a shift of QPmin
        //to the left for the frequencies
            for (int k = 0; k < _opt.homo - _opt.qpmin + _opt.rpamin; ++k) {
                //Now, we fill the k'th matrix in the adapted frequency tensor
                Eigen::MatrixXd ResFreqs = Eigen::MatrixXd::Zero(_opt.order, _opt.qptotal);
                //We temporarily move to arrays for component-wise operations
                Eigen::VectorXd DeltaE = frequencies.array() - _energies(k);
                Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                //the imaginary perturbation is positive for occupied states
                Eigen::VectorXd QuadPointsPlusEta = _quadpoints.array()+eta;
                Eigen::VectorXd QuadPointsPlusEtaSq = QuadPointsPlusEta.cwiseAbs2();
                //We cannot do a component-wise operation now, because we deal
                //with different dimensions for j and m in the matrix
                for (int m = 0; m < _opt.qptotal; ++m) {
                        for (int j = 0; j < _opt.order; ++j) {
                            ResFreqs(j, m) =
                                QuadPointsPlusEta(j) / (DeltaESq(m) + QuadPointsPlusEtaSq(j));
                            }
                        }
                //We make sure the M tensor is double-valued
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMx = _Mmn[ k ];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
                    #endif
                //We now pre-compute one matrix product
                Eigen::MatrixXd MMxXSumInvMinId = MMx * SummedDielInvMinId;
                //Now, we can add the resulting matrices to the result
                result +=
                    (AdapWeights.transpose() * ResFreqs).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose());
                result +=
                    (MMxXSumInvMinId * MMx.transpose()) * 
                    (AdapWeights.transpose() * ResFreqs).asDiagonal();
                }
        //Now, we handle the indices for the unoccupied states similarly
                    for (int k = _opt.homo - _opt.qpmin + _opt.rpamin; k < _opt.qptotal; ++k) {
                Eigen::MatrixXd ResFreqs =
                    Eigen::MatrixXd::Zero(_opt.order, _opt.qptotal);
                Eigen::VectorXd DeltaE = frequencies.array() - _energies(k);
                Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                Eigen::VectorXd QuadPointsMinusEta = _quadpoints.array() - eta;
                Eigen::VectorXd QuadPointsMinusEtaSq = QuadPointsMinusEta.cwiseAbs2();
                    for (int m = 0; m < _opt.qptotal; ++m) {
                        for (int j = 0; j < _opt.order; ++j) {
                            ResFreqs(j, m) =
                                QuadPointsMinusEta(j) / (DeltaESq(m) + QuadPointsMinusEtaSq(j));
                            }
                        }
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMx = _Mmn[ k ];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
                    #endif
                Eigen::MatrixXd MMxXSumInvMinId = MMx * SummedDielInvMinId;
                result +=
                    (AdapWeights.transpose() * ResFreqs).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose());
                result +=
                    (MMxXSumInvMinId * MMx.transpose()) * 
                    (AdapWeights.transpose() * ResFreqs).asDiagonal();    
                }
        return 0.5*result;
        }
        
    Eigen::VectorXd GaussianQuadrature::SigmaGQDiag(const Eigen::VectorXd&
        frequencies, const RPA& rpa)const{    
        Eigen::VectorXd result = Eigen::VectorXd::Zero(_opt.qptotal);
        Eigen::VectorXd AdapWeights = AdaptedWeights();
        Eigen::MatrixXd SummedDielInvMinId = SumDielInvMinId(rpa);
        const double eta = rpa.getEta();
        for (int k = 0; k < _opt.homo - _opt.qpmin + _opt.rpamin +1; ++k) {
                Eigen::MatrixXd ResFreqs =
                    Eigen::MatrixXd::Zero(_opt.order, _opt.qptotal);
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMx = _Mmn[ k ];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
                    #endif
                Eigen::MatrixXd MMxXSumInvMinId = MMx * SummedDielInvMinId;
                    for (int m = 0; m < _opt.qptotal; ++m) {
                        for (int j = 0; j < _opt.order; ++j) {
                            Eigen::VectorXd DeltaE =
                                frequencies.array() - _energies(k);
                            Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                            Eigen::VectorXd QuadPointsPlusEta = 
                                _quadpoints.array() + eta;
                            Eigen::VectorXd QuadPointsPlusEtaSq = 
                                QuadPointsPlusEta.cwiseAbs2();
                            ResFreqs(j, m) = 
                                QuadPointsPlusEta(j) / (DeltaESq(m) + QuadPointsPlusEtaSq(j));
                            }
                        }
                //There is a double contribution here, since m = n
                result += 
                    ((AdapWeights.transpose() * ResFreqs).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose())).diagonal();
                }
                    for (int k = _opt.homo - _opt.qpmin + _opt.rpamin +1; k < _opt.qptotal; ++k) {
                Eigen::MatrixXd ResFreqs =
                    Eigen::MatrixXd::Zero(_opt.order, _opt.qptotal);
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMx = _Mmn[ k ];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
                    #endif
                Eigen::MatrixXd MMxXSumInvMinId = MMx * SummedDielInvMinId;
                    for (int m = 0; m < _opt.qptotal; ++m) {
                        for (int j = 0; j < _opt.order; ++j) {
                            Eigen::VectorXd DeltaE =
                                frequencies.array() - _energies(k);
                            Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                            Eigen::VectorXd QuadPointsMinusEta = 
                                _quadpoints.array() - eta;
                            Eigen::VectorXd QuadPointsMinusEtaSq = 
                                QuadPointsMinusEta.cwiseAbs2();
                            ResFreqs(j, m) = 
                                QuadPointsMinusEta(j) / (DeltaESq(m) + QuadPointsMinusEtaSq(j));
                            }
                        }
                result += 
                    ((AdapWeights.transpose() * ResFreqs).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose())).diagonal();
                }
        //Because of the double contributions, there is no factor 1/2 here
        return result;
        }
    }
}
