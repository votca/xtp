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

namespace votca {
  namespace xtp {
      
    //Constructor  
    GaussianQuadrature::GaussianQuadrature(const Eigen::VectorXd& energies,
        const TCMatrix_gwbse& Mmn):_energies(energies),_Mmn(Mmn){
        //we first initialise the Gaussian quadrature weights and evaluation
        //points
        _quadpoints = Eigen::VectorXd::Zero(_order);
        _quadweights = Eigen::VectorXd::Zero(_order);
        //Now, we fill in the quadrature weights and points, which happens here
        //specifically for order n=12
        _quadpoints << 0.9815606342467191, 0.9041172563704749, 
            0.7699026741943047, 0.5873179542866174, 0.3678314989981801,
            0.1252334085114689, -0.1252334085114689, 0.3678314989981801,
            -0.5873179542866174, -0.7699026741943047, -0.9041172563704749,
            -0.9815606342467191;
        _quadweights << 0.04717533638651180, 0.1069393259953181,
            0.1600783285433461, 0.2031674267230658, 0.2334925365383548,
            0.2491470458134029, 0.2491470458134029, 0.2334925365383548,
            0.2031674267230658, 0.1600783285433461, 0.1069393259953181,
            0.04717533638651180;
        }
    
    //This function returns the adapted Gaussian quadrature weights
    Eigen::VectorXd GaussianQuadrature::AdaptedWeights() const{
        //We make temporarily arrays to enable scalar addition and component
        //-wise operations like multiplication and division, and taking the
        //squares
        Eigen::ArrayXd one_min_quadpoints_sqd = 1-(_quadpoints.array()).square();
        //Now, we eventually return to vector form
        return (_quadweights.array() * (one_min_quadpoints_sqd.inverse())).matrix();
    }

    //This function returns the coordinate transformation applied to the
    //quadrature points. These vector entries will serve as *frequencies*
    //for the dielectric matrix inverses; hence the name    
    Eigen::VectorXd GaussianQuadrature::CooTfFreq() const{
        //We make temporary arrays to enable scalar addition and component
        //-wise operations like multiplication and division, and taking the
        //standard logarithm
        Eigen::ArrayXd quadpoints_p1 = _quadpoints.array() + 1;
        Eigen::ArrayXd Inv_quadpoints_m1 = (1 - _quadpoints.array()).inverse();
        //Now, we eventually return to vector form
        return (quadpoints_p1 * Inv_quadpoints_m1).log().matrix();
        }

    //This function calculates the inverses of the microscopic dielectric
    //matrix we need, stored in a vector
    std::vector<Eigen::MatrixXd> GaussianQuadrature::CalcDielInvVector(const RPA& rpa)const{
        std::vector<Eigen::MatrixXd> result;
        //We load in the vector containing the frequencies for which we want
        //the dielectric matrix inverses, and we use a command from rpa.h
        Eigen::VectorXd CooTfFreqs = CooTfFreq();
        for ( int i = 0 ; i < _order ; i++ ) {
            result.push_back(rpa.calculate_epsilon_i(CooTfFreqs(i)).inverse());
            }
        return result;
        }

    //This function returns the "B matrix", which contains the sum of
    //the dielectric matrices evaluated at the translated frequency vector 
    //entries from aforementioned vector CooTfFreq, minus the identity
    Eigen::MatrixXd GaussianQuadrature::SumDielInvMinId(const RPA& rpa)const{
        //First, we intialise the result, which is a sum, at zero
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal,_qptotal);
        //We load in the vector containing the dielectric matrix inverses
        std::vector<Eigen::MatrixXd> DielInvVector = CalcDielInvVector(rpa);
        //Now, using a for loop, we get the sum of all those matrices
            for (int k = 0 ; k < _order ; ++k){
                result += DielInvVector[k];
                }
        //Lastly, we subtract the identity, yielding the B matrix
        result -= Eigen::MatrixXd::Identity(_qptotal,_qptotal);
        return result;
    }
    
    //This function returns the whole quadrature contribution matrix as part of
    //the correlation part of the self-energy
    Eigen::MatrixXd GaussianQuadrature::SigmaGQ(const Eigen::VectorXd&
        frequencies, const RPA& rpa)const{
        //We first initialise the quadrature contribution matrix, which appears
        //in the result and is a sum, at zero
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal,_qptotal);
    //Now, we will set up the adapted frequency tensor in the following
        //for loop, as an array of matrices. We split this tensor in a real and 
        //an imaginary part, and we will use the help matrix B from the function
        //SumDielInvMinId to this end
        Eigen::MatrixXd SummedDielInvMinId = SumDielInvMinId(rpa);
        Eigen::VectorXd AdapWeights = AdaptedWeights();
        Eigen::VectorXd CooTfFreqs = CooTfFreq();
        const double eta = rpa.getEta();
        //To obtain the part from the occupied states, where the energies are
        //negative, we take the indices that correspond to the DFT/RPA energies
        //being negative (where the indices are between 1 and HOMO): this 
        //results in a shift of QPmin to the left for the frequencies
            for (int k = 0; k < _homo - _qpmin; ++k) {
                //Now, we fill the k'th matrix in the array in the following for 
                //loops, where we will use the vector of frequencies from the
                //function CooTfFreq, and the adapted GQ weights
                Eigen::MatrixXd ResFreqs =
                    Eigen::MatrixXd::Zero(_order, _qptotal);
                //We move to an array to have a constant component-
                //wise subtraction, and then return to vector form
                Eigen::VectorXd DeltaE = frequencies.array() - _energies(k);
                //component-wise squaring for the denominator
                Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                //idem. Here, the imaginary perturbation is positive
                
                Eigen::VectorXd CooTfFreqsPlusEta = CooTfFreqs.array()+eta;
                Eigen::VectorXd CooTfFreqsPlusEtaSq = CooTfFreqsPlusEta.cwiseAbs2();
                // now, we cannot do a component-wise operation, since
                //we deal with different dimensions for j and m in 
                //the matrix. So, we fill the matrix entry via a for loop
                for (int m = 0; m < _qptotal; ++m) {
                        for (int j = 0; j < _order; ++j) {
                            ResFreqs(j, m) =
                                CooTfFreqsPlusEta(j) / (DeltaESq(m) + CooTfFreqsPlusEtaSq(j));
                            }
                        }
                //We make sure the M tensor is double-valued
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMatrix = _Mmn[ k ];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
                    #endif
                //We now pre-compute one matrix product, in order to have only 
                //one additional product before we take the diagonal
                Eigen::MatrixXd MMxXSumInvMinId = MMx * SummedDielInvMinId;
                //Now, we can add the resulting matrices to the quadrature
                //contribution matrix, which happens in parts
                result +=
                    (AdapWeights.transpose() * ResFreqs).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose());
                result +=
                    (MMxXSumInvMinId * MMx.transpose()) * 
                    (AdapWeights.transpose() * ResFreqs).asDiagonal();
                }
        
        //Now, the same holds for the indices for the unoccupied states
                    for (int k = _homo - _qpmin; k < _qptotal; ++k) {
                //Now, we fill the k'th matrix in the array in the following for 
                //loops, where we will use the vector of frequencies from the
                //function CooTfFreq.
                Eigen::MatrixXd ResFreqs =
                    Eigen::MatrixXd::Zero(_order, _qptotal);
                //We move to an array to have a constant component-
                //wise subtraction, and then return to vector form
                Eigen::VectorXd DeltaE = frequencies.array() - _energies(k);
                //component-wise squaring for the denominator
                Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                //idem. Here, the imaginary perturbation is negative
                Eigen::VectorXd CooTfFreqsMinusEta = CooTfFreqs.array()-eta;
                Eigen::VectorXd CooTfFreqsMinusEtaSq = CooTfFreqsMinusEta.cwiseAbs2();
                // now, we cannot do a component-wise operation, since
                //we deal with different dimensions for j and m in 
                //the matrix. So, we move to a for loop for filling in the matrix entry
                    for (int m = 0; m < _qptotal; ++m) {
                        for (int j = 0; j < _order; ++j) {
                            ResFreqs(j, m) =
                                CooTfFreqsMinusEta(j) / (DeltaESq(m) + CooTfFreqsMinusEtaSq(j));
                            }
                        }
                //We make sure the M tensor is double-valued
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMatrix = _Mmn[ k ];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
                    #endif
                //We now pre-compute one matrix product, in order to have only 
                //one additional product before we take the diagonal
                Eigen::MatrixXd MMxXSumInvMinId = MMx * SummedDielInvMinId;
                //Now, we can add the resulting matrices to the quadrature
                //contribution matrix, which happens in parts
                result +=
                    (AdapWeights.transpose() * ResFreqs).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose());
                result +=
                    (MMxXSumInvMinId * MMx.transpose()) * 
                    (AdapWeights.transpose() * ResFreqs).asDiagonal();
                }
        return result;
        }
        
    //This function only returns the diagonal of aforementioned matrix
    //SigmaGQ in vector form. All happens analogously
    Eigen::VectorXd GaussianQuadrature::SigmaGQDiag(const Eigen::VectorXd&
        frequencies, const RPA& rpa)const{    
        Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);
        Eigen::VectorXd CooTfFreqs = CooTfFreq();
        Eigen::VectorXd AdapWeights = AdaptedWeights();
        Eigen::MatrixXd SummedDielInvMinId = SumDielInvMinId(rpa);
        const double eta = rpa.getEta();
        for (int k = 0; k < _homo - _qpmin; ++k) {
                Eigen::MatrixXd ResFreqs =
                    Eigen::MatrixXd::Zero(_order, _qptotal);
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMatrix = _Mmn[ k ];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
                    #endif
                Eigen::MatrixXd MMxXSumInvMinId = MMx * SummedDielInvMinId;
                    for (int m = 0; m < _qptotal; ++m) {
                        for (int j = 0; j < _order; ++j) {
                            Eigen::VectorXd DeltaE =
                                frequencies.array() - _energies(k);
                            Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                            Eigen::VectorXd CooTfFreqsPlusEta = 
                                CooTfFreqs.array() + eta;
                            Eigen::VectorXd CooTfFreqsPlusEtaSq = 
                                CooTfFreqsPlusEta.cwiseAbs2();
                            ResFreqs(j, m) = 
                                CooTfFreqsPlusEta(j) / (DeltaESq(m) + CooTfFreqsPlusEtaSq(j));
                            }
                        }
                //There is a similar, double contribution, since m = n on the
                //diagonal: this will be evened out by the factor 2 at the end
                result += 
                    ((AdapWeights.transpose() * ResFreqs).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose())).diagonal();
                }
                    for (int k = _homo - _qpmin; k < _qptotal; ++k) {
                Eigen::MatrixXd ResFreqs =
                    Eigen::MatrixXd::Zero(_order, _qptotal);
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMatrix = _Mmn[ k ];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
                    #endif
                Eigen::MatrixXd MMxXSumInvMinId = MMx * SummedDielInvMinId;
                    for (int m = 0; m < _qptotal; ++m) {
                        for (int j = 0; j < _order; ++j) {
                            Eigen::VectorXd DeltaE =
                                frequencies.array() - _energies(k);
                            Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                            Eigen::VectorXd CooTfFreqsMinusEta = 
                                CooTfFreqs.array() - eta;
                            Eigen::VectorXd CooTfFreqsMinusEtaSq = 
                                CooTfFreqsMinusEta.cwiseAbs2();
                            ResFreqs(j, m) = 
                                CooTfFreqsMinusEta(j) / (DeltaESq(m) + CooTfFreqsMinusEtaSq(j));
                            }
                        }
                //There is a similar, double contribution, since m = n on the
                //diagonal: this will be evened out by the factor 2 at the end
                result += 
                    ((AdapWeights.transpose() * ResFreqs).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose())).diagonal();
                }
        //we get a factor two now, since m = n on the diagonal
        return 2*result;
        }
    }
}
