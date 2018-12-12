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
      
    //Default constructor  
    
    GaussianQuadrature::GaussianQuadrature(const Eigen::VectorXd& qpenergies,
        const TCMatrix_gwbse& Mmn):_qpenergies(qpenergies),_Mmn(Mmn){
    
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

    //This function returns the coordinate transformation applied to the
    //quadrature points. These vector entries will serve as frequencies
    //for the dielectric matrix inverses; hence the name    

    Eigen::VectorXd GaussianQuadrature::CooTfFreq() {
    
        //We make temporarily arrays to enable scalar addition and component
        //-wise operations like multiplication and division, and taking the
        //standard logarithm
        Eigen::ArrayXd Freq_p1 = _quadpoints.array() + 1;
        Eigen::ArrayXd Inv_Freq_m1 = (1 - _quadpoints.array()).inverse();
    
        //Now, we eventually return to vector form
        return (Freq_p1 * Inv_Freq_m1).log().matrix();

        }

    //This function calculates the inverse of the microscopic dielectric
    //matrix for given complex frequency and Kohn-Sham energies
    
    Eigen::MatrixXd GaussianQuadrature::CalcDielInv(double freqReal,
            double freqImag, RPA& rpa){
    
        //First, using commands from rpa.h, we set up the dielectric matrix
        //itself
        Eigen::VectorXd real = Eigen::VectorXd::Ones(1) * freqReal;
        Eigen::VectorXd imag = Eigen::VectorXd::Ones(1) * freqImag;
        rpa.setScreening(real,imag);
        rpa.calculate_epsilon(_qpenergies,_Mmn);
        const Eigen::MatrixXd& DielMx = rpa.GetEpsilon_i()[0];
    
        //Now, we just take the inverse of this matrix
        return DielMx.inverse();
    
        }

    //This function returns the B matrix, which contains the sum of
    //the dielectric matrices evaluated at the translated frequency vector 
    //entries from aforementioned vector CooTfFreq
    
    Eigen::MatrixXd GaussianQuadrature::SumDielInvMinId(RPA& rpa){
    
        //First, we intialise the result, which is a sum, at zero
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(noqplevels,noqplevels);
    
        //We load in the vector containing the frequencies for which we want
        //the dielectric matrix inverses
        Eigen::VectorXd CooTfFreqs = CooTfFreq();
    
        //Now, using a for loop, we get the sum of all those matrices
            for (int k = 0 ; k < _order ; ++k){
            
                result += CalcDielInv(0,CooTfFreqs(k),rpa);
            
                }
    
        //Lastly, we subtract the identity, yielding the B matrix
        int nobasisfcs = result.size();
        result -= Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
        return result;
    
    }
    
    //This function returns the residual contribution matrix

    Eigen::MatrixXd GaussianQuadrature::SigmaRes(RPA& rpa){
    
        //First, we initialise the result, which is a sum, at zero
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(noqplevels,noqplevels);
    
        //Then, using a triple for loop, we fill the matrix. We use symmetry,
        //so the second index runs up and including the first index, in order to
        //get the diagonal too.
            for (int m = 0 ; m < noqplevels ; ++m){
                
                for (int n = 0 ; n < m + 1 ; ++n){
                
                    //we split up the two sums over k and l in four ones
                    for (int k = 0 ; k < m - 1 ; ++k){
                    
                        //no else, because there is zero contribution then
                        if (_qpenergies(k) > 0){
                    
                            //making sure the M tensor is double-valued
                            #if (GWBSE_DOUBLE)
                                const Eigen::MatrixXd& MMatrix = _Mmn[k];
                            #else
                                const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
                            #endif
                    
                            //evaluating the dielectric matrix inverse at the right
                            //(real) frequency value
                            Eigen::MatrixXd DielInvMinId =
                                CalcDielInv(_qpenergies(k) - _qpenergies(m),0,rpa);
                    
                            //now, we subtract the identity
                            int nobasisfcs = DielInvMinId.size();        
                            DielInvMinId -=
                                Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                    
                            //now, we compute the matrix multiplication, which is 
                            //the contribution for this part of the loop
                            Eigen::MatrixXd ResPart = 
                                MMx*DielInvMinId*MMx.transpose();
                    
                            //the resulting matrix entry is the sum of all these
                            result(m,n) += ResPart(m,n);
                    
                            //QUESTION: IS THIS REALLY BETTER, JUST MAKING A FULL
                            //MATRIX AND THEN PICKING ONE ELEMENT OF IT? THE 
                            //OTHER OPTION WOULD CONTAIN THREE ADDITIONAL FOR LOOPS.
                            }
                
                        }
                
                    //the other part of the sum over k: all is similar
                    for (int k = m ; k < noqplevels ; ++k){
                    
                        if (_qpenergies(k) < 0){
                        
                            #if (GWBSE_DOUBLE)
                                const Eigen::MatrixXd& MMatrix = _Mmn[k];
                            #else
                                const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
                            #endif
                
                            Eigen::MatrixXd DielInvMinId =
                                CalcDielInv(_qpenergies(k) - _qpenergies(m),0,rpa);
                            int nobasisfcs=DielInvMinId.size();        
                            DielInvMinId -=
                                Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                            Eigen::MatrixXd ResPart =
                                MMx * DielInvMinId * MMx.transpose();
                            result(m,n) += ResPart(m,n);
                        
                            }
                    
                        }
                
                    //the other sum, now over l: all is similar
                    for (int l = 0 ; l < n - 1 ; ++l){
                    
                        if (_qpenergies(l) > 0){
                        
                            #if (GWBSE_DOUBLE)
                                const Eigen::MatrixXd& MMatrix = _Mmn[k];
                            #else
                                const Eigen::MatrixXd MMx = _Mmn[l].cast<double>();       
                            #endif
                
                            Eigen::MatrixXd DielInvMinId =
                                CalcDielInv(_qpenergies(l) - _qpenergies(m),0,rpa);
                            int nobasisfcs = DielInvMinId.size();        
                            DielInvMinId -= 
                                Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                            Eigen::MatrixXd ResPart = 
                                MMx * DielInvMinId * MMx.transpose();
                            result(m,n) += ResPart(m,n);
                
                            }
                    
                        }
                
                    for (int l = n ; l < noqplevels ; ++l){
                
                        if (_qpenergies(l) < 0){
                        
                            #if (GWBSE_DOUBLE)
                                const Eigen::MatrixXd& MMatrix = _Mmn[k];
                            #else
                                const Eigen::MatrixXd MMx = _Mmn[l].cast<double>();       
                            #endif

                            Eigen::MatrixXd DielInvMinId = 
                                CalcDielInv(_qpenergies(l) - _qpenergies(m),0,rpa);
                            int nobasisfcs = DielInvMinId.size();        
                            DielInvMinId -= 
                                Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                            Eigen::MatrixXd ResPart =
                                MMx * DielInvMinId * MMx.transpose();
                            result(m,n) += ResPart(m,n);
                        
                            }
                    
                        }
                
                    }
            
                }
        
        //Now, we just add the transpose matrix, to get the full, symmetric,
        //matrix. Note that we now get the diagonal doubled, so we must still
        //halve it to obtain the real matrix
        result+=result.transpose();
        result.diagonal()*=0.5;
        //There is a factor (-1) to be put in front, and we obtain the result
        return (-1)*result;
    
        }
    
    //This function returns the diagonal of aforementioned matrix SigmaRes
    //in vector form
    
    Eigen::VectorXd GaussianQuadrature::SigmaResDiag(RPA& rpa){
        
        //First, we initialise the result, which is a sum, at zero
        Eigen::VectorXd result=Eigen::VectorXd::Zero(noqplevels,noqplevels);
            
            //Now, everything is the same, but now we only have one sum to split
            //since m = n on the diagonal
            for (int m = 0 ; m < noqplevels ; ++m){
                
                for (int k = 0 ; k < m - 1 ; ++k){
                    
                    if (_qpenergies(k) > 0){
                    
                        #if (GWBSE_DOUBLE)
                            const Eigen::MatrixXd& MMatrix = _Mmn[k];
                        #else
                            const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
                        #endif

                        Eigen::MatrixXd DielInvMinId = 
                            CalcDielInv(_qpenergies(k) - _qpenergies(m),0,rpa);
                        int nobasisfcs = DielInvMinId.size();        
                        DielInvMinId -= 
                            Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                        Eigen::MatrixXd ResPart = 
                            MMx * DielInvMinId * MMx.transpose();
                        result(m) += ResPart(m,m);
                        
                        }
            
                    }
                
                for (int k = m ; k < noqplevels ; ++k){
                
                    if (_qpenergies(k) < 0){
                    
                        #if (GWBSE_DOUBLE)
                            const Eigen::MatrixXd& MMatrix = _Mmn[k];
                        #else
                            const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
                        #endif

                        Eigen::MatrixXd DielInvMinId = 
                            CalcDielInv(_qpenergies(k) - _qpenergies(m),0,rpa);
                        int nobasisfcs = DielInvMinId.size();        
                        DielInvMinId -= 
                            Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                        Eigen::MatrixXd ResPart = 
                            MMx * DielInvMinId * MMx.transpose();
                        result(m) += ResPart(m,m);
                
                        }
                
                    }
        
                }
        
        //Because m = n on the diagonal, we also get a double factor in front
        return (-2)*result;
        
        }

    //This function returns the whole self-energy expectation matrix
    //for given Kohn-Sham energies and M coefficients (hence "RPA")

    Eigen::MatrixXcd GaussianQuadrature::Sigma(RPA& rpa) {
        
        //We first initialise the quadrature contribution matrix, which appears
        //in the result and is a sum, at zero
        Eigen::MatrixXcd SigmaGQ = Eigen::MatrixXcd::Zero(noqplevels,noqplevels);

        //Now, we will set up the adapted frequency tensor in the following
        //for loop, as an array of matrices. We split this tensor in a real and 
        //an imaginary part, and we will use the help matrix B from the function
        //SumDielInvMinId to this end
        Eigen::MatrixXd SummedDielInvMinId = SumDielInvMinId(rpa);
        
            for (int k = 0; k < noqplevels; ++k) {
                
                //Now, we fill the k'th matrix in the array in the following for 
                //loops, where we will use the vector of frequencies from the
                //function CooTfFreq. Also, we initialise the real and imag part
                Eigen::VectorXd CooTfFreqs = CooTfFreq();
                Eigen::MatrixXd ResFreqsReal =
                    Eigen::MatrixXd::Zero(_order, noqplevels);
                Eigen::MatrixXd ResFreqsImag = 
                    Eigen::MatrixXd::Zero(_order, noqplevels);
                
                    for (int m = 0; m < noqplevels; ++m) {
                        
                        for (int j = 0; j < _order; ++j) {
                            
                            //We move to an array to have a constant component-
                            //wise subtraction, and then return to vector form
                            Eigen::VectorXd DeltaE = 
                                _qpenergies.array() - _qpenergies(k);
                            
                            //component-wise squaring for the denominator
                            Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                            
                            //idem
                            Eigen::VectorXd CooTfFreqsSq = 
                                CooTfFreqs.cwiseAbs2();
                            
                            // now, we cannot do a component-wise thing, since
                            //we deal with different dimensions for j and m in 
                            //the matrix
                            ResFreqsReal(j, m) =
                                CooTfFreqs(j) / (DeltaESq(m) + CooTfFreqsSq(j));
                            ResFreqsImag(j, m) = 
                                DeltaE(m) / (DeltaESq(m) + CooTfFreqsSq(j));
                            
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
                SigmaGQ.real() +=
                    (_quadweights.transpose() * ResFreqsReal).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose());
                SigmaGQ.real() +=
                    (MMxXSumInvMinId * MMx.transpose()) * 
                    (_quadweights.transpose() * ResFreqsReal).asDiagonal();
                SigmaGQ.imag() += 
                    (_quadweights.transpose() * ResFreqsImag).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose());
                SigmaGQ.imag() +=
                    (MMxXSumInvMinId * MMx.transpose()) * 
                    (_quadweights.transpose() * ResFreqsImag).asDiagonal();
                
                }
        
        //Now, we add the already computed residual expressions, and we have
        //our result
        return SigmaGQ;//+SigmaRes(rpa);
        
        }
     //QUESTION: HOW TO PROFIT FROM THE SYMMETRY OF THE MATRIX HERE???
     //WE DO NOT MAKE IT ELEMENT AFTER ELEMENT, SO HOW MUST THIS BE EXPLOITED
     //HERE? WE CANNOT TAKE TRIANGULAR PARTS FROM SOME MATRICES AND USE THE SAME
     //METHOD AS BEFORE, BECAUSE THEN THE MATRIXES ARE STILL FULLY MADE  
        
    //This function only returns the diagonal of aforementioned matrix
    //Sigma in vector form. All happens analogously
    
    Eigen::VectorXcd GaussianQuadrature::SigmaDiag(RPA& rpa) {
            
        Eigen::VectorXcd SigmaGQDiag =
            Eigen::VectorXcd::Zero(noqplevels);
        Eigen::VectorXd CooTfFreqs = CooTfFreq();
        Eigen::MatrixXd SummedDielInvMinId = SumDielInvMinId(rpa);

            for (int k = 0; k < noqplevels; ++k) {
                
                Eigen::MatrixXd ResFreqsReal =
                    Eigen::MatrixXd::Zero(_order, noqplevels);
                Eigen::MatrixXd ResFreqsImag = 
                    Eigen::MatrixXd::Zero(_order, noqplevels);
                
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMatrix = _Mmn[ k ];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
                    #endif

                Eigen::MatrixXd MMxXSumInvMinId = MMx * SummedDielInvMinId;
                
                    for (int m = 0; m < noqplevels; ++m) {
                        
                        for (int j = 0; j < _order; ++j) {
                            
                            Eigen::VectorXd DeltaE =
                                _qpenergies.array() - _qpenergies(k);
                            Eigen::VectorXd DeltaESq = DeltaE.cwiseAbs2();
                            Eigen::VectorXd CooTfFreqsSq =
                                CooTfFreqs.cwiseAbs2();
                            ResFreqsReal(j, m) = 
                                CooTfFreqs(j) / (DeltaESq(m) + CooTfFreqsSq(j));
                            ResFreqsImag(j, m) = 
                                DeltaE(m) / (DeltaESq(m) + CooTfFreqsSq(j));
                            
                            }
                        
                        }
                
                //There is a similar, double contribution, since m = n on the
                //diagonal: this will be evened out by the factor 2 at the end
                SigmaGQDiag.real() += 
                    ((_quadweights.transpose() * ResFreqsReal).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose())).diagonal();
                SigmaGQDiag.imag() += 
                    ((_quadweights.transpose() * ResFreqsImag).asDiagonal() * 
                    (MMxXSumInvMinId * MMx.transpose())).diagonal();
                
                }
        
        //we get a factor two now, since m = n on the diagonal
        return 2*SigmaGQDiag;//+SigmaResDiag(rpa);
        
        }


    //This function checks the answer using the alternative method, 
    //using a sextuple for loop
    Eigen::MatrixXd GaussianQuadrature::CheckFunction(RPA& rpa){
        Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(noqplevels,noqplevels);
        Eigen::MatrixXcd resultGQ=Eigen::MatrixXcd::Zero(noqplevels,noqplevels);
        Eigen::MatrixXd resultGQReal=Eigen::MatrixXd::Zero(noqplevels,noqplevels);
        Eigen::MatrixXd resultGQImag =Eigen::MatrixXd::Zero(noqplevels,noqplevels);
        Eigen::MatrixXd resultRes=Eigen::MatrixXd::Zero(noqplevels,noqplevels);
                
        Eigen::VectorXd CooTfFreqs = CooTfFreq();
        
            for  ( int k = 0 ; k < noqplevels ; ++k ) {
                
                #if (GWBSE_DOUBLE)
                    const Eigen::MatrixXd& MMatrix = _Mmn[ k ];
                #else
                    const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();
                #endif

                for  ( int m = 0 ; m < noqplevels ; ++m ) {
            
                    for  ( int n = 0 ; n < noqplevels ; ++n ) {
                    
                        for  ( int j = 0 ; j < _order ; ++j ) {
                            
                            Eigen::MatrixXd DielInvGQ=CalcDielInv(0,CooTfFreqs(j),rpa);
                            int nobasisfcs = DielInvGQ.size();
                            Eigen::MatrixXd DielInvResM=CalcDielInv(_qpenergies(k)-_qpenergies(m),0,rpa);
                            Eigen::MatrixXd DielInvResN=CalcDielInv(_qpenergies(k)-_qpenergies(n),0,rpa);
                            
                                
                                for  ( int mu = 0 ; mu < nobasisfcs ; ++mu ) {
            
                                    for  ( int nu = 0 ; nu < nobasisfcs ; ++nu ) {
                                    
                                        if ( mu == nu ) {
                                            
                                            if ( m == n ) {
                                            
resultGQReal(m,n) += 2 * _quadweights(j) * MMx(m,mu) * MMx(n,nu) * (DielInvGQ(mu,nu)-1) / 
        (std::pow(_qpenergies(m) - _qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        CooTfFreqs(j);
resultGQImag(m,n) += 2 * _quadweights(j) * MMx(m,mu) * MMx(n,nu) * (DielInvGQ(mu,nu)-1) / 
        (std::pow(_qpenergies(m)-_qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        (_qpenergies(m) - _qpenergies(k) );

                                                if ( _qpenergies(k) < 0 ) {
                                                    if ( _qpenergies(m) < _qpenergies(k) ) {
                                                
resultRes(m,n) -= 2 * MMx(m,mu) * MMx(n,nu) * (DielInvResM(mu,nu) - 1);

                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }
                                                else {
                                                    if ( _qpenergies(m) > _qpenergies(k) ) {
resultRes(m,n) -= 2 * MMx(m,mu) * MMx(n,nu) * (DielInvResM(mu,nu) - 1);   
                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }
                                            
resultGQReal(m,n) += _quadweights(j) * (MMx(m,mu) * MMx(n,nu) * (DielInvGQ(mu,nu)-1) / 
        (std::pow(_qpenergies(m) - _qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        CooTfFreqs(j) + MMx(m,mu) * MMx(n,nu) * (DielInvGQ(mu,nu)-1) / 
        (std::pow(_qpenergies(n) - _qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        CooTfFreqs(j));
resultGQImag(m,n) += _quadweights(j) * (MMx(m,mu) * MMx(n,nu) * (DielInvGQ(mu,nu)-1) / 
        (std::pow(_qpenergies(m)-_qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        (_qpenergies(m) - _qpenergies(k) ) + MMx(m,mu) * MMx(n,nu) * (DielInvGQ(mu,nu)-1) / 
        (std::pow(_qpenergies(n)-_qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        (_qpenergies(n) - _qpenergies(k)));
                                                if ( _qpenergies(k) < 0 ) {
                                                    if ( _qpenergies(m) < _qpenergies(k) ) {
                                                
resultRes(m,n) -= MMx(m,mu) * MMx(n,nu) * (DielInvResM(mu,nu) - 1);

                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }
                                                else {
                                                    if ( _qpenergies(m) > _qpenergies(k) ) {
resultRes(m,n) -= MMx(m,mu) * MMx(n,nu) * (DielInvResM(mu,nu) - 1);   
                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }

                                                if ( _qpenergies(k) < 0 ) {
                                                    if ( _qpenergies(n) < _qpenergies(k) ) {
                                                
resultRes(m,n) -= MMx(m,mu) * MMx(n,nu) * (DielInvResN(mu,nu) - 1);

                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }
                                                else {
                                                    if ( _qpenergies(n) > _qpenergies(k) ) {
resultRes(m,n) -= MMx(m,mu) * MMx(n,nu) * (DielInvResN(mu,nu) - 1);   
                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }
                                            
                                                }
                                                                                        
                                            }
                                        
                                        else {
                                        
                                        if ( m == n ) {
                                        
resultGQReal(m,n) += 2 * _quadweights(j) * MMx(m,mu) * MMx(n,nu) * DielInvGQ(mu,nu) / 
        (std::pow(_qpenergies(m) - _qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        CooTfFreqs(j);
resultGQImag(m,n) += 2 * _quadweights(j) * MMx(m,mu) * MMx(n,nu) * DielInvGQ(mu,nu) / 
        (std::pow(_qpenergies(m)-_qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        (_qpenergies(m) - _qpenergies(k) );

 if ( _qpenergies(k) < 0 ) {
                                                    if ( _qpenergies(m) < _qpenergies(k) ) {
                                                
resultRes(m,n) -= 2 * MMx(m,mu) * MMx(n,nu) * (DielInvResM(mu,nu) - 1);

                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }
                                                else {
                                                    if ( _qpenergies(m) > _qpenergies(k) ) {
resultRes(m,n) -= 2 * MMx(m,mu) * MMx(n,nu) * (DielInvResM(mu,nu) - 1);   
                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }

                                            }
                                        
                                        else {
                                        
resultGQReal(m,n) += _quadweights(j)  * (MMx(m,mu) * MMx(n,nu) * DielInvGQ(mu,nu) / 
        (std::pow(_qpenergies(m) - _qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        CooTfFreqs(j) + MMx(m,mu) * MMx(n,nu) * DielInvGQ(mu,nu) / 
        (std::pow(_qpenergies(n) - _qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        CooTfFreqs(j));
resultGQImag(m,n) += _quadweights(j) * (MMx(m,mu) * MMx(n,nu) * DielInvGQ(mu,nu) / 
        (std::pow(_qpenergies(m)-_qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        (_qpenergies(m) - _qpenergies(k)) + MMx(m,mu) * MMx(n,nu) * DielInvGQ(mu,nu) / 
        (std::pow(_qpenergies(n)-_qpenergies(k),2)+std::pow(CooTfFreqs(j),2)) * 
        (_qpenergies(n) - _qpenergies(k)));
                                            
 if ( _qpenergies(k) < 0 ) {
                                                    if ( _qpenergies(m) < _qpenergies(k) ) {
                                                
resultRes(m,n) -= MMx(m,mu) * MMx(n,nu) * (DielInvResM(mu,nu) - 1);

                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }
                                                else {
                                                    if ( _qpenergies(m) > _qpenergies(k) ) {
resultRes(m,n) -= MMx(m,mu) * MMx(n,nu) * (DielInvResM(mu,nu) - 1);   
                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }

                                                if ( _qpenergies(k) < 0 ) {
                                                    if ( _qpenergies(n) < _qpenergies(k) ) {
                                                
resultRes(m,n) -= MMx(m,mu) * MMx(n,nu) * (DielInvResN(mu,nu) - 1);

                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }
                                                else {
                                                    if ( _qpenergies(n) > _qpenergies(k) ) {
resultRes(m,n) -= MMx(m,mu) * MMx(n,nu) * (DielInvResN(mu,nu) - 1);   
                                                    }
                                                    else {resultRes(m,n)+=0;}
                                                }
                                            

                                            }
                                        
                                            }
                                        
                                        }
            
                                    }
                        
                            }
                    
                        }
                
                    }
            
            }
        
        resultGQ.real()= resultGQReal;
        resultGQ.imag()= resultGQImag;
                
        return resultRes;//(resultGQ+resultRes).diagonal();
    }
    
    
    }
  
}
