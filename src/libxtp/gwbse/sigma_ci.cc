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

    //This function calculates the inverse of the microscopic dielectric
    //matrix for given (real or purely imaginary) frequency and Kohn-Sham energies
          
    void Sigma_CI::PrepareScreening(){
         
        }
  
    //This function returns the whole self-energy expectation matrix
    //for given Kohn-Sham energies and M coefficients (hence "RPA")

    Eigen::MatrixXd Sigma_CI::CalcCorrelationOffDiag(const Eigen::VectorXd&
        frequencies)const{
        
        const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
        
    
        //First, we initialise the result, which is a sum, at zero
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal,_qptotal);
    
        //Then, using a triple for loop, we fill the matrix. We use symmetry,
        //so the second index runs up and including the first index, in order to
        //get the diagonal too.
            for (int k = 0 ; k < _qptotal ; ++k){
                
                //making sure the M tensor is double-valued
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMatrix = _Mmn[k];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
                    #endif
                
                if (energies(k) > 0){
            
                    for (int m = 0 ; m < _qptotal ; ++m){
                            
                        if ( frequencies(m) > energies(k) ){
                        //evaluating the dielectric matrix inverse at the right
                        //(real) frequency value
                        Eigen::MatrixXd DielInvMinId =
                            _rpa.calculate_epsilon_r(energies(k)
                                - frequencies(m)).inverse();
                        
                        //now, we subtract the identity
                        int nobasisfcs = DielInvMinId.size();        
                        DielInvMinId -=
                            Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                    
                        //now, we compute the matrix multiplication, which is 
                        //the contribution for this part of the loop
                        Eigen::MatrixXd ResPart = 
                            MMx*DielInvMinId*MMx.transpose();
                            
                        for (int n = 0 ; n < _qptotal ; ++n){
                        
                            //the resulting matrix entry is the sum of all these
                            result(m,n) += ResPart(m,n);
                    
                            }
                        
                            }
                        } 
                    
                    for (int n = 0 ; n < _qptotal ; ++n){
                        
                        if ( frequencies(n) > energies(k) ){
                
                        Eigen::MatrixXd DielInvMinId =
                            _rpa.calculate_epsilon_r(energies(k)
                                - frequencies(n)).inverse();
                        int nobasisfcs = DielInvMinId.size();        
                        DielInvMinId -= 
                            Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                        Eigen::MatrixXd ResPart = 
                            MMx * DielInvMinId * MMx.transpose();
                        
                            for (int m = 0 ; m < _qptotal ; ++m){
                        
                            //the resulting matrix entry is the sum of all these
                            result(m,n) += ResPart(m,n);
                    
                            }
                        
                        }
                    
                    }
                    }
                    
                    else  if (energies(k) < 0){
                        
                                for (int m = 0 ; m < _qptotal ; ++m){
                                      
                                    if ( frequencies(m) < energies(k) ){
                                        
                                    Eigen::MatrixXd DielInvMinId =
                                        _rpa.calculate_epsilon_r(energies(k)
                                            - energies(m)).inverse();
                                    int nobasisfcs=DielInvMinId.size();        
                                    DielInvMinId -=
                                        Eigen::MatrixXd::Identity(nobasisfcs,
                                        nobasisfcs);
                                    Eigen::MatrixXd ResPart =
                                        MMx * DielInvMinId * MMx.transpose();
                                        
                                        
                                        for (int n = 0 ; n < _qptotal ; ++n){
                        
                                            //the resulting matrix entry is the 
                                            //sum of all these
                                            result(m,n) += ResPart(m,n);
                    
                                            }
                        
                                    }
                                }
                                
                                for (int n = 0 ; n < _qptotal ; ++n){
                                
                                    Eigen::MatrixXd DielInvMinId = 
                                        _rpa.calculate_epsilon_r(energies(k)
                                            - energies(n)).inverse();
                                    int nobasisfcs = DielInvMinId.size();        
                                    DielInvMinId -= 
                                        Eigen::MatrixXd::Identity(nobasisfcs,
                                        nobasisfcs);
                                    Eigen::MatrixXd ResPart =
                                        MMx * DielInvMinId * MMx.transpose();
                            
                                        for (int m = 0 ; m < _qptotal ; ++n){
                        
                                            //the resulting matrix entry is
                                            //the sum of all these
                                            result(m,n) += ResPart(m,n);    
                                            
                                            }
                                
                                    }
                                
                                }
                            
                            }
                    
                
                
        //Now, we just add the transpose matrix, to get the full, symmetric,
        //matrix
        //result+=result.transpose();
        //There is a factor (-1) to be put in front, and we obtain the result
        return (-1)*result+_gq.SigmaGQ(frequencies,_rpa);
    
        }
    
    //This function returns the diagonal of aforementioned matrix SigmaRes
    //in vector form
    
    Eigen::VectorXd Sigma_CI::CalcCorrelationDiag(const Eigen::VectorXd& 
    frequencies)const{
        
        const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
        

        //First, we initialise the result, which is a sum, at zero
        Eigen::VectorXd result=Eigen::VectorXd::Zero(_qptotal,_qptotal);
            
            //Now, everything is the same, but now we only have one sum to split
            //since m = n on the diagonal
            for (int k = 0 ; k < _qptotal ; ++k){
                
                    #if (GWBSE_DOUBLE)
                        const Eigen::MatrixXd& MMatrix = _Mmn[k];
                    #else
                        const Eigen::MatrixXd MMx = _Mmn[k].cast<double>();       
                    #endif
                
                if (energies(k) > 0){
                    
                        for (int m = 0 ; m < _qptotal ; ++m){
                            
                            if ( frequencies(m) > energies(k) ){
                
                            Eigen::MatrixXd DielInvMinId = 
                                _rpa.calculate_epsilon_r(energies(k) -
                                    frequencies(m)).inverse();
                            int nobasisfcs = DielInvMinId.size();        
                            DielInvMinId -= 
                                Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                            Eigen::MatrixXd ResPart = 
                                MMx * DielInvMinId * MMx.transpose();
                            result(m) += ResPart(m,m);
                        
                            }
                        
                        }
            
                    }
                
                else if (energies(k) < 0){
                    
                        for (int m = 0 ; m < _qptotal ; ++m){
                            
                            if ( frequencies(m) < energies(k) ){
                           
                            Eigen::MatrixXd DielInvMinId = 
                                _rpa.calculate_epsilon_r(energies(k)
                                    - frequencies(m)).inverse();
                            int nobasisfcs = DielInvMinId.size();        
                            DielInvMinId -= 
                                Eigen::MatrixXd::Identity(nobasisfcs,nobasisfcs);
                            Eigen::MatrixXd ResPart = 
                                MMx * DielInvMinId * MMx.transpose();
                            result(m) += ResPart(m,m);
                
                            }
                
                        }
        
                    }
                        
            }
                
        //Because m = n on the diagonal, we also get a double factor in front
        return (-2)*result+_gq.SigmaGQDiag(frequencies,_rpa);
        
        }
    
   
  }

}

