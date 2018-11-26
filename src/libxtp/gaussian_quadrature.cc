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
namespace votca {
  namespace xtp {
GaussianQuadrature::GaussianQuadrature(const Eigen::VectorXd& qpenergies,const TCMatrix_gwbse& Mmn):_qpenergies(qpenergies),_Mmn(Mmn){
    //initialise weights and points
    _integrationpoints=Eigen::VectorXd::Zero(_order);
    _integrationweights=Eigen::VectorXd::Zero(_order);
    //fill in weights and points
    _integrationpoints << 0.9815606342467191, 0.9041172563704749, 0.7699026741943047,
  0.5873179542866174, 0.3678314989981801, 0.1252334085114689, -0.1252334085114689,
  0.3678314989981801, -0.5873179542866174, -0.7699026741943047, -0.9041172563704749,
  -0.9815606342467191;
    _integrationweights << 0.04717533638651180, 0.1069393259953181, 0.1600783285433461,
  0.2031674267230658, 0.2334925365383548, 0.2491470458134029, 0.2491470458134029,
  0.2334925365383548, 0.2031674267230658, 0.1600783285433461, 0.1069393259953181, 
  0.04717533638651180;
}


Eigen::MatrixXcd GaussianQuadrature::Integrate(RPA& rpa){
    Eigen::VectorXd qpenergies;
    int AmountOfGaussianBasisFunctions = qpenergies.size();
    const std::complex<double> I(0.0,1.0);       
    //initialise the resulting matrix
    Eigen::MatrixXcd result;
    result=Eigen::MatrixXd::Zero(AmountOfGaussianBasisFunctions,AmountOfGaussianBasisFunctions);
    
   //fill the frequency vector
    Eigen::VectorXd frequencies;
    for (int k = 0 ; k < _order ; ++k){
        frequencies(k) = std::log((1+_integrationpoints(k))/(1-_integrationpoints(k)));
    }
    
    std::cout<<_Mmn[0](0,0)<<std::endl;
   //fill the matrix
    for (int n = 0; n < AmountOfGaussianBasisFunctions; ++n){
        for (int m = 0; m < AmountOfGaussianBasisFunctions; ++m){
            for (int j = 0; j < _order ; ++j){
                Eigen::MatrixXd DielectricMatrixInverse 
                = CalcInverse(frequencies(j),rpa);
                int AmountOfAuxiliaryBasisFunctions = DielectricMatrixInverse.rows();
                for (int mu = 0; mu < AmountOfAuxiliaryBasisFunctions; ++mu){
                    for (int nu = 0; nu < AmountOfAuxiliaryBasisFunctions; ++nu){
                        for (int k = 0; k < AmountOfGaussianBasisFunctions; ++k){
                            //First, we consider the off-diagonal terms of the self energy matrix
                           if (m != n){
                                result(n,m)+=0.5*(_Mmn[k](m,mu)*_Mmn[k](n,mu)*_integrationweights(j)*(2.0*I/(1-std::pow(_integrationpoints(j),2)))*
                                ((DielectricMatrixInverse(mu,mu)-1)/(std::pow(qpenergies(m)-qpenergies(k),2)+pow(frequencies(j),2)))
                                *(frequencies(j)+I*(qpenergies(m)-qpenergies(k)))+
                                _Mmn[k](m,mu)*_Mmn[k](n,mu)*_integrationweights(j)*(2.0*I/(1-std::pow(_integrationpoints(j),2)))*
                                ((DielectricMatrixInverse(mu,mu)-1)/(std::pow(qpenergies(n)-qpenergies(k),2)+pow(frequencies(j),2)))
                                *(frequencies(j)+I*(qpenergies(n)-qpenergies(k)))); 
                                if (mu != nu){
                                result(n,m)+=0.5*(_Mmn[k](m,mu)*_Mmn[k](n,nu)*_integrationweights(j)*(2.0*I/(1-std::pow(_integrationpoints(j),2)))*
                                (DielectricMatrixInverse(mu,nu)/(std::pow(qpenergies(m)-qpenergies(k),2)+pow(frequencies(j),2)))
                                *(frequencies(j)+I*(qpenergies(m)-qpenergies(k)))+
                                _Mmn[k](m,mu)*_Mmn[k](n,nu)*_integrationweights(j)*(2.0*I/(1-std::pow(_integrationpoints(j),2)))*
                                (DielectricMatrixInverse(mu,nu)/(std::pow(qpenergies(n)-qpenergies(k),2)+pow(frequencies(j),2)))
                                *(frequencies(j)+I*(qpenergies(n)-qpenergies(k)))); 
                                }
                                }
                                    // Now, m = n (diagonal terms of the self-energy matrix)
                                 if (mu != nu){
                                     //contributing terms without Kronecker delta, when mu != nu and m = n
                                result(m,m)+=_Mmn[k](m,mu)*_Mmn[k](m,nu)*_integrationweights(j)*(2.0*I/(1-std::pow(_integrationpoints(j),2)))*
                                 (DielectricMatrixInverse(mu,nu)/(std::pow(qpenergies(m)-qpenergies(k),2)+pow(frequencies(j),2)))
                                 *(frequencies(j)+I*(qpenergies(m)-qpenergies(k))); 
                                }
                                 //contributing term with Kronecker delta, when mu = nu and m = n
                                 result(m,m)+=_Mmn[k](m,mu)*_Mmn[k](m,mu)*_integrationweights(j)*(2.0*I/(1-std::pow(_integrationpoints(j),2)))*
                                 ((DielectricMatrixInverse(mu,mu)-1)/(std::pow(qpenergies(m)-qpenergies(k),2)+pow(frequencies(j),2)))
                                 *(frequencies(j)+I*(qpenergies(m)-qpenergies(k)));
                                 
                            }                  
                        }
                    }
                }
            }
        }
 return result;    
}

Eigen::MatrixXd GaussianQuadrature::CalcInverse(double frequency, RPA& rpa){
    Eigen::VectorXd real=Eigen::VectorXd::Zero(0);
    Eigen::VectorXd imag=Eigen::VectorXd::Ones(1)*frequency;
    rpa.setScreening(real,imag);
    rpa.calculate_epsilon(_qpenergies,_Mmn);
    const Eigen::MatrixXd& dielectric_matrix=rpa.GetEpsilon_i()[0];
    return dielectric_matrix.inverse();
}
    
  }
}
