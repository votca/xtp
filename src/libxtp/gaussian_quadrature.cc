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


Eigen::VectorXd GaussianQuadrature::TranslateFrequencies() {
    Eigen::VectorXd result(_order);
       for (int k = 0 ; k < _order ; ++k){
        result(k) = std::log((1+_integrationpoints(k))/(1-_integrationpoints(k)));
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

Eigen::MatrixXd GaussianQuadrature::SumInversesMinusOne(RPA& rpa){
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(numofqplevels,numofqplevels);
    Eigen::VectorXd TranslatedFrequencies = TranslateFrequencies();
    for (int k = 0 ; k < _order ; ++k){
        result+=CalcInverse(TranslatedFrequencies(k),rpa);
    }
    result-=Eigen::MatrixXd::Identity(numofqplevels,numofqplevels);
    return result;
}

Eigen::MatrixXcd GaussianQuadrature::Integrate(RPA& rpa){
const std::complex<double> I(0.0,1.0);
Eigen::MatrixXd resultRealDiagonal=Eigen::MatrixXd::Zero(numofqplevels,numofqplevels);
Eigen::MatrixXd resultImagDiagonal=Eigen::MatrixXd::Zero(numofqplevels,numofqplevels);
Eigen::MatrixXd resultRealOffDiagonal=Eigen::MatrixXd::Zero(numofqplevels,numofqplevels);
Eigen::MatrixXd resultImagOffDiagonal=Eigen::MatrixXd::Zero(numofqplevels,numofqplevels);
Eigen::VectorXd TranslatedFrequencies = TranslateFrequencies();
Eigen::MatrixXd SummedInversesMinusOne=SumInversesMinusOne(rpa);

for (int k = 0 ; k < numofqplevels ; ++k){
    Eigen::MatrixXd AdaptedFrequencyReal=Eigen::MatrixXd::Zero(_order,numofqplevels);
    Eigen::MatrixXd AdaptedFrequencyImag=Eigen::MatrixXd::Zero(_order,numofqplevels);
    #if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& MMatrix = _Mmn[  k ];
#else
                const Eigen::MatrixXd MMatrix = _Mmn[  k ].cast<double>();       
#endif
                
    for (int m = 0 ; m < numofqplevels ; ++m){
                for (int j = 0 ; j < _order ; ++j){
                    AdaptedFrequencyReal(j,m) = TranslatedFrequencies(j)/(std::pow(_qpenergies(m)-_qpenergies(k),2)+std::pow(TranslatedFrequencies(j),2));
                    AdaptedFrequencyImag(j,m) = (_qpenergies(m)-_qpenergies(k))/(std::pow(_qpenergies(m)-_qpenergies(k),2)+std::pow(TranslatedFrequencies(j),2));
                }
        }
resultRealDiagonal+=2*_integrationweights.transpose()*AdaptedFrequencyReal*((MMatrix*SummedInversesMinusOne*MMatrix.transpose()).diagonal()).asDiagonal();
resultImagDiagonal+=2*_integrationweights.transpose()*AdaptedFrequencyImag*((MMatrix*SummedInversesMinusOne*MMatrix.transpose()).diagonal()).asDiagonal();
resultRealOffDiagonal+=(_integrationweights.transpose()*AdaptedFrequencyReal).asDiagonal()*(MMatrix*(SummedInversesMinusOne+SummedInversesMinusOne.transpose())*MMatrix.transpose());
resultImagOffDiagonal+=(_integrationweights.transpose()*AdaptedFrequencyImag).asDiagonal()*(MMatrix*(SummedInversesMinusOne+SummedInversesMinusOne.transpose())*MMatrix.transpose());
}
resultRealDiagonal=resultRealDiagonal.asDiagonal();
resultImagDiagonal=resultImagDiagonal.asDiagonal();
resultRealOffDiagonal-= (resultRealOffDiagonal.diagonal()).asDiagonal();
resultImagOffDiagonal-= (resultImagOffDiagonal.diagonal()).asDiagonal();
return resultRealDiagonal+resultRealOffDiagonal+I*(resultImagOffDiagonal+resultImagDiagonal);
}
}
}
