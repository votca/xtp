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

#ifndef _VOTCA_XTP_RPA_H
#define _VOTCA_XTP_RPA_H

#include <votca/xtp/eigen.h>
#include <vector>
#include <votca/xtp/vc2index.h>
#include <votca/xtp/rpa_eigensolution.h>

namespace votca
{
namespace xtp
{
class TCMatrix_gwbse;

class RPA
{
public:

    RPA(const Eigen::VectorXd& energies, const TCMatrix_gwbse& Mmn):
        _energies(energies),_Mmn(Mmn),_vc2index(0, 0, 0){}; // TODO: Pass vmin, cmin, ctotal as input?

    void configure(int homo, int rpamin, int rpamax){
        _homo = homo;
        _rpamin = rpamin;
        _rpamax = rpamax;
        _vc2index = vc2index(rpamin, homo + 1, rpamax - (homo + 1));
    }

    Eigen::MatrixXd calculate_epsilon_i(double frequency)const{
        return calculate_epsilon<true>(frequency);
    }

    Eigen::MatrixXd calculate_epsilon_r(double frequency)const{
        return calculate_epsilon<false>(frequency);
    }
    
    rpa_eigensolution calculate_eigenvalues();

    //calculates full RPA vector of energies from gwa and dftenergies and qpmin
    //RPA energies have three parts, lower than qpmin: dftenergies,between qpmin and qpmax:gwa_energies,above:dftenergies+homo-lumo shift
    static Eigen::VectorXd UpdateRPAInput(const Eigen::VectorXd& dftenergies,const Eigen::VectorXd& gwaenergies,int qpmin, int homo);

private:

    int _homo; // HOMO index
    int _rpamin;
    int _rpamax;

    const Eigen::VectorXd& _energies;
    const TCMatrix_gwbse& _Mmn;

    template< bool imag>
    Eigen::MatrixXd calculate_epsilon(double frequency)const;
    
    vc2index _vc2index;
    
    Eigen::VectorXd calculate_spectral_AmB();
    Eigen::MatrixXd calculate_spectral_ApB();
    Eigen::MatrixXd calculate_spectral_C(Eigen::VectorXd& AmB, Eigen::MatrixXd& ApB);
    rpa_eigensolution diag_C(Eigen::VectorXd& AmB, Eigen::MatrixXd& C);

};
}
}

#endif /* _VOTCA_RPA_RPA_H */
