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

#include <votca/xtp/rpa.h>
#include <votca/xtp/aomatrix.h>
#include "votca/xtp/threecenter.h"

namespace votca {
  namespace xtp {
      
  Eigen::VectorXd RPA::UpdateRPAInput(const Eigen::VectorXd& dftenergies,const Eigen::VectorXd& gwaenergies,int qpmin, int homo){
        int dftsize=dftenergies.size();
        Eigen::VectorXd rpaenergies=dftenergies;
        int gwsize=gwaenergies.size();
        int lumo=homo+1;

        int qpmax=qpmin+gwsize-1;
        rpaenergies.segment(qpmin,gwsize)=gwaenergies;
        double DFTgap = dftenergies(lumo) - dftenergies(homo);
        double QPgap = gwaenergies(lumo-qpmin) - gwaenergies(homo-qpmin);
        double shift=QPgap - DFTgap;
        rpaenergies.segment(qpmax+1,dftsize-qpmax-1).array()+=shift;
        return rpaenergies;
    }

 template< bool imag>
    Eigen::MatrixXd RPA::calculate_epsilon(double frequency)const{
        const int size = _Mmn.auxsize(); // size of gwbasis
        Eigen::MatrixXd result = Eigen::MatrixXd::Identity(size, size);
        const int lumo = _homo + 1;
        const int n_occ = lumo - _rpamin;
        const int n_unocc = _rpamax - _homo;
        const double freq2 = frequency*frequency;

#pragma omp parallel for
        for (int m_level = 0; m_level < n_occ; m_level++)        {
            const double qp_energy_m = _energies(m_level + _rpamin);
#if (GWBSE_DOUBLE)
            const Eigen::MatrixXd Mmn_RPA = _Mmn[ m_level].block(n_occ, 0, n_unocc, size);
#else
            const Eigen::MatrixXd Mmn_RPA = _Mmn[ m_level].block(n_occ, 0, n_unocc, size).cast<double>();
#endif
            const Eigen::ArrayXd deltaE=_energies.segment(lumo,n_unocc).array()-qp_energy_m;
            Eigen::VectorXd denom;
            if (imag){
                denom=4*deltaE/(deltaE.square()+freq2);
            }else{
                denom=2.0*((deltaE-frequency).inverse()+(deltaE+frequency).inverse());
            }
            auto temp=Mmn_RPA.transpose() *denom.asDiagonal();
            Eigen::MatrixXd tempresult = temp* Mmn_RPA;

#pragma omp critical
            {
                result += tempresult;
            }
        }
        return result;
    }


 template Eigen::MatrixXd RPA::calculate_epsilon<true>(double frequency)const;
 template Eigen::MatrixXd RPA::calculate_epsilon<false>(double frequency)const;

        rpa_eigensolution RPA::calculate_eigenvalues() const {

            Eigen::VectorXd AmB = calculate_spectral_AmB();
            Eigen::MatrixXd ApB = calculate_spectral_ApB();
            
            Eigen::MatrixXd C = calculate_spectral_C(AmB, ApB);
            
            return diag_C(AmB, C);
        }
        
        Eigen::VectorXd RPA::calculate_spectral_AmB() const {
            const int rpasize = (_homo - _rpamin + 1) * (_rpamax - (_homo + 1) + 1);
            Eigen::VectorXd AmB = Eigen::VectorXd::Zero(rpasize);
            
            // TODO: Symmetries
            
            for (int i = 0; i < rpasize; i++) {
                AmB(i) = _energies(_vc2index.c(i)) - _energies(_vc2index.v(i));
            } // Composite index i
            
            return AmB;
        }
        
        Eigen::MatrixXd RPA::calculate_spectral_ApB() const {
            const int rpasize = (_homo - _rpamin + 1) * (_rpamax - (_homo + 1) + 1);
            const int auxsize = _Mmn.auxsize(); // Size of gwbasis
            Eigen::MatrixXd ApB = Eigen::VectorXd::Zero(rpasize, rpasize);
            
            // TODO: Symmetries
            
            for (int i = 0; i < rpasize; i++) {
                ApB(i, i) = _energies(_vc2index.c(i)) - _energies(_vc2index.v(i));
            } // Composite index i

            for (int i_1 = 0; i_1 < rpasize; i_1++) {

                int v_1 = _vc2index.v(i_1);
                int c_1 = _vc2index.c(i_1);
                
                for (int i_2 = 0; i_2 < rpasize; i_2++) {
                    
                    int v_2 = _vc2index.v(i_2);
                    int c_2 = _vc2index.c(i_2);
                    
                    double fourcenter = 0.0;
                    
                    for (int i_aux = 0; i_aux < auxsize; ++i_aux) {
                        
                        VectorXfd tc_vc_1 = _Mmn[v_1].col(i_aux);
                        VectorXfd tc_vc_2 = _Mmn[v_2].col(i_aux);
                        
                        fourcenter += tc_vc_1(c_1) * tc_vc_2(c_2);
                        
                    } // Auxiliary basis function

                    ApB(i_1, i_2) -= 2 * fourcenter; // Fill (A + B)
                    
                } // Composite index i_2
            } // Composite index i_1
            
            return ApB;
        }
        
        Eigen::MatrixXd RPA::calculate_spectral_C(Eigen::VectorXd& AmB, Eigen::MatrixXd& ApB) const {
            
            return AmB.cwiseSqrt().asDiagonal() * ApB * AmB.cwiseSqrt().asDiagonal();
        }
        
        rpa_eigensolution RPA::diag_C(Eigen::VectorXd& AmB, Eigen::MatrixXd& C) const {
            const int size = _Mmn.auxsize(); // size of gwbasis
            
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);

            Eigen::VectorXd omega = es.eigenvalues().cwiseSqrt();
            Eigen::MatrixXd Z = es.eigenvectors();

            Eigen::VectorXd AmB_sqrt = AmB.cwiseSqrt();
            Eigen::VectorXd AmB_sqrt_inv = AmB_sqrt.cwiseInverse();
            Eigen::VectorXd Omega_sqrt = omega.cwiseSqrt();
            
            Eigen::MatrixXd XpY = Eigen::MatrixXd(size, size);

            for (int s = 0; s < size; s++) {

                Eigen::VectorXd lhs = (1 / Omega_sqrt(s)) * AmB_sqrt;
                Eigen::VectorXd rhs = (1 * Omega_sqrt(s)) * AmB_sqrt_inv;

                XpY.col(s) = 0.50 * (
                        (lhs + rhs).cwiseProduct(Z.col(s)) + // X
                        (lhs - rhs).cwiseProduct(Z.col(s))); // Y

            }

            // TODO: More efficient way to solve eq. 36 without sqrt is
            // described in section 6.2.
            
            rpa_eigensolution sol;
            sol._Omega = omega;
            sol._XpY = XpY;

            return sol;

        }

  }}
