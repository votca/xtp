/*
 *            Copyright 2009-2017 The VOTCA Development Team
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
#include <votca/xtp/basisset.h>
#include <votca/xtp/threecenter.h>

namespace votca {
    namespace xtp {

        class RPA {
        public:

            void configure(unsigned homo, unsigned rpamin, unsigned rpamax) {
                _homo = homo;
                _rpamin = rpamin;
                _rpamax = rpamax;
            }
            
            void setScreening(Eigen::VectorXd& _screen_freq_r,Eigen::VectorXd& _screen_freq_i){
                screen_freq_r = _screen_freq_r;
                screen_freq_i = _screen_freq_i;
                _epsilon_r.resize(screen_freq_r.size());
                _epsilon_i.resize(screen_freq_i.size());
                
            }

            const Eigen::VectorXd& GetScreening_freq_r() const {
                return screen_freq_r;
            }
            
            const Eigen::VectorXd& GetScreening_freq_i() const {
                return screen_freq_i;
            }

            const std::vector<Eigen::MatrixXd>& GetEpsilon_r() const {
                return _epsilon_r;
            }
            
            const std::vector<Eigen::MatrixXd>& GetEpsilon_i() const {
                return _epsilon_i;
            }

            void prepare_threecenters(const TCMatrix_gwbse& _Mmn_full);

            void calculate_epsilon(const Eigen::VectorXd& qp_energies);

            void FreeMatrices() {
                _Mmn_RPA.Cleanup();
                for (Eigen::MatrixXd & matrix:_epsilon_r){
                    matrix.resize(0,0);
                }
                for (Eigen::MatrixXd & matrix:_epsilon_i){
                    matrix.resize(0,0);
                }
            
            }

        private:
            TCMatrix_gwbse _Mmn_RPA;

            unsigned _homo; // HOMO index
            unsigned _rpamin;
            unsigned _rpamax;
            double _shift; // pre-shift of DFT energies

            // container for the epsilon matrix
            std::vector<Eigen::MatrixXd > _epsilon_r;
            
            std::vector<Eigen::MatrixXd > _epsilon_i;
           
            
            // We cannot calculate screening at complex frequencies only at real or imaginary points
            Eigen::VectorXd screen_freq_r; //real screening frequencies
            Eigen::VectorXd screen_freq_i;//imaginary screening frequencies

           


        };
    }
}

#endif /* _VOTCA_RPA_RPA_H */
