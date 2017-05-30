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

#ifndef BULKESP_H
#define BULKESP_H

#include <votca/xtp/espfit.h>
#include <votca/xtp/orbitals.h>



using namespace votca::tools;

namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;

        /**
         * \brief Wraps a vector to periodic boundary conditions
         * 
         * @param vector r to wrap
         * @param vector box specifying the periodic box
         * @return the wrapped vector
         */
        tools::vec Wrap(const tools::vec r, const tools::vec box);

        class Bulkesp : public Espfit {
        public:

            struct Bond {
                ctp::QMAtom* a;
                ctp::QMAtom* b;
                tools::vec ba;
                int a_indx, b_indx;
            };

            struct Molecule {
                std::vector< ctp::QMAtom* > atoms;
                std::vector< int > atomIndeces;
            };


        public:

            Bulkesp(ctp::Logger *log) : Espfit(log) {
                periodic = false;
                boxLen[0] = 0;
                boxLen[1] = 0;
                boxLen[2] = 0;

                dipolesLog = new std::ofstream();
            }

            ~Bulkesp() {
                delete dipolesLog;
            }

            void setBox(double b[3]) {
                periodic = true;
                boxLen[0] = b[0];
                boxLen[1] = b[1];
                boxLen[2] = b[2];
            }

            std::vector<Bulkesp::Molecule> BreakIntoMolecules(std::vector< ctp::QMAtom* > a, double scale);

            ub::vector<double> ComputeESP(std::vector< ctp::QMAtom* > & _global_atomlist,
                    std::vector< ctp::QMAtom* > & _local_atomlist, std::vector<int> _local_atomIndeces,
                    ub::matrix<double> &_global_dmat, AOBasis &_global_basis, BasisSet &bs, string gridsize, Grid &_grid, double &netcharge);

            void Evaluate(std::vector< ctp::QMAtom* >& _atomlist, ub::matrix<double> &_global_dmat, Orbitals& _orbitals, ub::matrix<double> _global_MO_Coeffs,
                    AOBasis &_basis, BasisSet &bs, string gridsize, double maxBondScale,
                    std::string _state, std::string _spin, int _state_no);

            void FillElement2NBF(std::vector< ctp::QMAtom* >& _atomlist, BasisSet &bs);

            std::vector<double> FitPartialCharges(std::vector< tools::vec >& _fitcenters, Grid& _grid, ub::vector<double>& _potential, double& _netcharge);

        private:

            std::map<std::string, int> _element2NBF; //Number of Basis Functions for each element in the basis set
            list<std::string> _elements; //list of all elements in the QMatoms vector
            std::string fn_prefix; //prefix for output files containing potentials
            bool periodic; //is the box periodic for the purposes of assigning atoms to molecules?
            double boxLen[3]; //dimensions of the box, assume cuboid shape
            std::ofstream* dipolesLog; //file to log dipoles of all the molecules

            std::map<ctp::QMAtom*, int> MapAtom2MOCoefIndex(std::vector< ctp::QMAtom* >& _atomlist);

            ub::matrix<double> BuildDenMat(Orbitals &_orb, std::string _state, std::string _spin, int _state_no);

            ub::matrix<double> BuildOverlapMat(Orbitals &_molOrb, Orbitals &_globalOrb);
        };
}}




#endif /* BULKESP_H */

