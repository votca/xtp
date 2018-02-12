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

#ifndef __XTP_NUMERICAL_INTEGRATION_PERIODIC__H
#define	__XTP_NUMERICAL_INTEGRATION_PERIODIC__H


#include <votca/xtp/numerical_integrations.h>
#include <votca/xtp/grid.h>

namespace votca { namespace xtp {

    namespace ub = boost::numeric::ublas;
    
    class NumericalIntegrationPeriodic: public NumericalIntegration {
        public: 
            
            ~NumericalIntegrationPeriodic(){
                if(_expanded_basis!=_basis)
                    delete _expanded_basis;
                for(const auto& imgatom:_toclean_atoms){
                    delete imgatom;
                }
                
            };
            
            virtual void GridSetup(std::string type, BasisSet* bs , std::vector<ctp::QMAtom* > _atoms,AOBasis* basis  );
            void SetRelevantAtomIds(std::vector< unsigned >& relevantAtomIndeces){
                _relevant_atomids = relevantAtomIndeces;
                return;
            };
            
            virtual double IntegrateDensity(const ub::matrix<double>& _density_matrix);
            void IntegratePotential_w_PBC_gromacs_like(Grid &eval_grid, ub::vector<double>& _ESPatGrid);
            double CalcDipole_w_PBC(vec rvector);
            void findAlpha(double Rc, double dtol);
            void PrepKspaceDensity_gromacs_like(double ext_alpha, std::vector< ctp::QMAtom* > & _local_atomlist, bool ECP, Grid &eval_grid, int nK);
            inline void setBox(vec box){
                boxLen=box;
                return;
            }
            
        protected:
            virtual void FindSignificantShells();
            virtual void SortGridpointsintoBlocks(std::vector< std::vector< GridContainers::integration_grid > >& grid);
            std::vector<double> SSWpartition(int igrid, int ncenters, int nexpandedcenters,  std::vector< std::vector<double> >& rq );
                    
        private:
            void ExpandBasis(vector<ctp::QMAtom*> _atoms);
            void FindCenterCenterDist(vector<ctp::QMAtom*> _atoms);
            std::vector< std::vector<double> > FindGridpointCenterDist(vector<ctp::QMAtom*> _atoms, std::vector< GridContainers::integration_grid > _atomgrid);
            
            std::complex<double>* Rho_k; //density in k-space, used for Ewald summation of potential in periodic systems
            std::vector< std::vector< std::vector< std::complex<double> > > > eikR;  //gromacs-like storage for exp(k*R) -> where to evaluate
            std::vector< std::vector< std::vector< std::complex<double> > > > eikr;  //gromacs-like storage for exp(k*r) -> charge distribution
            vec lll;
            vec boxLen;    //in Bohr
            int numK[3];   //number of k-vectors along each axis
            double alpha;  //inverse length in Ewald summation
            double E_rspace;
            double E_kspace;
            double E_erfc;
            GridBox _periodicGridBox;
            
            int _nExpantionCells; //expand basis and atoms to include atoms in this many periodic cells away
            AOBasis* _expanded_basis;
            vector<ctp::QMAtom*> _expanded_atoms;
            vector<ctp::QMAtom*> _toclean_atoms;
            std::vector<unsigned> _relevant_atomids;
            
            int nFuncInMol;
            std::vector<ub::range> global_mol_aoranges;
            std::vector<ub::range> global_mol_inv_aoranges;
            std::vector<AOShell*> central_cell_shells;
            
    };
        
    tools::vec WrapPoint(const tools::vec r, const tools::vec box);
    tools::vec WrapDisplacement(const tools::vec a, const tools::vec b, const tools::vec box);

}}
#endif	/* NUMERICAL_INTEGRATION_PERIODIC__H */
