/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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

#ifndef __XTP_NUMERICAL_INTEGRATION__H
#define	__XTP_NUMERICAL_INTEGRATION__H

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/xtp/basisset.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/grid_containers.h>

#include <votca/ctp/qmatom.h>


namespace votca { namespace xtp {

    namespace ub = boost::numeric::ublas;
    
    

        class NumericalIntegration {
        public: 
            
            NumericalIntegration():density_set(false) {};

            void GridSetup(std::string type, BasisSet* bs , std::vector<ctp::QMAtom* > _atoms,AOBasis* basis  );
            
            //void FindsignificantAtoms2(AOBasis* basis);
            //used for test purposes
            double StupidIntegrate( std::vector<double>& _data );
            
            std::vector<const vec*> getGridpoints();
            
            
            
            double IntegrateDensity_Atomblock(const ub::matrix<double>& _density_matrix, AOBasis* basis);
            double IntegrateDensity_Molecule(ub::matrix<double>& _density_matrix, AOBasis* basis, std::vector<int> AtomIndeces);
            double IntegratePotential(const vec& rvector);
            
            double IntegrateField(const std::vector<double>& externalfield);

            double IntegratePotential_w_PBC(vec rvector, vec boxLen);
            void IntegratePotential_w_PBC_gromacs_like(Grid &eval_grid, vec boxLen, ub::vector<double>& _ESPatGrid);
            double IntegrateEnergy_w_PBC(vec rvector, vec boxLen);
            double CalcDipole_w_PBC(vec rvector, vec boxLen);
            void findAlpha(double Rc, double dtol);
            void PrepKspaceDensity(vec boxLen, double ext_alpha, std::vector< ctp::QMAtom* > & _local_atomlist, bool ECP, int nK);
            void PrepKspaceDensity_gromacs_like(vec boxLen, double ext_alpha, std::vector< ctp::QMAtom* > & _local_atomlist, bool ECP, Grid &eval_grid, int nK);
            void FreeKspace(void);
            std::vector< std::vector< GridContainers::integration_grid > > _Madelung_grid;
            void FillMadelungGrid(vec boxLen, int natomsonside);

            double getExactExchange(const std::string _functional);
            // in principle a symmetric matrix would be nicer but we calculate whole vxc matrix because of numerics and symmetrize explicitly 
            ub::matrix<double> IntegrateVXC_Atomblock (const ub::matrix<double>& _density_matrix,const std::string _functional);
            //ub::matrix<double> IntegrateVXC_Atomblock2 (const ub::matrix<double>& _density_matrix, AOBasis* basis,const std::string _functional);
            ub::matrix<double> IntegrateExternalPotential_Atomblock(const std::vector<double>& Potentialvalues);
         
            
            // this gives int (e_xc-V_xc)*rho d3r
            double getTotEcontribution(){return EXC;}
          
            
        private:
            
            AOBasis* _basis;
            
           void FindsignificantAtoms();
           double erf1c(double x);
            double erfcc(double x);
            std::vector<double> SSWpartition(int igrid, int ncenters ,  std::vector< std::vector<double> >& rq );
            
            
            std::vector<double> Rij;

            double  _totalgridsize;
            std::vector< std::vector< GridContainers::integration_grid > > _grid;
            double EXC;
            bool density_set;
            std::vector< std::vector< std::vector<int> > > _significant_atoms;
            std::vector < int > _startIdx;
            std::vector < int > _blocksize;
            std::vector< std::vector< AOBasis::AOShellIterator > > _atomshells;
            std::vector< AOBasis::AOShellIterator > _singleatom;
            std::vector< std::vector< ub::matrix<double> > >dmat_vector;
            std::vector< std::vector< std::vector< ub::matrix<double> > > > xcmat_vector_thread;
            std::vector< std::vector< ub::matrix<double> > > xcmat_vector;
           //vector< vector<int> > _atomsforshells;
        
                    
        public:
            std::complex<double>* Rho_k; //density in k-space, used for Ewald summation of potential in periodic systems
            //std::complex<double>**** eikR;  //gromacs-like storage for exp(k*R) -> where to evaluate
            //std::complex<double>*** eikr;  //gromacs-like storage for exp(k*r) -> charge distribution
            std::vector<std::vector<std::vector< std::complex<double> > > > eikR;  //gromacs-like storage for exp(k*R) -> where to evaluate
            std::vector<std::vector<std::vector< std::vector<std::complex<double> > > > > eikr;  //gromacs-like storage for exp(k*r) -> charge distribution
            vec lll;
            int numK[3];   //number of k-vectors along each axis
            double alpha;  //inverse length in Ewald summation
            double *Kcoord;//k-values
            double *prefactor; //prefactors for k-space components of potential before summation
            double E_rspace;
            double E_kspace;
            double E_erfc;
        };

    }}
#endif	/* NUMERICAL_INTEGRATION_H */
