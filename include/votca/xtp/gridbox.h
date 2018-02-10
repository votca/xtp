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

#ifndef __XTP_GRIDBOX__H
#define	__XTP_GRIDBOX__H

#include <votca/tools/vec.h>

#include <votca/tools/linalg.h>
#include <votca/xtp/grid_containers.h>
#include <votca/xtp/aoshell.h>

namespace votca { namespace xtp {

    namespace ub = boost::numeric::ublas;

        class GridBox {
            
        public: 
            
            
            
            const std::vector<tools::vec>& getGridPoints() const{return grid_pos;}
            
            const std::vector<double>& getGridWeights() const{return weights;}
            
            const std::vector<const AOShell* >& getShells() const{return significant_shells;}
            const std::vector<const AOShell* >& getShells_perMolecule() const{return mol_shells;}
            
            const std::vector<ub::range>& getAOranges() const{return aoranges;}
            const std::vector<ub::range>& getAOranges_perMolecule() const{return mol_aoranges;}
            
            unsigned size() const{return grid_pos.size();}
            
            unsigned Shellsize() const{return significant_shells.size();}
            
            unsigned Matrixsize() const{return matrix_size;}
            unsigned Matrixsize_perMolecule() const{return mol_matrix_size;}
            
            
            
            void addGridBox(const GridBox& box){
                const std::vector<tools::vec>& p=box.getGridPoints();
                const std::vector<double>& w=box.getGridWeights();
                for (unsigned i=0;i<w.size();++i){
                    grid_pos.push_back(p[i]);
                    weights.push_back(w[i]);
                }
                return;
            }

            void addGridPoint(const GridContainers::integration_grid& point){
                grid_pos.push_back(point.grid_pos);
                weights.push_back(point.grid_weight);
            };
            
            void addShell(const AOShell* shell){
              significant_shells.push_back(shell);  
              matrix_size+=shell->getNumFunc();
            };
            
            void prepareDensity(){
                densities.reserve(grid_pos.size());
            }
            
            void addDensity(double density){
                densities.push_back(density);
            }
            
            const std::vector<double>& getGridDensities() const{return densities;}
            
            void PrepareForIntegration();
            void PrepareForIntegration_perMolecule(std::vector<unsigned> relevant_atomids, std::vector<AOShell*> central_cell_shells);
            
            ub::matrix<double> ReadFromBigMatrix(const ub::matrix<double>& bigmatrix);
            ub::matrix<double> ReadFromBigMatrix_perMolecule(const ub::matrix<double>& bigmatrix);
            
            void AddtoBigMatrix(ub::matrix<double>& bigmatrix,const ub::matrix<double>& smallmatrix);
            void AddtoBigMatrix_perMolecule(ub::matrix<double>& bigmatrix, const ub::matrix<double>& smallmatrix);
            
            void setIndexoffirstgridpoint(unsigned indexoffirstgridpoint){_indexoffirstgridpoint=indexoffirstgridpoint;}
            unsigned getIndexoffirstgridpoint() const{return _indexoffirstgridpoint;}
            
			static bool compareGridboxes(GridBox& box1,GridBox& box2 ){
                if(box1.Matrixsize()!=box2.Matrixsize()){
                    return false;
                }
                if(box1.Shellsize()!=box2.Shellsize()){
                    return false;
                }
                for(unsigned i=0;i<box1.significant_shells.size();++i){
                    if(box1.significant_shells[i]!=box2.significant_shells[i]){
                        return false;
                    }
                    
                }
                return true;
            }


            void appendBoxData(GridBox& other){
                grid_pos.insert(grid_pos.end() , other.getGridPoints().begin(), other.getGridPoints().end());
                weights.insert(weights.end() , other.getGridWeights().begin(), other.getGridWeights().end());
                densities.insert(densities.end() , other.getGridDensities().begin(), other.getGridDensities().end());
            };
            
            void wrapPositions(tools::vec boxLen){
                std::vector< tools::vec >::iterator ip;
                for(ip=grid_pos.begin(); ip!=grid_pos.end(); ++ip) {
                    for(int k=0; k<3; k++)
                    {
                        (*ip)[k] = fmod( (*ip)[k], boxLen[k]);
                    }
                }
            }
            
            
        private:
            
                unsigned _indexoffirstgridpoint;
                unsigned matrix_size=0;
                std::vector<ub::range> aoranges;
                std::vector<ub::range> ranges;
                std::vector<ub::range> inv_ranges;
                std::vector< tools::vec > grid_pos;
                std::vector<const AOShell* > significant_shells;
                std::vector< double > weights;
                std::vector< double > densities;
                std::vector< ub::matrix<double> > dens_grad;
                
                //these are the ranges associated only with the shells of the target molecule.
                unsigned mol_matrix_size;
                std::vector<ub::range> mol_aoranges;
                std::vector<ub::range> mol_ranges;
                std::vector<ub::range> mol_inv_ranges;
                std::vector<const AOShell* >mol_shells;
            };

       

    }}
#endif	/* NUMERICAL_INTEGRATION_H */
