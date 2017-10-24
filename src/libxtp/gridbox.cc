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

#include <votca/xtp/gridbox.h>





namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
        void GridBox::AddtoBigMatrix(ub::matrix<double>& bigmatrix, const ub::matrix<double>& smallmatrix) {


            for (unsigned i = 0; i < ranges.size(); i++) {
                for (unsigned j = 0; j < ranges.size(); j++) {
                    ub::project(bigmatrix, ranges[i], ranges[j]) += ub::project(smallmatrix, inv_ranges[i], inv_ranges[j]);
                }
            }
            return;
        }

            

        ub::matrix<double> GridBox::ReadFromBigMatrix(const ub::matrix<double>& bigmatrix) {
            
            ub::matrix<double> _matrix = ub::zero_matrix<double>(matrix_size);
            for (unsigned i = 0; i < ranges.size(); i++) {
                for (unsigned j = 0; j < ranges.size(); j++) {
                    
                    ub::project(_matrix, inv_ranges[i], inv_ranges[j]) = ub::project(bigmatrix, ranges[i], ranges[j]);
                }
            }
            return _matrix;
        }

        void GridBox::PrepareForIntegration() {
            matrix_size = 0;

            std::vector<unsigned> start;
            std::vector<unsigned> end;

            for (unsigned i=0;i< significant_shells.size();++i) {
                const AOShell* shell=significant_shells[i];
                aoranges.push_back(ub::range(matrix_size, matrix_size+shell->getNumFunc()));
                matrix_size += shell->getNumFunc();
                start.push_back(shell->getStartIndex());
                end.push_back(shell->getStartIndex() + shell->getNumFunc());
            }
            std::vector<unsigned> startindex;
            std::vector<unsigned> endindex;
            
            if(start.size()>1){
                startindex.push_back(start[0]);

                for (unsigned i = 0; i < start.size() - 1; ++i) {

                    if(end[i]!=start[i+1]){
                        startindex.push_back(start[i+1]);
                        endindex.push_back(end[i]);
                    }
                }
                endindex.push_back(end[end.size()-1]);
            }
            else{
                startindex=start;
                endindex=end;
             }
            unsigned shellstart = 0;
            for (unsigned i = 0; i < startindex.size(); ++i) {
                ranges.push_back(ub::range(startindex[i], endindex[i]));
                
                unsigned size = endindex[i]-startindex[i];
                
                inv_ranges.push_back(ub::range(shellstart, shellstart + size));
                shellstart += size;
            }

            
            return;
        }
    
        
    
    //This is not only per molecule, but also used for periodic systems.
    //ReadFromBigMatrix_perMolecule needs to output all periodically expanded
    //signifficant shells in rows, and only the relevant shells from atoms in the first periodic cell.
    //For periodicity, just treat all atoms in the first periodic cell as the molecule.
    void GridBox::PrepareForIntegration_perMolecule(std::vector<unsigned> relevant_atomids, std::vector<AOShell*> central_cell_shells){
            //prepare the ranges for rows
            PrepareForIntegration();
            
            //and now do it for colums
            mol_matrix_size = 0;

            std::vector<unsigned> start;
            std::vector<unsigned> end;

            for (unsigned i=0;i< significant_shells.size();++i) {
                //Expanded shells have the same atom index as shells from first periodic cell.
                //So can check if a significant shell belongs to the molecule by comparing these atom indeces.
                const AOShell* shell=significant_shells[i];
                //skip any significant shells not in this molecule
                if(std::find(relevant_atomids.begin(), relevant_atomids.end(), shell->getIndex()) != relevant_atomids.end()){
                    //also check if the shell is from the first periodic box (unexpanded basis)
                    //don't want projections onto periodic images counted as part of this molecule's charge
                    if(std::find(central_cell_shells.begin(), central_cell_shells.end(), shell) != central_cell_shells.end())
                    {
                        mol_aoranges.push_back(ub::range(mol_matrix_size, mol_matrix_size+shell->getNumFunc()));
                        mol_matrix_size += shell->getNumFunc();
                        start.push_back(shell->getStartIndex());
                        end.push_back(shell->getStartIndex() + shell->getNumFunc());
                        mol_shells.push_back(shell);
                    }
                }
            }
            std::vector<unsigned> startindex;
            std::vector<unsigned> endindex;
            
            if(start.size()>1){
                startindex.push_back(start[0]);

                for (unsigned i = 0; i < start.size() - 1; ++i) {

                    if(end[i]!=start[i+1]){
                        startindex.push_back(start[i+1]);
                        endindex.push_back(end[i]);
                    }
                }
                endindex.push_back(end[end.size()-1]);
            }
            else{
                startindex=start;
                endindex=end;
             }
            unsigned shellstart = 0;
            for (unsigned i = 0; i < startindex.size(); ++i) {
                mol_ranges.push_back(ub::range(startindex[i], endindex[i]));
                
                unsigned size = endindex[i]-startindex[i];
                
                mol_inv_ranges.push_back(ub::range(shellstart, shellstart + size));
                shellstart += size;
            }

            
            return;
    }
    
    
    void GridBox::AddtoBigMatrix_perMolecule(ub::matrix<double>& bigmatrix, const ub::matrix<double>& smallmatrix) {


            for (unsigned i = 0; i < ranges.size(); i++) {
                for (unsigned j = 0; j < mol_ranges.size(); j++) {
                    ub::project(bigmatrix, ranges[i], mol_ranges[j]) += ub::project(smallmatrix, inv_ranges[i], mol_inv_ranges[j]);
                }
            }
            return;
        }
    
    
    ub::matrix<double> GridBox::ReadFromBigMatrix_perMolecule(const ub::matrix<double>& bigmatrix) {

        //the returned matrix is not square. Only Shells from the molecule appear in the columns
        ub::matrix<double> _matrix = ub::zero_matrix<double>(matrix_size, mol_matrix_size);
        for (unsigned i = 0; i < ranges.size(); i++) {
            for (unsigned j = 0; j < mol_ranges.size(); j++) {

                ub::project(_matrix, inv_ranges[i], mol_inv_ranges[j]) = ub::project(bigmatrix, ranges[i], mol_ranges[j]);
            }
        }
        return _matrix;
    }

    
}}