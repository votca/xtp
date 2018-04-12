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
// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>
//for libxc
#include <votca/xtp/votca_config.h>

#include <votca/xtp/numerical_integrations_periodic.h>
#include <boost/math/constants/constants.hpp>
#include <votca/xtp/radial_euler_maclaurin_rule.h>
#include <votca/xtp/sphere_lebedev_rule.h>
#include <votca/xtp/aoshell.h>
#include <votca/tools/constants.h>


#include <votca/xtp/aomatrix.h>
#include <votca/xtp/elements.h>
#include <fstream>
#include <iterator>
#include <string>






namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;
        
        ///Creates an expanded AOBasis and vector of QMAtoms* to include periodic images of atoms up to _nExpantionCells away from the original cell.
        void NumericalIntegrationPeriodic::ExpandBasis(vector<ctp::QMAtom*> _atoms){
            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }
#ifdef DEBUG            
            cout<<"boxLen = "<<boxLen<<" Bohr = "<<boxLen*tools::conv::bohr2ang << "A"<<endl<<flush;
#endif            
            
            
            
            _expanded_basis =  new AOBasis;
            vector< ctp::QMAtom* > ::iterator ait;
            if(_nExpantionCells>0){ //more than one box,  need to expand


#ifdef DEBUG
                //debug: list original shells
    //            for (AOBasis::AOShellIterator _row = _basis->firstShell(); _row != _basis->lastShell(); _row++) {
    //                AOShell* _store=(*_row);
    //                vec orgpos = _store->getPos();
    //                cout<<"\tfrom: "<<_store->getType()<<" at "<<orgpos.getX()*tools::conv::bohr2ang<<" "<<orgpos.getY()*tools::conv::bohr2ang<<" "<<orgpos.getZ()*tools::conv::bohr2ang <<endl<<flush;
    //            }
#endif


                //need to expand
#ifdef DEBUG
                cout<<"\nexpanding atom list into nearby periodic boxes."<<endl<<flush;
                cout<<"_nExpantionCells = "<< _nExpantionCells<<endl<<flush;
#endif

                for(int cx= -_nExpantionCells; cx<=_nExpantionCells; cx++){
                    for(int cy= -_nExpantionCells; cy<=_nExpantionCells; cy++){
                        for(int cz= -_nExpantionCells; cz<=_nExpantionCells; cz++){
                            if(cx==0 && cy==0 && cz==0){ //cell 0 0 0 needs to be at the very end of _expanded_atoms
                                continue;
                            }

                            tools::vec shift = tools::vec(boxLen.getX()*cx, boxLen.getY()*cy, boxLen.getZ()*cz); //Bohr

                            //loop over shells (in Bohr)
                            for (AOBasis::AOShellIterator _row = _basis->firstShell(); _row != _basis->lastShell(); _row++) {
                                AOShell* _store=(*_row);
                                vec imgpos = shift + _store->getPos();
                                AOShell* newShell = _expanded_basis->addShell(_store->getType(), _store->getLmax(), _store->getLmin(), _store->getScale(), _store->getNumFunc(),
                                                          _store->getStartIndex(), _store->getOffset(), imgpos, _store->getName(), _store->getIndex());
                                newShell->copyGaussians(_store);
                                newShell->CalcMinDecay();
                            }                        

                            shift *= tools::conv::bohr2ang; //boxLen is in Bohr, AOBasis positions are in Bohr, QMAtom positions are in Angstroms

                            //loop over atoms (in Angstroms)
                            //for (ait = _atoms.begin() + 1; ait != _atoms.end(); ++ait) { //why no image for first atom in molecule?
                            for (ait = _atoms.begin(); ait != _atoms.end(); ++ait) {
                                vec imgpos = shift + (*ait)->getPos();
                                ctp::QMAtom* imgatom = new ctp::QMAtom((*ait)->type, imgpos.getX(), imgpos.getY(), imgpos.getZ(), (*ait)->charge, (*ait)->from_environment);
#ifdef DEBUG
                                //cout<<"Image atom: "<<imgatom->type<<" at "<<imgatom->x<<" "<<imgatom->y<<" "<<imgatom->z<<" with charge "<<imgatom->charge<<endl<<flush;
                                //cout<<"\tfrom: "<<(*ait)->type<<" at "<<(*ait)->x<<" "<<(*ait)->y<<" "<<(*ait)->z<<" with charge "<<(*ait)->charge<<endl<<flush;
#endif
                                _expanded_atoms.push_back(imgatom);
                                _toclean_atoms.push_back(imgatom);
                            }
                        }
                    }
                }
            
            
            }//more than one box in expansion
            
            //cell 0 0 0 needs to be at the very end of _expanded_atoms
            //atoms
#ifdef DEBUG
            cout<<"adding original atoms to expanded atom list."<<endl<<flush;
            for (ait = _atoms.begin(); ait != _atoms.end(); ++ait) {
                //cout<<"Original atom: "<<(*ait)->type<<" at "<<(*ait)->x<<" "<<(*ait)->y<<" "<<(*ait)->z<<" with charge "<<(*ait)->charge<<"\t ait="<<(*ait)<<endl<<flush;
            }
#endif
            for (ait = _atoms.begin(); ait != _atoms.end(); ++ait) {
                _expanded_atoms.push_back((*ait));
            }
            //shells
            central_cell_shells.clear();
            for (AOBasis::AOShellIterator _row = _basis->firstShell(); _row != _basis->lastShell(); _row++) {
                AOShell* _store=(*_row);
                AOShell* newShell = _expanded_basis->addShell(_store->getType(), _store->getLmax(), _store->getLmin(), _store->getScale(), _store->getNumFunc(),
                                          _store->getStartIndex(), _store->getOffset(), _store->getPos(), _store->getName(), _store->getIndex());
                newShell->copyGaussians(_store);
                newShell->CalcMinDecay();
                central_cell_shells.push_back(newShell);
            }
            

//#ifdef DEBUG
//            cout<<"Expanded Basis:"<<endl<<flush;
//            unsigned int lastInd=-1;
//            for (auto& shell:_expanded_basis->getShells()) {
//                if(shell->getIndex() != lastInd)
//                    cout<<shell->getName()<<" "<<shell->getType()<<" "<<shell->getPos()<< " index:"<<shell->getIndex()<<endl<<flush;
//                lastInd = shell->getIndex();
//            }
//#endif            
            
            
            return;
        }

        
        
        /**\brief Prepares a grid of points around the atoms to evaluate electron density on.
         * 
         * Creates atoms and associated shells to represent periodic images,
         * makes grids around them, filters out irrelevant gridpoints and shells,
         * and groups points into blocks for easy evaluation.
         */
        void NumericalIntegrationPeriodic::GridSetup(string type, BasisSet* bs, vector<ctp::QMAtom*> _atoms,AOBasis* global_basis) {
            
            if(_relevant_atomids.size()==0)
            {
                throw std::runtime_error("_relevant_atomids not set. Run NumericalIntegrationPeriodic::SetRelevantAtomIds() first!"); 
            }
            
            _nExpantionCells=2; //expand basis and atoms to include atoms in this many periodic cells away (2-> 5 cells wide)
            _basis=global_basis;
            ExpandBasis(_atoms);
#ifdef DEBUG            
//            int Ngp_generated=0;
#endif            
            
            
            std::vector< std::vector< GridContainers::integration_grid > > grid;
            const double pi = boost::math::constants::pi<double>();
            // get GridContainer
            GridContainers initialgrids;

            // get radial grid per element
            EulerMaclaurinGrid _radialgrid;
            _radialgrid.getRadialGrid(bs, _atoms, type, initialgrids); // this checks out 1:1 with NWChem results! AWESOME

     
           map<string, GridContainers::radial_grid>::iterator it;

            LebedevGrid _sphericalgrid;
         
            for (it = initialgrids._radial_grids.begin(); it != initialgrids._radial_grids.end(); ++it)
            {
               _sphericalgrid.getSphericalGrid(_atoms, type, initialgrids);
            }

            
            // for the partitioning, we need all inter-center distances later, stored in one-directional list
            int ij = 0;
            vector< ctp::QMAtom* > ::iterator ait;
            vector< ctp::QMAtom* > ::iterator bit;
            int i = 1;
            
#ifdef DEBUG
//            cout<<"\nPre Rij\n"<<flush;
//            for (ait = _atoms.begin(); ait != _atoms.end(); ++ait)
//            {
//                cout<<"ait="<<(*ait)<<endl<<flush;
//            }
//            cout<<"\n\n\n"<<flush;
//            for (bit = _expanded_atoms.begin(); bit != _expanded_atoms.end(); ++bit)
//            {
//                cout<<"bit="<<(*bit)<<endl<<flush;
//            }
#endif

            for (ait = _atoms.begin(); ait != _atoms.end(); ++ait)
            {
                // get center coordinates in Bohr
                vec pos_a = (*ait)->getPos() * tools::conv::ang2bohr;
                int j = 0;
                for (bit = _expanded_atoms.begin(); (*bit) != (*ait); ++bit)
                {
                    ij++;
                    // get center coordinates in Bohr
                    vec pos_b = (*bit)->getPos() * tools::conv::ang2bohr;                   
                    Rij.push_back(1.0 / abs(pos_a-pos_b));
                    j++;
                } // atoms
                Rij.push_back(0.0); // self-distance again
                i++;
            } // atoms
            
            int i_atom = 0;
            _totalgridsize = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait)
            {
                // get center coordinates in Bohr
                std::vector< GridContainers::integration_grid > _atomgrid;
                const vec atomA_pos =(*ait)->getPos() * tools::conv::ang2bohr;
             
                string name = (*ait)->type;
                
                // get radial grid information for this atom type
                GridContainers::radial_grid _radial_grid = initialgrids._radial_grids.at(name);
                
                // get spherical grid information for this atom type
                GridContainers::spherical_grid _spherical_grid = initialgrids._spherical_grids.at(name);

                // maximum order (= number of points) in spherical integration grid
                int maxorder = _sphericalgrid.Type2MaxOrder(name,type);
                int maxindex = _sphericalgrid.getIndexFromOrder(maxorder);

                // for pruning of integration grid, get interval boundaries for this element
                std::vector<double> PruningIntervals = _radialgrid.getPruningIntervals( name );
#ifdef DEBUG
                //cout << " Pruning Intervals: " << PruningIntervals[0] << " " << PruningIntervals[1] << " " << PruningIntervals[2] << " " << PruningIntervals[3] << endl;
#endif
                
                int current_order = 0;
                // get spherical grid
                std::vector<double> _theta;
                std::vector<double> _phi;
                std::vector<double> _weight;
                
                // for each radial value
                for (unsigned _i_rad = 0; _i_rad < _radial_grid.radius.size(); _i_rad++)
                {
                    double r = _radial_grid.radius[_i_rad];
                    int order;
                    // which Lebedev order for this point?
                    if ( maxindex == 1 ) {
                        // smallest possible grid anyway, nothing to do
                        order = maxorder;
                    } else if ( maxindex == 2 ) {
                        // only three intervals
                        if ( r < PruningIntervals[0] ) {
                            order = _sphericalgrid.getOrderFromIndex(1);//1;
                        } else if ( ( r >= PruningIntervals[0] ) && ( r < PruningIntervals[3] )   ){
                            order = _sphericalgrid.getOrderFromIndex(2);
                        } else {
                            order = _sphericalgrid.getOrderFromIndex(1);
                        } // maxorder == 2
                    } else {
                        // five intervals
                        if ( r < PruningIntervals[0] ) {
                            order = _sphericalgrid.getOrderFromIndex(int(2));
                        } else if ( ( r >= PruningIntervals[0]) && ( r < PruningIntervals[1] ) ) {
                            order = _sphericalgrid.getOrderFromIndex(4);
                        } else if ( ( r >= PruningIntervals[1]) && ( r < PruningIntervals[2] ) ) {
                            order = _sphericalgrid.getOrderFromIndex(max(maxindex-1, 4));
                        } else if ( (r >= PruningIntervals[2]) && ( r < PruningIntervals[3] ) ) {
                            order = maxorder;
                        } else {
                            order = _sphericalgrid.getOrderFromIndex(max(maxindex-1,1));
                        }
                    }                        


                    // get new spherical grid, if order changed
                    if ( order != current_order )
                    {
                        _theta.clear();
                        _phi.clear();
                        _weight.clear();
                        
                        _sphericalgrid.getUnitSphereGrid(order,_theta,_phi,_weight);
                        current_order = order;
                    }
                    
                  

                    for (unsigned _i_sph = 0; _i_sph < _phi.size(); _i_sph++)
                    {
                        double p   = _phi[_i_sph] * pi / 180.0; // back to rad
                        double t   = _theta[_i_sph] * pi / 180.0; // back to rad
                        double ws  = _weight[_i_sph];
                        
                        const vec s = vec(sin(p) * cos(t), sin(p) * sin(t),cos(p));

                        //We only need densities in the first periodic cell, but the wrapping isn't done. Ewald summation takes care of it anyway.
                        //Wrapping positions here leads to wrong number of electrons later on.
                        vec ppos = atomA_pos+r*s;
                        GridContainers::integration_grid _gridpoint;
                        _gridpoint.grid_pos = ppos;
                        _gridpoint.grid_weight = _radial_grid.weight[_i_rad] * ws;
                        _atomgrid.push_back(_gridpoint);

                    } // spherical gridpoints
                } // radial gridpoint
                
                
                
                //If wrapping, do it here, so rq has wrapped gridpoint coordinates.
                //But wrapping gridpoints around already wrapped atoms causes the
                //wrapped parts to fully overlap orbitals of image atoms.
                //So DON'T WRAP GRIDPOINTS!
//                tools::vec maxofset(0.0);
//                tools::vec minofset(0.0);
//                for (auto& git: _atomgrid){
////                    cout<<git.grid_pos[0]<<"\t"<<git.grid_pos[1]<<"\t"<<git.grid_pos[2]<<endl;
//                    git.grid_pos = WrapPoint(git.grid_pos, boxLen); //seems to fix the total density
//                }
                
                
//                cout<<flush;
//                cout<<"around atom at: "<<atomA_pos[0]<<"\t"<<atomA_pos[1]<<"\t"<<atomA_pos[2]<<" Bohr"<<endl<<flush;
//                cout<<"max ofset: "<<maxofset[0]<<"\t"<<maxofset[1]<<"\t"<<maxofset[2]<<" Bohr"<<endl<<flush;
//                cout<<"min ofset: "<<minofset[0]<<"\t"<<minofset[1]<<"\t"<<minofset[2]<<" Bohr"<<endl<<flush;
//                
//                cout<<"All extended atoms (also in Bohr):"<<endl;
//                for (auto& b: _expanded_atoms) {
//                    tools::vec posb = b->getPos() * tools::conv::ang2bohr;
//                    cout <<b->type <<"\t" <<posb[0]<<"\t"<<posb[1]<<"\t"<<posb[2]<<" Bohr"<<endl;
//                }
//                cout <<endl<<flush;
////                exit(-1);
                
                
                // get all distances from grid points to centers
                std::vector< std::vector<double> > rq;
                // for each center
                for (bit = _expanded_atoms.begin(); bit < _expanded_atoms.end(); ++bit) {
                    // get center coordinates
                    const vec atom_pos = (*bit)->getPos() * tools::conv::ang2bohr;
                    std::vector<double> temp;
                    
                    // for each gridpoint
                    for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end(); ++git) 
                    {
                        temp.push_back(abs(git->grid_pos-atom_pos));
                    } // gridpoint of _atomgrid
                    rq.push_back(temp); 

                } // centers
                
#ifdef DEBUG
                //cout << " Calculated all gridpoint distances to centers for atom " << i_atom << endl;
                //Ngp_generated+=_atomgrid.size();
#endif
                
//                // find nearest-neighbor of this atom
//                double distNN = 1e10;
//                vector< ctp::QMAtom* > ::iterator NNit;
//               
//                // now check all other centers
//                int i_b =0;
//                for (bit = _expanded_atoms.begin(); bit != _expanded_atoms.end(); ++bit) {
//                    if ((*bit) != (*ait)) {
//                        // get center coordinates
//                        const vec atomB_pos=(*bit)->getPos() * tools::conv::ang2bohr;
//                        double distSQ = (atomA_pos-atomB_pos)*(atomA_pos-atomB_pos);
//
//                        // update NN distance and iterator
//                        if ( distSQ < distNN )
//                        {
//                            distNN = distSQ;
//                            NNit = bit;          
//                        }
//
//                    } // if ( ait != bit) 
//                    i_b++;
//                }// bit centers
                
                
                /*
                 * Periodic weights:
                 * Points are around original atoms in the first periodic cell.
                 * Points aren't wrapped.
                 * Only the atom whose points these are can receive weights.
                 * 
                 * 
                 * Procedure:
                 *  for atom_i:
                 *      for atom_j:
                 *          for atom_j_image:
                 *              p[atom_i] *= p(rq(atom_i,point), rq(atom_j_image, point), Rij(atom_i, atom_j_image))
                 * 
                 * The case of atom_i atom_i_image should also be counted, because
                 * the image contribution goes to the weight of a different point on atom_i
                 * 
                 * Variables:
                 * p[expanded atoms]
                 * rq[expanded atoms][point] not wrapped
                 * Rij[original atoms*expanded atoms/2] not wrapped
                 */
                for ( unsigned i_grid = 0; i_grid < _atomgrid.size() ; i_grid++){
                    // call some shit called grid_ssw0 in NWChem
                    std::vector<double> _p = SSWpartition( i_grid, _atoms.size(), _expanded_atoms.size(), rq);
                    
                    double _p_parent = 0; //summed p of the parent atom where the relative weights have been put onto the original (not expanded) atom
                    for (unsigned i = 0 ; i < _p.size(); i++ )
                    {
                        if( i % _atoms.size() == i_atom){
                            _p_parent += _p[i];
                        }
                    }
                    
                    // check weight sum
                    double wsum = 0.0;
                    //int nadded=(_expanded_atoms.size() - _atoms.size()); //number of centers added by expansion
                    //for (unsigned i = nadded ; i < _p.size(); i++ ) //only sum the weights of the original atoms. Periodic copies don't have grids around them.
                    for (unsigned i =0 ; i < _p.size(); i++ )
                    {
                        wsum += _p[i];
                    }

                    if ( wsum != 0.0 ){
                        
                        // update the weight of this grid point
                        //_p is indexed by expanded_atoms
                        //_atomgrid[i_grid].grid_weight = _atomgrid[i_grid].grid_weight * _p[i_atom + nadded]/wsum;
                        _atomgrid[i_grid].grid_weight = _atomgrid[i_grid].grid_weight * _p_parent/wsum;
                    } else {
                        
                       cerr << "\nSum of partition weights of grid point " << i_grid << " of atom " << i_atom << " is zero! ";
//                       for (unsigned i =0 ; i < _p.size(); i++ )
//                       {
//                           cerr <<"\n _p["<< i <<"] = "<<_p[i]<< "\t\t atom image pos: "<<_expanded_atoms[i]->x<<"\t"<<_expanded_atoms[i]->y<<"\t"<<_expanded_atoms[i]->z;
//                       }
//                       cerr <<"\n number of atoms added by expansion: "<< nadded;
//                       cerr <<"\n atom index = "<< i_atom << "\t so the parent atom index is "<< nadded+i_atom;
//                       cerr <<"\n point position: = "<< _atomgrid[i_grid].grid_pos <<endl;
                       throw std::runtime_error("\nThis should never happen!");                   
                    }
                } // partition weight for each gridpoint
               
                // now remove points from the grid with negligible weights
                for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end();)
                {
                    if (git->grid_weight < 1e-13 )
                    {
                        git = _atomgrid.erase(git);
                    }
                    else
                    {
                        ++git;
                    }
                }                
                _totalgridsize += _atomgrid.size() ;
                grid.push_back(_atomgrid);
                i_atom++;
            } // atoms

            
            SortGridpointsintoBlocks(grid);
#ifdef DEBUG
            //cout<<"grid boxes after SortGridpointsintoBlocks: "<<_grid_boxes.size() <<endl<<flush;
#endif
            FindSignificantShells();
#ifdef DEBUG
            //cout<<"grid boxes after FindSignificantShells: "<<_grid_boxes.size() <<endl<<flush;
#endif
            
            
            //build a vector of ranges to what elements of DMAT & Overlap matrix belong to this molecule
            nFuncInMol=0;
            unsigned aoFuncCounter=0;
            i=0;
            global_mol_aoranges.clear();
            global_mol_inv_aoranges.clear();
            for (AOBasis::AOShellIterator _row = _basis->firstShell(); _row != _basis->lastShell(); _row++) {
                //Check if a  shell belongs to the molecule
                const AOShell* shell=*(_row);
                if(std::find(_relevant_atomids.begin(), _relevant_atomids.end(), shell->getIndex()) != _relevant_atomids.end()){
                    global_mol_aoranges.push_back(ub::range(aoFuncCounter, aoFuncCounter+shell->getNumFunc()));
                    global_mol_inv_aoranges.push_back(ub::range(nFuncInMol, nFuncInMol+shell->getNumFunc()));
                    nFuncInMol += shell->getNumFunc();
                }
                aoFuncCounter += shell->getNumFunc();
            }
            
            
            
#ifdef DEBUG
            //cout<<"Number of Grid points generated = "<<Ngp_generated<<"\tand used = "<< _totalgridsize <<endl<<flush;
#endif
            
            return;
        }
        
        
        std::vector<double> NumericalIntegrationPeriodic::SSWpartition(int igrid, int ncenters, int nexpandedcenters, std::vector< std::vector<double> >& rq)
        {
            const double ass = 0.725;
            // initialize partition vector to 1.0
            std::vector<double> p(nexpandedcenters,1.0);
            
            int nadded=(nexpandedcenters-ncenters); //number of centers added by expansion
#ifdef DEBUG
            //cout<<"nadded= "<< nadded<< endl<< flush;
            //cout<<"igrid= "<< igrid<< endl<< flush;
#endif
            
            const double tol_scr = 1e-10;
            const double leps    = 1e-6; 
            // go through centers
            for ( int ibase = 0; ibase < ncenters; ibase++ ){
                
                int ij = (ibase*(ibase+1)/2 -1)+(ibase*nadded); // indexing magic
                int i = ibase + nadded;
                double rag = rq[i][igrid] ;
                
                // through all other centers (one-directional)
                for (int j = 0; j < i ; j++ ){
                    
                    ij++;
                    if ( ( std::abs(p[i]) > tol_scr  ) || ( std::abs(p[j]) > tol_scr  ) ){
                        
                        double mu = ( rag - rq[j][igrid] )*Rij[ij]; 
                        if ( mu > ass ) {
                            p[i] = 0.0;
                        } else if ( mu < -ass ) {
                            p[j] = 0.0;
                        } else{
                            double sk;
                            if (std::abs(mu) < leps ) {
                                sk = -1.88603178008*mu + 0.5;
                            } else {
                                sk = erf1c(mu); 
                            }
                            if ( mu > 0.0 ) sk = 1.0 - sk;
                            p[j] = p[j] * sk;
                            p[i] = p[i] * (1.0-sk);              
                        }   
                    }  
                }

            }
            
            return p;
        }
        
        /// Filters out shells with insignificant contributions to density
        void NumericalIntegrationPeriodic::FindSignificantShells(){

            for (unsigned i=0;i<_grid_boxes.size();++i){
                GridBox & box=_grid_boxes[i];
                for (AOBasis::AOShellIterator _row = _expanded_basis->firstShell(); _row != _expanded_basis->lastShell(); _row++) {
                      AOShell* _store=(*_row);
                      const double decay=(*_row)->getMinDecay();
                      const tools::vec& shellpos=(*_row)->getPos();
                      
                      for(const auto& point : box.getGridPoints()){
                          tools::vec dist=shellpos-point;
                          double distsq=dist*dist;
#ifdef DEBUG 
                          //cout<<"\tdist: "<<dist.getX()<<"\t"<<dist.getY()<<"\t"<<dist.getZ()<<"\tdistsq= "<<distsq<<"\tdecay= "<<decay<<endl<<flush;
#endif
                          // if contribution is smaller than -ln(1e-10), add atom to list
                        if ( (decay * distsq) < 20.7 ){
#ifdef DEBUG 
                            //cout<< "distsq= "<<distsq<<"\tdecay= "<<decay<<"\t(decay * distsq)="<<(decay * distsq)<<endl<<flush;
#endif
                            box.addShell(_store);
                            break;
                        }
                      }
                }
#ifdef DEBUG 
//                if(i==2365){
//                    cout<<"box "<<i<<"\t has "<<box.Shellsize()<<" shells \t and "<<box.size()<<" points."<<endl;
//                    cout<<"the points are:\n";
//                    for(auto& p:box.getGridPoints()){
//                        cout<<"\t"<<p[0]<<"\t"<<p[1]<<"\t"<<p[2]<<endl;
//                    }
//                    cout<<"the shells are:\n";
//                    for(auto& s:box.getShells()){
//                        cout<<"\t"<<s->getName()<<"\t"<<s->getType()<<"\tmin decay:"<<s->getMinDecay()<<"\tpos:"<<s->getPos()[0]<<"\t"<<s->getPos()[1]<<"\t"<<s->getPos()[2]<<endl;
//                    }
//                }
//                
//                
//                for(const auto& point : box.getGridPoints()){
//                    if(point.getX()>9.75 && point.getX()<9.85 && point.getY()>4.95 && point.getY()<5.05 && point.getZ()>4.06 && point.getZ()<4.16){
//                        cout<<"BOX "<<i<<" has acceptable points"<<endl;
//                        break;
//                    }
//                }
#endif
            }
            std::vector< GridBox > _grid_boxes_copy=_grid_boxes;
            
            std::vector<unsigned> sizes;
            sizes.reserve(_grid_boxes_copy.size());
            for(auto& box: _grid_boxes_copy){
                sizes.push_back(box.size());
            }
           
            
            std::vector<unsigned> indexes=std::vector<unsigned>(sizes.size());
            iota(indexes.begin(), indexes.end(), 0);
            std::sort(indexes.begin(), indexes.end(),[&sizes](unsigned i1, unsigned i2) {return sizes[i1] > sizes[i2];});
            _grid_boxes.resize(0);
            unsigned indexoffirstgridpoint=0;
            for(unsigned& index: indexes){
                if(_grid_boxes_copy[index].Shellsize()>0){
                    GridBox newbox=_grid_boxes_copy[index];
                    newbox.setIndexoffirstgridpoint(indexoffirstgridpoint);
                    indexoffirstgridpoint+=newbox.size();
                    newbox.PrepareForIntegration_perMolecule(_relevant_atomids, central_cell_shells);
                    _grid_boxes.push_back(newbox);
                }   
            }
            
            return;
        }
        
/*
        /// Sorts points into smaller boxes for easy parallelization during density integration.
        void NumericalIntegrationPeriodic::SortGridpointsintoBlocks(std::vector< std::vector< GridContainers::integration_grid > >& grid){
            const double boxsize=1.0;
            
            std::vector< std::vector< std::vector< std::vector< GridContainers::integration_grid* > > > >  boxes;
            
            tools::vec min=vec(std::numeric_limits<double>::max());
            tools::vec max=vec(std::numeric_limits<double>::min());
                   
            min = vec(0.0);
            max = boxLen;

            //Allow for grid positions that extend from atoms in the first periodic box into nearby periodic cells
            //so that there are boxes for gridpoints outside the first periodic cell.
            //These belonging to (wrapped) atoms near the periodic boundaries.
            for ( unsigned i = 0 ; i < grid.size(); i++){
                for ( unsigned j = 0 ; j < grid[i].size(); j++){
                    const tools::vec& pos= grid[i][j].grid_pos;
                    if(pos.getX()>max.getX()){
                        max.x()=pos.getX();
                    }
                    else if(pos.getX()<min.getX()){
                        min.x()=pos.getX();
                    }
                    if(pos.getY()>max.getY()){
                        max.y()=pos.getY();
                    }
                    else if(pos.getY()<min.getY()){
                        min.y()=pos.getY();
                    }
                    if(pos.getZ()>max.getZ()){
                        max.z()=pos.getZ();
                    }
                    else if(pos.getZ()<min.getZ()){
                        min.z()=pos.getZ();
                        }
                    }
                }
            
            
            vec molextension=(max-min);
            vec numberofboxes=molextension/boxsize;
            vec roundednumofbox=vec(std::ceil(numberofboxes.getX()),std::ceil(numberofboxes.getY()),std::ceil(numberofboxes.getZ()));
#ifdef DEBUG
            cout<<"min:\t"<<min[0]<<"\t"<<min[1]<<"\t"<<min[2]<<endl<<flush;
            cout<<"max:\t"<<max[0]<<"\t"<<max[1]<<"\t"<<max[2]<<endl<<flush;
#endif
            
            //creating temparray
            for (unsigned i=0;i<unsigned(roundednumofbox.getX());i++){
                std::vector< std::vector< std::vector< GridContainers::integration_grid* > > > boxes_yz;
                for (unsigned j=0;j<unsigned(roundednumofbox.getY());j++){
                    std::vector< std::vector< GridContainers::integration_grid* > >  boxes_z;
                    for (unsigned k=0;k<unsigned(roundednumofbox.getZ());k++){
                        std::vector< GridContainers::integration_grid* >  box;
                        box.reserve(100);
                        boxes_z.push_back(box);
                    }
                    boxes_yz.push_back(boxes_z);
            }
                boxes.push_back(boxes_yz);
            }
#ifdef DEBUG
            cout<<"boxes size:\t"<<boxes.size()<<"\t"<<boxes[0].size()<<"\t"<<boxes[0][0].size()<<endl<<flush;
#endif            
             for ( auto & atomgrid : grid){
                for ( auto & gridpoint : atomgrid){
                    tools::vec pos = gridpoint.grid_pos - min;
                    tools::vec index=pos/boxsize;
                    int i_x=int(index.getX());
                    int i_y=int(index.getY());
                    int i_z=int(index.getZ());
                    boxes[i_x][i_y][i_z].push_back(&gridpoint);
                }
             }
            
            for ( auto& boxes_xy : boxes){
                for( auto& boxes_z : boxes_xy){
                    for ( auto& box : boxes_z){      
                        if( box.size()<1){
                            continue;
                        }
                        GridBox gridbox;
                        
                        for(const auto&point:box){
                            gridbox.addGridPoint(*point);
                        }
                        _grid_boxes.push_back(gridbox);
                    }
                }
            }
            
            return;
        }
*/                
        
        
        
        double NumericalIntegrationPeriodic::IntegrateDensity(const ub::matrix<double>& _density_matrix, Orbitals& _orbitals){
            
            double N = 0;
            
            unsigned nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif
               
            std::vector<double> N_thread=std::vector<double>(nthreads,0.0);
            
            
            
//            ub::matrix<double> num_AO_bigmat = ub::zero_matrix<double>(_density_matrix.size1(), _density_matrix.size2());
               
            
            
            //need to project every relevant shell onto every relevant shell in this molecule
            #pragma omp parallel for
            for (unsigned thread=0;thread<nthreads;++thread){
                for (unsigned i = thread; i < _grid_boxes.size(); i+=nthreads) {

                    double N_box=0.0;
                    GridBox& box = _grid_boxes[i];

                    //gcc + MKL does some weird memory optimizations with parts of DMAT_here
                    //resulting in garbage matrix multiplications. icc + MKL works, though.
                    const ub::matrix<double>  DMAT_here=box.ReadFromBigMatrix_perMolecule(_density_matrix);

                    const std::vector<tools::vec>& points=box.getGridPoints();
                    const std::vector<double>& weights=box.getGridWeights();

                    ub::range one=ub::range(0,1);

                    ub::matrix<double> ao=ub::matrix<double>(1,box.Matrixsize());
                    ub::matrix<double> ao_mol=ub::matrix<double>(1,box.Matrixsize_perMolecule());

                    box.prepareDensity();
#ifdef DEBUG
                    //cout<<endl<<"iterate over gridpoints"<<endl<<flush;
#endif
                    //iterate over gridpoints
                    for(unsigned p=0;p<box.size();p++){
                        //for row vector: use all significant shells
                        ao=ub::zero_matrix<double>(1,box.Matrixsize());
                        const std::vector<ub::range>& aoranges=box.getAOranges();
                        const std::vector<const AOShell* > shells=box.getShells();
                        for(unsigned j=0;j<box.Shellsize();++j){
                            const AOShell* shell=shells[j];
                            ub::matrix_range< ub::matrix<double> > aoshell=ub::project(ao,one,aoranges[j]);

                            shell->EvalAOspace(aoshell,points[p]);
                        }
                        
                        
                        ub::matrix<double> _temp=ub::prod( ao, DMAT_here);
                        
                        //For the column vector: only significant shells in the relevant molecule in the first periodic cell
                        ao_mol=ub::zero_matrix<double>(1,box.Matrixsize_perMolecule());
                        const std::vector<ub::range>& aoranges_mol=box.getAOranges_perMolecule();
                        const std::vector<const AOShell* > shells_mol=box.getShells_perMolecule(); //should only include shells in the first periodic box
                        for(unsigned j=0;j<shells_mol.size();++j){
                            const AOShell* shell=shells_mol[j];
                            ub::matrix_range< ub::matrix<double> > aoshell=ub::project(ao_mol,one,aoranges_mol[j]);

                            shell->EvalAOspace(aoshell,points[p]);
                        }


                        ub::matrix<double> tr = ub::trans( ao_mol);                        
                        ub::matrix<double> pr = ub::prod(_temp, tr );
                        double rho = pr(0,0);
                        box.addDensity(rho);
                        N_box+=rho*weights[p];
                        
//                        ub::matrix<double> local_AO = ub::prod(tr, ao ) * weights[p];
//                        #pragma omp critical
//                        {
//                            box.AddtoBigMatrix(num_AO_bigmat, local_AO);
//                        }
                        
                        
//                        tools::vec targetP1(2.8288,-2.36288,0.0259489);
//                        tools::vec targetP2(1.96248,3.83502,4.39304);
//                        if(abs(points[p]-targetP1)<0.2 || abs(points[p]-targetP2)<0.1){
//                            cout<<points[p][0]<<" "<<points[p][1]<<" "<<points[p][2]<<" rho="<<rho<<" w="<<weights[p]<<endl<<flush;
//                        }
                        
                    }
                    N_thread[thread]+=N_box;

                }
            }   
            for(unsigned int i=0;i<nthreads;++i){
                N+=N_thread[i];
            }
            
               
               
            //find the number of electrons directly from Density and Overlap matrices
            AOOverlapPeriodic overlap;
            overlap.setBox(boxLen); //in Bohr
            cout<< "Setting up periodic overlap with box: "<<boxLen[0]<<"\t"<<boxLen[1]<<"\t"<<boxLen[2]<< " in Bohr."<<endl<<flush;
            //this is the global, non expanded basis
            overlap.Fill(*_basis);   //AOOverlapPeriodic will build an overlap matrix taking periodicity into account here
            //_density_matrix is the _global_dmat
            ub::matrix<double>& AO=overlap.Matrix();
            double N_comp=0.0;
            
            //make the submatrices
            //Only Shells from the molecule appear in the columns
            ub::matrix<double> DMAT_submat = ub::zero_matrix<double>(AO.size1(), nFuncInMol);
            ub::matrix<double> AO_submat = ub::zero_matrix<double>(AO.size1(), nFuncInMol);
            ub::range everything = ub::range(0, AO.size1());
            for (unsigned i = 0; i < global_mol_aoranges.size(); i++) {
                ub::project(DMAT_submat, everything, global_mol_inv_aoranges[i]) = ub::project(_density_matrix, everything, global_mol_aoranges[i]);
                ub::project(AO_submat  , everything, global_mol_inv_aoranges[i]) = ub::project(AO             , everything, global_mol_aoranges[i]);
            }
            
            //sum up N_comp
            for (unsigned i = 0; i < AO_submat.size1(); i++) {
                for (unsigned j = 0; j < AO_submat.size2(); j++) {
                    N_comp += DMAT_submat(i,j)*AO_submat(i,j);
                }
            }
            
            
            //make the submatrices
            //Only Shells from the molecule appear in the columns
            AOOverlap overlap_np;
            overlap_np.Fill(*_basis);
            AO=overlap_np.Matrix();
            double N_comp_non_periodic=0.0;
            AO_submat = ub::zero_matrix<double>(AO.size1(), nFuncInMol);
            for (unsigned i = 0; i < global_mol_aoranges.size(); i++) {
                ub::project(AO_submat  , everything, global_mol_inv_aoranges[i]) = ub::project(AO             , everything, global_mol_aoranges[i]);
            }
            
            //sum up N_comp
            for (unsigned i = 0; i < AO_submat.size1(); i++) {
                for (unsigned j = 0; j < AO_submat.size2(); j++) {
                    N_comp_non_periodic += DMAT_submat(i,j)*AO_submat(i,j);
                }
            }
            
            
            
            
            
//#ifdef DEBUG
//            double N_direct=0.0;
//            for (unsigned i = 0; i < _density_matrix.size1(); i++) {
//                for (unsigned j = 0; j < _density_matrix.size2(); j++) {
//                    N_direct += _density_matrix(i,j)*AO(i,j);
//                }
//            }
//            cout << "N="<<N<<"\tN_comp="<<N_comp<<"\tN_direct="<<N_direct<<endl<<endl<<flush;
//            
//            ub::matrix<double> MO = _orbitals.MOCoefficients();
//            cout << "MO size:"<< MO.size1() <<" x "<< MO.size2() <<endl<<flush;
//            
//            cout << "MO linear independence:"<<endl;
//            ub::range all_basis_funcs = ub::range(0, MO.size2());
//            for (unsigned i = 0; i < MO.size1(); i++) {
//                ub::range ri = ub::range(i, i+1);
//                ub::matrix<double> Ci = ub::project(MO, ri, all_basis_funcs);
//                for (unsigned j = 0; j < MO.size1(); j++) {
//                    ub::range rj = ub::range(j, j+1);
//                    ub::matrix<double> Cj = ub::trans(ub::project(MO, rj, all_basis_funcs));
//                    ub::matrix<double> SCj = ub::prod(AO,Cj);
//                    cout << ub::prod(Ci,SCj)(0,0) <<"\t";
//                }
//                cout<<endl<<endl;
//            }
//            cout<<endl<<flush;
//            
//            cout<<"MO coefficients:"<<endl;
//            for (unsigned i = 0; i < MO.size1(); i++) {
//                for (unsigned j = 0; j < MO.size1(); j++) {
//                    cout<<MO(i,j)<<"\t";
//                }
//                cout<<endl;
//            }
//            cout<<endl<<flush;
//            
//            //exit(0);
//#endif
            
            cout << "N_comp from AO & DMAT is: " << N_comp << "\t and N from numerical density integration is: " << N << endl << flush;
            cout << "N_comp_non_periodic from AO & DMAT is: " << N_comp_non_periodic << endl << flush;
            //check if the numbers of electrons are the same
            if(std::abs(N-N_comp)>0.005){
                cout << "N_comp from AO & DMAT is: " << N_comp << "\t and N from numerical density integration is: " << N << endl << flush;
#ifdef DEBUG
//                cout << "N_direct is: "<<N_direct << endl << flush;
#endif
                cout <<"=======================" << endl << flush; 
                cout <<"WARNING: Calculated Densities at Numerical Grid with boxes, Number of electrons "<< N <<" is far away from the the real value "<< N_comp<<", you should increase the accuracy of the integration grid." << endl << flush; 
                cout <<"=======================" << endl << flush; 
                cout <<"#of Boxes=" << _grid_boxes.size() << endl << flush; 
                exit(-1);
            }
               
            density_set=true;
            return N;
        }
        

        
        
        ///Computes the dipole moment of a periodic density distribution
        double NumericalIntegrationPeriodic::CalcDipole_w_PBC(vec rvector) {
            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }
            //as a check, compute the dipole moment of the density distribution.
            //center it on the first atom, for example.
            

            double z[3] = {0};
            vec dip(z);
            if (density_set) {
                const std::vector<double>& _densities = _periodicGridBox.getGridDensities();
                const std::vector<double>& _weights = _periodicGridBox.getGridWeights();
                const std::vector<tools::vec>& _positions = _periodicGridBox.getGridPoints();
                
                #pragma omp parallel
                {
                    vec thread_local_dip(0.0);  //avoids false sharing
                    double q;
                    #pragma omp for
                    for (unsigned i = 0; i < _densities.size(); i++) {
                        q = -_densities[i] * _weights[i]; //density is neg of charge
                        vec dif = WrapDisplacement(_positions[i], rvector, boxLen);
                        thread_local_dip += dif * q;
                    }//i
                    
                    #pragma omp critical
                    dip += thread_local_dip;
                }
            } else {
                throw std::runtime_error("Density not calculated");
            }
            return (abs(dip)); //in Bohr*elementary charge
        }

        
        
        // Prepares inverse space representation of density and the evaluation grid for Ewald Summation. Algorithm adapted from gromacs 4.6.
        void NumericalIntegrationPeriodic::PrepKspaceDensity_gromacs_like(double ext_alpha, std::vector< ctp::QMAtom* > & _local_atomlist, bool ECP, Grid &eval_grid, int nK) {
            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }
            
            
            if (!density_set) {
                throw std::runtime_error("Density not calculated");
            }

            //check if _periodicGridBox is already set up, if not fill it
            if(_periodicGridBox.size()==0){
#ifdef DEBUG
                cout<<"Filling up _periodicGridBox"<<endl<<flush;
#endif
                //fill the _periodicGridBox from all the smaller grid boxes
                //The coordinates of all the gridpoints in the smaller boxes are sometimes outside the main periodic cell.
                //So, need to wrap the positions of the densities into the first periodic cell.
                for (unsigned i = 0; i < _grid_boxes.size(); i++) {
                    _periodicGridBox.appendBoxData(_grid_boxes[i]);
                }
                
                //already wrapped
                //wrap
                //_periodicGridBox.wrapPositions(boxLen); // wraping positions here works better than wraping when they are first created.
#ifdef DEBUG                
                cout<<"_periodicGridBox.size() due to electrons: "<<_periodicGridBox.size()<<endl<<flush;
#endif
                //add nuclear charges to _periodicGridBox to avoid having to do periodic calculations on them separately;
                Elements _elements;
                for (std::vector< ctp::QMAtom* >::iterator it = _local_atomlist.begin(); it != _local_atomlist.end(); ++it) {
                    GridContainers::integration_grid el;
                    ctp::QMAtom* atom = *it;
                    if (ECP) {
                        _periodicGridBox.addDensity(-_elements.getNucCrgECP(atom->type));
                    } else {
                        _periodicGridBox.addDensity(-_elements.getNucCrg(atom->type));
                    }
                    el.grid_weight = 1.0;
                    el.grid_pos = vec(atom->x, atom->y, atom->z) * tools::conv::ang2bohr;
                    _periodicGridBox.addGridPoint(el);
                }
            }
            
            const std::vector<double>& _densities = _periodicGridBox.getGridDensities();
            const std::vector<tools::vec>& _positions = _periodicGridBox.getGridPoints();
            const std::vector<double>& _weights = _periodicGridBox.getGridWeights();
            
            double totQ=0;
            double totQweighted=0;
            for(unsigned int i=0; i<_densities.size(); i++){
                totQ+=_densities[i];
                totQweighted+=_densities[i]*_weights[i];
            }
#ifdef DEBUG
            cout<<"totQ= "<<totQ<<"\ttotQweighted= "<<totQweighted<<endl<<flush;
#endif            

            alpha = ext_alpha;

            //this is going to be slow, as points we have density for are not on a periodic grid,
            //so will use simple Ewald summation, not any of the FFT methods.

            //find number of k-vectors in each lattice direction
            for (int i = 0; i < 3; i++)
            {
                numK[i] = nK;
            }
            
#ifdef DEBUG
            cout << "numK={" << numK[0] << ", " << numK[1] << ", " << numK[2] << "}" << endl<<flush;
#endif


            //allocate space for eikr
            eikr.resize(3); //3 dim
            for (unsigned m = 0; m < 3; m++) { //dimension
                eikr[m].resize(numK[m]);
                
                //#pragma omp parallel for
                for (int k = 0; k < numK[m]; k++) { //k vectors
                    eikr[m][k].resize(_densities.size());    //density points
                }

                //also compute the inverse box dimensions
                lll[m] = 2.0 * boost::math::constants::pi<double>() / boxLen[m]; //units are 1/Bohr
            }
            //eikr indexing is [DIM][k][i][j]
#ifdef DEBUG
            cout << "lll={" << lll[0] << ", " << lll[1] << ", " << lll[2] << "}" << endl<<flush;
#endif
            
            //based on gromacs 4.6 tabulate_eir()
            //#pragma omp parallel for
            for (unsigned i = 0; i < _densities.size(); i++) { //density points
                vec r = _positions[i]; //grid is in Bohr, need it in reduced units (*2pi/box))
                for (unsigned m = 0; (m < 3); m++) { //dimensions
                    eikr[m][0][i] = std::complex<double>(1, 0);
                    eikr[m][1][i] = std::complex<double>(cos(r[m] * lll[m]), sin(r[m] * lll[m])); //this is where reduced units are applied
                    for (int k = 2; k < numK[m]; k++) { //k vectors
                        eikr[m][k][i] = eikr[m][k - 1][i] * eikr[m][1][i];
                    }//k
                }//m
            }//i


            //now do the same for the eval_grid (eikR), the points were potential is to be evaluated
            //allocate space for eikR
            eikR.resize(3); //3 dim
            for (unsigned m = 0; m < 3; m++) { //dimensions
                eikR[m].resize(numK[m]);
                for (int k = 0; k < numK[m]; k++) {
                    eikR[m][k].resize(eval_grid.getGrid().size());
                }
            }
            //eikR indexing is [DIM][k][i]

            //based on gromacs 4.6 tabulate_eir()
            //#pragma omp parallel for
            for (unsigned i = 0; i < eval_grid.getGrid().size(); i++) {
                vec R = eval_grid.getGrid()[i] * tools::conv::nm2bohr; //this is in Bohr, eval_grid is in nm
                for (unsigned m = 0; (m < 3); m++) {
                    eikR[m][0][i] = std::complex<double>(1, 0);
                    eikR[m][1][i] = std::complex<double>(cos(R[m] * lll[m]), sin(R[m] * lll[m]));
                    for (int k = 2; k < numK[m]; k++) {
                        eikR[m][k][i] = eikR[m][k - 1][i] * eikR[m][1][i];
                    }//k
                }//m
            }//i
        }

        /// Integrates the potential from a periodic charge distribution using Ewald Summation. Algorithm adapted from gromacs 4.6.
        void NumericalIntegrationPeriodic::IntegratePotential_w_PBC_gromacs_like(Grid &eval_grid, ub::vector<double>& _ESPatGrid) {
            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }
            
            double cutoff = min(min(boxLen[0], boxLen[1]), boxLen[2]) / 2.0; //Bohr
            double vol = boxLen[0] * boxLen[1] * boxLen[2]; //Bohr^3
            //eval_grid is in nm


            ub::vector<double> k_pot;
            k_pot.resize(_ESPatGrid.size());
            for (unsigned p = 0; p < k_pot.size(); p++) {
                k_pot[p] = 0;
            }


            if (density_set) {
                const std::vector<double>& _densities = _periodicGridBox.getGridDensities();
                const std::vector<double>& _weights = _periodicGridBox.getGridWeights();
                const std::vector<tools::vec>& _positions = _periodicGridBox.getGridPoints();
                
                //real space part
                #pragma omp parallel for
                for (unsigned p = 0; p < _ESPatGrid.size(); p++) {
                    _ESPatGrid[p] = 0;
                    vec rvector = eval_grid.getGrid()[p] * tools::conv::nm2bohr; //Bohr
                    for (unsigned i = 0; i < _densities.size(); i++) {

                            //charge at this point
                            double q = -_weights[i] * _densities[i]; //density is neg of charge

                            //r-space sum
                            vec dif = WrapDisplacement(_positions[i], rvector, boxLen);
                            double dist = abs(dif); //in bohr
                            
                            double potR = q * ((erfc(alpha * dist) / dist) - tools::conv::Pi / (vol * vol * alpha * alpha)); //erfc and shift average pot to 0

                            if (dist < 1.0e-12) { //point is in the same spot as we are evaluating potential
                                _ESPatGrid[p] -= 2.0 * (alpha / sqrt(tools::conv::Pi)) * q; //self correction
                            } else if (dist <= cutoff) {
                                //double potR=q/dist;
                                _ESPatGrid[p] += potR;
                            }
                    }//i
                }//p



                int ix, iy, iz;
                double tmp, cs, ss, ak, mx, my, mz, m2;
                double factor = -1.0 / (4 * alpha * alpha);
                double scaleRecip = 4.0 * boost::math::constants::pi<double>() / (vol); //final units of potential are Hartree/e (atomic units)
                std::vector<std::complex<double> > tab_xy;
                std::vector<std::complex<double> > tab_qxyz; //charge distrib.
                std::vector<std::complex<double> > tab_R_xyz; //where to evaluate
                std::vector<std::complex<double> > tab_R_xy;
                //allocate space
                tab_xy.resize(_densities.size());
                tab_qxyz.resize(_densities.size());
                tab_R_xy.resize(eval_grid.getsize());
                tab_R_xyz.resize(eval_grid.getsize());

                //inverse space part
                unsigned int lowiy = 0;
                unsigned int lowiz = 1;
                for (ix = 0; ix < numK[0]; ix++) {
                    mx = ix * lll[0];
                    for (iy = lowiy; iy < numK[1]; iy++) {
                        my = iy * lll[1];
                        if (iy >= 0) {
                            for (unsigned i = 0; i < _densities.size(); i++) //charge distrib.
                            {
                                tab_xy[i] = eikr[0][ix][i] * eikr[1][iy][i]; //eikr indexing is [DIM][k][i]
                            }//i
                            for (unsigned n = 0; n < eval_grid.getsize(); n++) //where to evaluate
                            {
                                tab_R_xy[n] = eikR[0][ix][n] * eikR[1][iy][n];
                            }//n
                        } else //negative k component
                        {
                            for (unsigned i = 0; i < _densities.size(); i++) //charge distrib.
                            {
                                tab_xy[i] = eikr[0][ix][i] * std::conj(eikr[1][-iy][i]); //eikr indexing is [DIM][k][i][j]
                            }//i
                            for (unsigned n = 0; n < eval_grid.getsize(); n++) //where to evaluate
                            {
                                tab_R_xy[n] = eikR[0][ix][n] * std::conj(eikR[1][-iy][n]);
                            }//n
                        }
                        
                        //#pragma omp parallel for
                        for (iz = lowiz; iz < numK[2]; iz++) {
                            mz = iz * lll[2];
                            m2 = mx * mx + my * my + mz * mz;
                            ak = exp(m2 * factor) / m2;
                            if (iz >= 0) {
                                for (unsigned i = 0; i < _densities.size(); i++) //charge distrib.
                                {
                                    tab_qxyz[i] = -_weights[i] * _densities[i] * tab_xy[i] * eikr[2][iz][i]; //eikr indexing is [DIM][k][i]
                                }//i
                                for (unsigned n = 0; n < eval_grid.getsize(); n++) //where to evaluate
                                {
                                    tab_R_xyz[n] = tab_R_xy[n] * eikR[2][iz][n];
                                }//n
                            } else {
                                for (unsigned i = 0; i < _densities.size(); i++) //charge distrib.
                                {
                                    tab_qxyz[i] = -_weights[i] * _densities[i] * tab_xy[i] * std::conj(eikr[2][-iz][i]); //eikr indexing is [DIM][k][i][j]
                                }//i
                                for (unsigned n = 0; n < eval_grid.getsize(); n++) //where to evaluate
                                {
                                    tab_R_xyz[n] = tab_R_xy[n] * std::conj(eikR[2][-iz][n]);
                                }//n
                            }

                            cs = ss = 0;
                            for (unsigned i = 0; i < _densities.size(); i++) {
                                cs += tab_qxyz[i].real();
                                ss += tab_qxyz[i].imag();
                            }
                            for (unsigned int n = 0; n < eval_grid.getsize(); n++) {
                                tmp = ak * (cs * tab_R_xyz[n].real() + ss * tab_R_xyz[n].imag());
                                _ESPatGrid(n) += tmp * 2 * scaleRecip;
                                k_pot(n) += tmp * 2 * scaleRecip;
                            }
                        }//iz
                        lowiz = 1 - numK[2];                       
                    }//iy
                    lowiy = 1 - numK[1];
                }//ix

            } else {
                throw std::runtime_error("Density not calculated");
            }
        }
  
        
        /**
         * \brief Wraps a vector to periodic boundary conditions. This eventually needs to be moved to tools::vec.
         * 
         * @param vector r to wrap
         * @param vector box specifying the periodic box
         * @return the wrapped vector
         */
        tools::vec WrapPoint(const tools::vec r, const tools::vec box) {
            tools::vec ret;
            for(int k=0; k<3; k++)
            {
                ret[k] = fmod(r[k], box[k]);
                if(ret[k]<0) ret[k] = box[k]+ret[k];
            }
            return(ret);
        }
        
        /**
         * \brief Finds the shortest displacement between images of two vectors due to periodic boundary conditions.
         * This eventually needs to be moved to tools::vec.
         * 
         * @param vector a
         * @param vector b
         * @param vector box specifying the periodic box
         * @return the shortest displacement
         */
        tools::vec WrapDisplacement(const tools::vec a, const tools::vec b , const tools::vec box) {
            tools::vec ret = a-b;
            for(int k=0; k<3; k++)
            {
                ret[k] = fmod(ret[k], box[k]); //if ret <0, then fmod<0 
                if(ret[k] < 0){
                    ret[k] += box[k];
                }
                if(ret[k] > box[k]*0.5){
                    ret[k] = box[k]-ret[k];
                }
                
            }
            return(ret);
        }
        

    }//xtp
}//votca
