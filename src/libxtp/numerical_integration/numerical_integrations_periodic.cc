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
        
        
        void NumericalIntegrationPeriodic::ExpandBasis(vector<ctp::QMAtom*> _atoms){
            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }
            
            cout<<"boxLen = "<<boxLen<<" Bohr = "<<boxLen*tools::conv::bohr2ang << "A"<<endl<<flush;
            
            if(_nExpantionCells==0){ //just one box, no need to expand
                _expanded_basis=_basis;
                _expanded_atoms=_atoms;
                //no _toclean_atoms
                return;
            }
            
//            //debug: list original shells
//            for (AOBasis::AOShellIterator _row = _basis->firstShell(); _row != _basis->lastShell(); _row++) {
//                AOShell* _store=(*_row);
//                vec orgpos = _store->getPos();
//                cout<<"\tfrom: "<<_store->getType()<<" at "<<orgpos.getX()*tools::conv::bohr2ang<<" "<<orgpos.getY()*tools::conv::bohr2ang<<" "<<orgpos.getZ()*tools::conv::bohr2ang <<endl<<flush;
//            }
            
            
            
            //need to expand
            cout<<"\nexpanding atom list into nearby periodic boxes."<<endl<<flush;
//            cout<<"_nExpantionCells = "<< _nExpantionCells<<endl<<flush;
            _expanded_basis =  new AOBasis;
            vector< ctp::QMAtom* > ::iterator ait;
            for(int cx= -_nExpantionCells; cx<=_nExpantionCells; cx++){
                for(int cy= -_nExpantionCells; cy<=_nExpantionCells; cy++){
                    for(int cz= -_nExpantionCells; cz<=_nExpantionCells; cz++){
//                        cout<<"box "<<cx<<""<<cy<<""<<cz<<endl<<flush;
                        if(cx==0 && cy==0 && cz==0){ //cell 0 0 0 needs to be at the very end of _expanded_atoms
                            continue;
                        }
                        
                        tools::vec shift = tools::vec(boxLen.getX()*cx, boxLen.getY()*cy, boxLen.getZ()*cz); //Bohr
                        
                        //loop over shells (in Bohr)
                        for (AOBasis::AOShellIterator _row = _basis->firstShell(); _row != _basis->lastShell(); _row++) {
                            AOShell* _store=(*_row);
                            vec imgpos = shift + _store->getPos();
//                            vec orgpos = _store->getPos();
//                            cout<<"Image shell: "<<_store->getType()<<" at "<<imgpos.getX()<<" "<<imgpos.getY()<<" "<<imgpos.getZ()<<" with start index "<<_store->getStartIndex()<<endl<<flush;
//                            cout<<"\tfrom: "<<_store->getType()<<" at "<<orgpos.getX()<<" "<<orgpos.getY()<<" "<<orgpos.getZ() <<endl<<flush;
                            AOShell* newShell = _expanded_basis->addShell(_store->getType(), _store->getLmax(), _store->getLmin(), _store->getScale(), _store->getNumFunc(),
                                                      _store->getStartIndex(), _store->getOffset(), imgpos, _store->getName(), _store->getIndex());
                            newShell->copyGaussians(_store);
                            newShell->CalcMinDecay();
                        }                        
                        
                        shift *= tools::conv::bohr2ang; //boxLen is in Bohr, AOBasis positions are in Bohr, QMAtom positions are in Angstroms

                        //loop over atoms (in Angstroms)
                        for (ait = _atoms.begin() + 1; ait != _atoms.end(); ++ait) {
                            vec imgpos = shift + (*ait)->getPos();
                            ctp::QMAtom* imgatom = new ctp::QMAtom((*ait)->type, imgpos.getX(), imgpos.getY(), imgpos.getZ(), (*ait)->charge, (*ait)->from_environment);
//                            cout<<"Image atom: "<<imgatom->type<<" at "<<imgatom->x<<" "<<imgatom->y<<" "<<imgatom->z<<" with charge "<<imgatom->charge<<endl<<flush;
//                            cout<<"\tfrom: "<<(*ait)->type<<" at "<<(*ait)->x<<" "<<(*ait)->y<<" "<<(*ait)->z<<" with charge "<<(*ait)->charge<<endl<<flush;
                            _expanded_atoms.push_back(imgatom);
                            _toclean_atoms.push_back(imgatom);
                        }
                    }
                }
            }
            
            //cell 0 0 0 needs to be at the very end of _expanded_atoms
            //atoms
//            cout<<"adding original atoms to expanded atom list."<<endl<<flush;
            for (ait = _atoms.begin(); ait != _atoms.end(); ++ait) {
                _expanded_atoms.push_back((*ait));
                cout<<"Original atom: "<<(*ait)->type<<" at "<<(*ait)->x<<" "<<(*ait)->y<<" "<<(*ait)->z<<" with charge "<<(*ait)->charge<<"\t ait="<<(*ait)<<endl<<flush;
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
            
//            exit(0);
            return;
        }

        
        
        
        void NumericalIntegrationPeriodic::GridSetup(string type, BasisSet* bs, vector<ctp::QMAtom*> _atoms,AOBasis* global_basis) {
            
            if(_relevant_atomids.size()==0){
                throw std::runtime_error("_relevant_atomids not set. Run NumericalIntegrationPeriodic::SetRelevantAtomIds() first!"); 
            }
            
            _nExpantionCells=2; //expand basis and atoms to include atoms in this many periodic cells away (2-> 5 cells wide)
            _basis=global_basis;
            ExpandBasis(_atoms);
            
            
            
            std::vector< std::vector< GridContainers::integration_grid > > grid;
            const double pi = boost::math::constants::pi<double>();
            // get GridContainer
            GridContainers initialgrids;

            // get radial grid per element
            EulerMaclaurinGrid _radialgrid;
            _radialgrid.getRadialGrid(bs, _atoms, type, initialgrids); // this checks out 1:1 with NWChem results! AWESOME

     
           map<string, GridContainers::radial_grid>::iterator it;

            LebedevGrid _sphericalgrid;
         
            for (it = initialgrids._radial_grids.begin(); it != initialgrids._radial_grids.end(); ++it) {
               _sphericalgrid.getSphericalGrid(_atoms, type, initialgrids);
       
            }

            
            // for the partitioning, we need all inter-center distances later, stored in one-directional list
            int ij = 0;
            //Rij.push_back(0.0); // 1st center "self-distance"
//            cout<<"\nPre Rij\n"<<flush;
            vector< ctp::QMAtom* > ::iterator ait;
            vector< ctp::QMAtom* > ::iterator bit;
            int i = 1;
//            for (ait = _atoms.begin(); ait != _atoms.end(); ++ait) {
//            cout<<"ait="<<(*ait)<<endl<<flush;
//            }
//            cout<<"\n\n\n"<<flush;
//            for (bit = _expanded_atoms.begin(); bit != _expanded_atoms.end(); ++bit) {
//            cout<<"bit="<<(*bit)<<endl<<flush;
//            }

            for (ait = _atoms.begin(); ait != _atoms.end(); ++ait) {
                // get center coordinates in Bohr
                vec pos_a = (*ait)->getPos() * tools::conv::ang2bohr;
//                cout<<"\natom; "<<"#_expanded_atoms: "<<_expanded_atoms.size()<<endl<<flush;
                int j = 0;
                for (bit = _expanded_atoms.begin(); (*bit) != (*ait); ++bit) {
//                    cout<<"expanded_atom "<<j<<"\t "<<(*bit)->type<<"\t"<<(*ait)<<"\t"<<(*bit)<<endl<<flush;
                    ij++;
                    // get center coordinates in Bohr
                    vec pos_b = (*bit)->getPos() * tools::conv::ang2bohr;
                    
                   
                    Rij.push_back(1.0 / abs(pos_a-pos_b));
                                        
                    j++;
                } // atoms
                Rij.push_back(0.0); // self-distance again
//                cout<<"\tSelf dist (Rij) at slot"<<Rij.size()<<endl<<flush;
                i++;
//                cout<<"Rij\n"<<flush;
            } // atoms
//            cout<<"\nGridSetup Debug1: Post Rij\n"<<flush;
            
            int i_atom = 0;
            _totalgridsize = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
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
              //  cout << " Pruning Intervals: " << PruningIntervals[0] << " " << PruningIntervals[1] << " " << PruningIntervals[2] << " " << PruningIntervals[3] << endl;
                
                int current_order = 0;
                // get spherical grid
                std::vector<double> _theta;
                std::vector<double> _phi;
                std::vector<double> _weight;
                
                // for each radial value
                for (unsigned _i_rad = 0; _i_rad < _radial_grid.radius.size(); _i_rad++) {
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
                    if ( order != current_order ){
                        _theta.clear();
                        _phi.clear();
                        _weight.clear();
                        
                        _sphericalgrid.getUnitSphereGrid(order,_theta,_phi,_weight);
                        current_order = order;
                    }
                    
                  

                    for (unsigned _i_sph = 0; _i_sph < _phi.size(); _i_sph++) {

                        double p   = _phi[_i_sph] * pi / 180.0; // back to rad
                        double t   = _theta[_i_sph] * pi / 180.0; // back to rad
                        double ws  = _weight[_i_sph];

                        const vec s = vec(sin(p) * cos(t), sin(p) * sin(t),cos(p));
                     

                        //we only need densities in the first periodic cell
                        vec ppos = atomA_pos+r*s;
//                        if(ppos.getX()>=0 && ppos.getY()>=0 && ppos.getZ()>=0 &&
//                           ppos.getX()<boxLen.getX() && ppos.getY()<boxLen.getY() && ppos.getZ()<boxLen.getZ()){

                            GridContainers::integration_grid _gridpoint;
                            //wrap position into the periodic box. Pos is in Bohr.
                            //Need to be wrapped into first periodic cell for potential evaluation from these points to make sense later.
                            //_gridpoint.grid_pos = WrapPoint(ppos, boxLen);
                            _gridpoint.grid_pos = ppos;
                            _gridpoint.grid_weight = _radial_grid.weight[_i_rad] * ws;
                            _atomgrid.push_back(_gridpoint);
//                        }


                    } // spherical gridpoints
                } // radial gridpoint
                
                // get all distances from grid points to centers
                std::vector< std::vector<double> > rq;
                // for each center
                for (bit = _expanded_atoms.begin(); bit < _expanded_atoms.end(); ++bit) {
                    // get center coordinates
                   const vec atom_pos = (*bit)->getPos() * tools::conv::ang2bohr;


                    std::vector<double> temp;
                    // for each gridpoint
                    for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end(); ++git) {

                        temp.push_back(abs(git->grid_pos-atom_pos));

                    } // gridpoint of _atomgrid
                    rq.push_back(temp); 

                } // centers
                // cout << " Calculated all gridpoint distances to centers for " << i_atom << endl;
                
                // find nearest-neighbor of this atom
                double distNN = 1e10;

                vector< ctp::QMAtom* > ::iterator NNit;
                //int i_NN;
               
                // now check all other centers
                int i_b =0;
                for (bit = _expanded_atoms.begin(); bit != _expanded_atoms.end(); ++bit) {

                    if ((*bit) != (*ait)) {
                        // get center coordinates
                       
                        const vec atomB_pos=(*bit)->getPos() * tools::conv::ang2bohr;
                        double distSQ = (atomA_pos-atomB_pos)*(atomA_pos-atomB_pos);

                        // update NN distance and iterator
                        if ( distSQ < distNN ) {
                            distNN = distSQ;
                            NNit = bit;
                           
                        }

                    } // if ( ait != bit) 
                    i_b++;
                }// bit centers
//                cout<<"\nGridSetup Debug2 with "<<_atomgrid.size()<<" grid points \n"<<flush;
                for ( unsigned i_grid = 0; i_grid < _atomgrid.size() ; i_grid++){
//                    cout<<"i_atom="<<i_atom<<"\t_atomgrid["<<i_grid<<"].grid_weight="<<_atomgrid[i_grid].grid_weight<<endl<<flush;
                    
                    // call some shit called grid_ssw0 in NWChem
                    std::vector<double> _p = SSWpartition( i_grid, _atoms.size(), _expanded_atoms.size(), rq);
//                    cout<<"Post SSW"<<endl<<flush;
                 
                    // check weight sum
//                    cout<<"\t _p:"<<endl;
//                    int width=0;
                    double wsum = 0.0;
                    for (unsigned i =0 ; i < _p.size(); i++ ){
//                        cout<<"\t"<<_p[i];
//                        width++;
//                        if(width==10){
//                            cout<<endl<<flush;
//                            width=0;
//                        }
                        wsum += _p[i];
                    }
//                    cout<<endl<<flush;
//                    cout << " sum of partition weights " << wsum << endl;
                    if ( wsum != 0.0 ){
                        
                        // update the weight of this grid point
                            //_p is indexed by expanded_atoms
                        _atomgrid[i_grid].grid_weight = _atomgrid[i_grid].grid_weight * _p[i_atom + (_expanded_atoms.size()-_atoms.size())]/wsum;
                        //cout << " adjusting gridpoint weight "  << endl;
                    } else {
                        
                       cerr << "\nSum of partition weights of grid point " << i_grid << " of atom " << i_atom << " is zero! ";
                       throw std::runtime_error("\nThis should never happen!"); 
                        
                    }
                    
//                    cout<<"\t post ssw, grid_weight="<<_atomgrid[i_grid].grid_weight<<endl<<flush;

                } // partition weight for each gridpoint
               
                // now remove points from the grid with negligible weights
                
                for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end();) {
                    if (git->grid_weight < 1e-13 ) {
                        git = _atomgrid.erase(git);
                    } else {
                        ++git;
                    }
                }
                
                _totalgridsize += _atomgrid.size() ;

                grid.push_back(_atomgrid);
                
                i_atom++;
                
            } // atoms
            cout<<" NumericalIntegrationPeriodic::GridSetup() almost done.\n"<<flush;
            SortGridpointsintoBlocks(grid);
//            cout<<"grid boxes after SortGridpointsintoBlocks: "<<_grid_boxes.size() <<endl<<flush;
            FindSignificantShells();
//            cout<<"grid boxes after FindSignificantShells: "<<_grid_boxes.size() <<endl<<flush;
            
            
            
            
            
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
            
            
            
            
            return;
        }
        
        
        std::vector<double> NumericalIntegrationPeriodic::SSWpartition(int igrid, int ncenters, int nexpandedcenters, std::vector< std::vector<double> >& rq){
            const double ass = 0.725;
            // initialize partition vector to 1.0
            std::vector<double> p(nexpandedcenters,1.0);
            
            int nadded=(nexpandedcenters-ncenters); //number of centers added by expansion
//            cout<<"nadded= "<< nadded<< endl<< flush;
//            cout<<"igrid= "<< igrid<< endl<< flush;
            
            
            const double tol_scr = 1e-10;
            const double leps    = 1e-6; 
            // go through centers
            for ( int ibase = 0; ibase < ncenters; ibase++ ){
                
                int ij = (ibase*(ibase+1)/2 -1)+(ibase*nadded); // indexing magic
                int i = ibase + nadded;
                double rag = rq[i][igrid] ;
                
//                cout<<"\tibase= "<< ibase<< "\ti= "<< i<< "\tij_init= "<< ij+1<< endl<< flush;
                
                // through all other centers (one-directional)
                for (int j = 0; j < i ; j++ ){
                    
                    ij++;
                    if ( ( std::abs(p[i]) > tol_scr  ) || ( std::abs(p[j]) > tol_scr  ) ){
                        
                      
                        
                        double mu = ( rag - rq[j][igrid] )*Rij[ij]; 
//                        cout<<"\t\tj= "<< j<< "\tmax= "<< i<<"\tmu= "<<mu<< endl<< flush;
                        if ( mu > ass ) {
                            p[i] = 0.0;
                        } else if ( mu < -ass ) {
                            p[j] = 0.0;
                        } else {
                            
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
//                          cout<<"\tdist: "<<dist.getX()<<"\t"<<dist.getY()<<"\t"<<dist.getZ()<<"\tdistsq= "<<distsq<<"\tdecay= "<<decay<<endl<<flush;
                          // if contribution is smaller than -ln(1e-10), add atom to list
                        if ( (decay * distsq) < 20.7 ){
//                            cout<< "distsq= "<<distsq<<"\tdecay= "<<decay<<"\t(decay * distsq)="<<(decay * distsq)<<endl<<flush;
                            box.addShell(_store);
                            break;
                        }
                      }
                }
//                cout<<"box "<<i<<"\t has "<<box.Shellsize()<<" shells \t and "<<box.size()<<" points."<<endl;
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
        
        
        
        void NumericalIntegrationPeriodic::SortGridpointsintoBlocks(std::vector< std::vector< GridContainers::integration_grid > >& grid){
            const double boxsize=50;
            
            std::vector< std::vector< std::vector< std::vector< GridContainers::integration_grid* > > > >  boxes;
            
            tools::vec min=vec(std::numeric_limits<double>::max());
            tools::vec max=vec(std::numeric_limits<double>::min());
                   
            min = vec(0.0);
            max = boxLen;

            //allow for grid positions that extend from atoms in the first periodic box into nearby boxes            
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
            
             for ( auto & atomgrid : grid){
                for ( auto & gridpoint : atomgrid){
                    //gridpoint positions was already wrapped into the periodic box when it was assigned to the gridpoint
                    tools::vec pos= gridpoint.grid_pos - min;
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
        
        
        
        
        double NumericalIntegrationPeriodic::IntegrateDensity(const ub::matrix<double>& _density_matrix){
            
            double N = 0;
            
            unsigned nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif
               
               std::vector<double> N_thread=std::vector<double>(nthreads,0.0);
               
            
            
            //need to project every relevant shell onto every relevant shell in this molecule
            
               
            #pragma omp parallel for
            for (unsigned thread=0;thread<nthreads;++thread){
                for (unsigned i = thread; i < _grid_boxes.size(); i+=nthreads) {

                    double N_box=0.0;
                    GridBox& box = _grid_boxes[i];

                    //compilers do some weird memory optimizations with parts of DMAT_here
                    //resulting in garbage matrix multiplications.
                    const ub::matrix<double>  DMAT_here=box.ReadFromBigMatrix_perMolecule(_density_matrix);

                    const std::vector<tools::vec>& points=box.getGridPoints();
                    const std::vector<double>& weights=box.getGridWeights();

                    ub::range one=ub::range(0,1);

                    ub::matrix<double> ao=ub::matrix<double>(1,box.Matrixsize());
                    ub::matrix<double> ao_mol=ub::matrix<double>(1,box.Matrixsize_perMolecule());

                    box.prepareDensity();

    //                cout<<endl<<"iterate over gridpoints"<<endl<<flush;
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
    //                    cout<<"\t_temp("<< _temp.size1() <<"," << _temp.size2() <<")\tao_mol("<< ao_mol.size1() <<","<< ao_mol.size2() <<")"<<endl<<flush;
                        ub::matrix<double> pr = ub::prod(_temp, tr );
                        double rho = pr(0,0);
                        box.addDensity(rho);
                        N_box+=rho*weights[p];
                    }
                    N_thread[thread]+=N_box;

                }
            }   
            for(int i=0;i<nthreads;++i){
                N+=N_thread[i];
            }
            
               
               
            //find the number of electrons directly from Density and Overlap matrices
//            AOOverlap overlap;
            AOOverlapPeriodic overlap;
            overlap.setBox(boxLen); //in Bohr
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
//                    cout << "\t" << N_comp << "\t+=\t" << DMAT_submat(i,j)*AO_submat(i,j) <<endl<<flush;
                    N_comp += DMAT_submat(i,j)*AO_submat(i,j);
                }
            }
            
//            double N_direct=0.0;
//            for (unsigned i = 0; i < _density_matrix.size1(); i++) {
//                for (unsigned j = 0; j < _density_matrix.size2(); j++) {
//                    N_direct += _density_matrix(i,j)*AO(i,j);
//                }
//            }
            
            //check if the numbers of electrons are the same
            if(std::abs(N-N_comp)>0.005){
                cout << "N_comp from AO & DMAT is: " << N_comp << "\t and N from numerical density integration is: " << N << endl << flush;
//                cout << "N_direct is: "<<N_direct << endl << flush;
                cout <<"=======================" << endl << flush; 
                cout <<"WARNING: Calculated Densities at Numerical Grid with boxes, Number of electrons "<< N <<" is far away from the the real value "<< N_comp<<", you should increase the accuracy of the integration grid." << endl << flush; 
                cout <<"=======================" << endl << flush; 
                cout <<"#of Boxes=" << _grid_boxes.size() << endl << flush; 
                exit(-1);
            }
            
               
            density_set=true;
            return N;
        }
        

        
        

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
                
                //                #pragma omp parallel for reduction(+:result)
                for (unsigned i = 0; i < _densities.size(); i++) {
                    double q = -_densities[i] * _weights[i]; //density is neg of charge
                    vec dif = WrapDisplacement(_positions[i], rvector, boxLen);
                    dip += dif * q;
                }//i
            }//density
            else {
                throw std::runtime_error("Density not calculated");
            }
            return (abs(dip)); //in Bohr*elementary charge
        }

        
        
        
        void NumericalIntegrationPeriodic::PrepKspaceDensity_gromacs_like(double ext_alpha, std::vector< ctp::QMAtom* > & _local_atomlist, bool ECP, Grid &eval_grid, int nK) {
            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }
            
            
            if (!density_set) {
                throw std::runtime_error("Density not calculated");
            }

            //check if _periodicGridBox is already set up, if not fill it
            if(_periodicGridBox.size()==0){
                cout<<"Filling up _periodicGridBox"<<endl<<flush;
                //fill the _periodicGridBox from all the smaller grid boxes
                //The coordinates of all the gridpoints in the smaller boxes are sometimes outside the main periodic cell.
                //So, need to wrap the positions of the densities into the first periodic cell.
                for (unsigned i = 0; i < _grid_boxes.size(); i++) {
                    _periodicGridBox.appendBoxData(_grid_boxes[i]);
//                    cout<<"grid box "<<i<< " size: "<<_grid_boxes[i].size() <<endl<<flush;
                }
                
                //wrap
                _periodicGridBox.wrapPositions(boxLen);
                
                
                
                cout<<"_periodicGridBox.size() due to electrons: "<<_periodicGridBox.size()<<endl<<flush;
//                exit(0);

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
            for(int i=0; i<_densities.size(); i++){
                totQ+=_densities[i];
                totQweighted+=_densities[i]*_weights[i];
            }
            cout<<"totQ= "<<totQ<<"\ttotQweighted= "<<totQweighted<<endl<<flush;
//            exit(0);
            

            alpha = ext_alpha;
            double fourasq = 4.0 * alpha*alpha;

            //this is going to be slow, as points we have density for are not on a periodic grid,
            //so will use simple Ewald summation, not any of the FFT methods.

            //find number of k-vectors in each lattice direction
            double minSq = 1.0e15;
            for (int i = 0; i < 3; i++) {
                /*
                //numK[i]=1+(boxLen[i]/Kspacing);
                //if(numK[i]%2 == 0) numK[i]++; //keep it odd

                //find maxK
                int maxK;
                double maxKsq;
                double err;
                double twoPiL=2.0*tools::conv::Pi/boxLen[i];
                for(maxK=2; true; maxK++){
                    maxKsq=maxK*maxK*twoPiL*twoPiL;
                    err=std::exp(-maxKsq/fourasq)/maxKsq;
                    if(err<1.0e-7){
                            break;
                    }
                }

                maxK=max(16, maxK);
                
                
                numK[i]=maxK;
                minSq=min(maxKsq, minSq);
                 */
                numK[i] = nK;
            }
            cout << "numK={" << numK[0] << ", " << numK[1] << ", " << numK[2] << "}" << endl<<flush;
            cout << "rel err of k-sum ~ " << std::exp(-minSq / fourasq) / minSq << endl<<flush;


            //allocate space for eikr
            eikr.resize(3); //3 dim
            for (unsigned m = 0; m < 3; m++) { //dimension
                eikr[m].resize(numK[m]);
//#pragma omp parallel for
                for (unsigned k = 0; k < numK[m]; k++) { //k vectors
                    eikr[m][k].resize(_densities.size());    //density points
                }

                //also compute the inverse box dimensions
                lll[m] = 2.0 * boost::math::constants::pi<double>() / boxLen[m]; //units are 1/Bohr
            }
            //eikr indexing is [DIM][k][i][j]
            cout << "lll={" << lll[0] << ", " << lll[1] << ", " << lll[2] << "}" << endl<<flush;

            //based on gromacs 4.6 tabulate_eir()
//#pragma omp parallel for
            for (unsigned i = 0; i < _densities.size(); i++) { //density points
                vec r = _positions[i]; //grid is in Bohr, need it in reduced units (*2pi/box))
                for (unsigned m = 0; (m < 3); m++) { //dimensions
                    eikr[m][0][i] = std::complex<double>(1, 0);
                    eikr[m][1][i] = std::complex<double>(cos(r[m] * lll[m]), sin(r[m] * lll[m])); //this is where reduced units are applied
                    for (unsigned k = 2; k < numK[m]; k++) { //k vectors
                        eikr[m][k][i] = eikr[m][k - 1][i] * eikr[m][1][i];
                    }//k
                }//m
            }//i


            //now do the same for the eval_grid (eikR), the points were potential is to be evaluated
            //allocate space for eikR
            eikR.resize(3); //3 dim
            for (unsigned m = 0; m < 3; m++) { //dimensions
                eikR[m].resize(numK[m]);
                for (unsigned k = 0; k < numK[m]; k++) {
                    eikR[m][k].resize(eval_grid.getGrid().size());
                }
            }
            //eikR indexing is [DIM][k][i]

            //based on gromacs 4.6 tabulate_eir()
//#pragma omp parallel for
            //printf("eval grid size = %d\n", eval_grid.getGrid().size());
            //printf("lll_x = %f\n", lll[0]);
            for (unsigned i = 0; i < eval_grid.getGrid().size(); i++) {
                vec R = eval_grid.getGrid()[i] * tools::conv::nm2bohr; //this is in Bohr, eval_grid is in nm
                for (unsigned m = 0; (m < 3); m++) {
//                    cout<<"m= "<<m<<endl<<flush;
                    eikR[m][0][i] = std::complex<double>(1, 0);
                    eikR[m][1][i] = std::complex<double>(cos(R[m] * lll[m]), sin(R[m] * lll[m]));
                    //unsigned k=0;
                    //printf("%d\t%f+i%f\t R_x=%f\n", k, eikR[m][k][i].real(),eikR[m][k][i].imag(), R[m]);
                    //k=1;
                    //printf("%d\t%f+i%f\t R_x=%f\n", k, eikR[m][k][i].real(),eikR[m][k][i].imag(), R[m]);
                    for (unsigned k = 2; k < numK[m]; k++) {
                        eikR[m][k][i] = eikR[m][k - 1][i] * eikR[m][1][i];
//                        printf("\t%d\t%f+i%f\t R_x=%f\n", k, eikR[m][k][i].real(),eikR[m][k][i].imag(), R[m]);
                    }//k
                    //exit(0);
                }//m
//                exit(0);
            }//i
//            exit(0);
        }

//        void NumericalIntegrationPeriodic::FillMadelungGrid(vec box, int natomsonside) {
//            //fill _Madelung_grid;
//            _Madelung_grid.clear();
//            std::vector< GridContainers::integration_grid > _Mad;
//            boxLen = box;
//            double a = boxLen[0]; //in bohr
//            for (int l = 0; l < 2; l++) {
//                for (int m = 0; m < 2; m++) {
//                    for (int n = 0; n < 2; n++) {
//                        GridContainers::integration_grid el;
//                        el.grid_density = std::pow(-1.0, l + m + n);
//                        el.grid_weight = 1.0;
//                        el.grid_pos = vec(l, m, n) * (a / 2);
//                        _Mad.push_back(el);
//                        //cout<< "q= "<< el.grid_density << "\t @ " << el.grid_x << " " << el.grid_y << " " << el.grid_z << "\t box size:"<< boxLen[0] << endl;
//                    }
//                }
//            }
//            _Madelung_grid.push_back(_Mad);
//            _grid = _Madelung_grid;
//            //exit(0);
//        }
//
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
                
                
                #pragma omp parallel for
                for (unsigned p = 0; p < _ESPatGrid.size(); p++) {
                    _ESPatGrid[p] = 0;
                    vec rvector = eval_grid.getGrid()[p] * tools::conv::nm2bohr; //Bohr
                    for (unsigned i = 0; i < _densities.size(); i++) {

                            //charge at this point
                            double q = -_weights[i] * _densities[i]; //density is neg of charge
                            //double q = -_Madelung_grid[i][j].grid_weight * _Madelung_grid[i][j].grid_density; //density is neg of charge

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
                E_rspace = _ESPatGrid(0);
                E_erfc = erfc(alpha);


                //adapted from gromacs 4.6 do_ewald_pot()

                //zero pot
                //                for (unsigned p = 0; p < _ESPatGrid.size(); p++) {
                //                    _ESPatGrid[p] = 0;
                //                }



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


                //printf("factor = %f\n", factor);
                //printf("scaleRecip = %f\n", scaleRecip);
                //printf("vol = %f\n", vol);
                //printf("actual box for Madelung: %f\t%f\t%f\n", boxLen[0], boxLen[1], boxLen[2]);
                //exit(0);


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
                                //                                    if(ix==4 && iy==4 && iz==4){
                                //                                        printf("cs+= %f = %f * (%f+i%f) * (%f+i%f)\n", tab_qxyz[i][j].real(), -_grid[i][j].grid_weight * _grid[i][j].grid_density, tab_xy[i][j].real(), tab_xy[i][j].imag(), eikr[2][iz][i][j].real(), eikr[2][iz][i][j].imag());
                                //                                        printf("\teikr=(%f+i%f) * (%f+i%f) * (%f+i%f)\n",eikr[0][ix][i][j].real(), eikr[0][ix][i][j].imag(),eikr[1][iy][i][j].real(), eikr[1][iy][i][j].imag(),eikr[2][iz][i][j].real(), eikr[2][iz][i][j].imag());
                                //                                        printf("\tr=%f\t%f\t%f\n", _grid[i][j].grid_pos[0],_grid[i][j].grid_pos[1], _grid[i][j].grid_pos[2]);
                                //                                        std::complex<double> krx = std::complex<double>(0, -_grid[i][j].grid_pos[0]*mz);
                                //                                        krx=exp(krx);
                                //                                        printf("\tkr_x=exp(-i%f)=%f+i%f\n",_grid[i][j].grid_pos[0]*mz, krx.real(), krx.imag());
                                //                                    }
                            }
                            //exit(0);
                            //                            if(abs(cs)>1e-6 || abs(ss)>1e-6)
                            //                            if(ix==4 && iy==4 && iz==4)
                            //                                printf("k=[%d,%d,%d]\tcs=%f\t ss=%f \tak=%f=exp(%f * %f) / %f\n", ix, iy, iz, cs, ss, ak, m2, factor, m2);
                            for (int n = 0; n < eval_grid.getsize(); n++) {
                                tmp = ak * (cs * tab_R_xyz[n].real() + ss * tab_R_xyz[n].imag());
                                //                                if((abs(cs)>1e-6 || abs(ss)>1e-6) && n==0){
                                //                                    printf("ak=%f=exp(%f * %f) / %f\n", ak, m2, factor, m2);
                                //                                    printf("tab_R_xyz=%f+i%f= (%f+i%f) * (%f+i%f) * (%f+i%f)\n", tab_R_xyz[n].real(), tab_R_xyz[n].imag(), eikR[0][ix][n].real(), eikR[0][ix][n].imag(), eikR[1][abs(iy)][n].real(), eikR[1][abs(iy)][n].imag(), eikR[2][abs(iz)][n].real(), eikR[2][abs(iz)][n].imag());
                                //                                    printf("\ttmp=%f=%f*(%f*%f + %f*%f)\t m2=%f\n", tmp, ak, cs, tab_R_xyz[n].real(), ss, tab_R_xyz[n].imag(), m2);
                                //                                    printf("\tk=[%d,%d,%d] of [%d,%d,%d]\ttmp=%f \t sum=%f\n", ix, iy, iz, numK[0], numK[1], numK[2], tmp, k_pot(n));
                                //                                }
                                _ESPatGrid(n) += tmp * 2 * scaleRecip;
                                k_pot(n) += tmp * 2 * scaleRecip;

                            }

                        }//iz
                        lowiz = 1 - numK[2];
                        //                        if(ix==1)
                        //                            exit(0);                        
                    }//iy
                    lowiy = 1 - numK[1];
                }//ix
                //printf("%f=%f+%f\n", _ESPatGrid(0), E_rspace, k_pot(0));
                //                exit(0);


            } else {
                throw std::runtime_error("Density not calculated");
            }
            E_kspace = k_pot(0);
            //E_kspace = _ESPatGrid(1);
            //_ESPatGrid(1)+=E_rspace;
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
