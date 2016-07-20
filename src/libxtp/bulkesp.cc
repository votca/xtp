
#include <votca/xtp/bulkesp.h>
#include <votca/xtp/aomatrix.h>
#include <votca/tools/linalg.h>
//#include <boost/progress.hpp>
#include <votca/xtp/numerical_integrations.h>
#include <math.h> 
#include <votca/tools/constants.h>

#include "votca/xtp/orbitals.h"

using namespace votca::tools;


namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
    std::vector<Bulkesp::Molecule> Bulkesp::BreakIntoMolecules(std::vector< QMAtom* > a, double scale){
        
        std::vector<Bulkesp::Molecule> mols;
        std::list<Bulkesp::Bond> bonds;
        
        Elements _elements;
        
        LOG(logDEBUG, *_log) << " BreakIntoMolecules(): locating bonds.\n" << flush;
        
        //find all the bonds;
        for (vector<QMAtom*>::iterator i = a.begin(); i != a.end(); ++i){
            for(vector<QMAtom*>::iterator j = i+1; j != a.end(); ++j){
                double dif[3];
                dif[0] = ((*i)->x - (*j)->x);
                dif[1] = ((*i)->y - (*j)->y);
                dif[2] = ((*i)->z - (*j)->z);
                for(int k=0; k<3; k++){
                    if(periodic && std::abs(dif[k])>boxLen[k]*0.5) //correct for for bond crossing PBC, if it exists
                        if(dif[k]>0)    //i.x>j.x
                            dif[k]-=boxLen[k];
                        else            //i.x<j.x
                            dif[k]+=boxLen[k];
                }
                vec v(dif);
                double distSq = v*v;
                
                double acceptDist=_elements.getCovRad((*i)->type) + _elements.getCovRad((*j)->type);
                if(distSq<=acceptDist*acceptDist*scale){
                    Bond nb;        //bond goes
                    nb.a=(*i);      //from a
                    nb.b=(*j);      //to b
                    nb.ba=v;        //and has a (PBC corrected) vector ba
                    bonds.push_back(nb);
                }
            }
        }
        LOG(logDEBUG, *_log) << " BreakIntoMolecules(): "<< bonds.size() <<" bonds found.\n" << flush;
  
        
        //now add all the atoms that have no bonds as a separate molecule
        Bulkesp::Molecule leftover;
        for (vector<QMAtom*>::iterator i = a.begin(); i != a.end(); ++i){
            bool found = false; //does this atom have bonds?
            for (std::list<Bond>::iterator b = bonds.begin(); b != bonds.end(); ++b){
                if(b->a==*i || b->b==*i){
                    found=true;
                    break;
                }
            }
            if(!found) //atom has no bonds
                leftover.atoms.push_back(*i);
        }
        mols.push_back(leftover);
        
        
        
        //assign bonds to molecules
        while(bonds.size()>0){
            //init the new molecule with the last bond we have
            Bulkesp::Molecule m;
            m.atoms.push_back(bonds.end()->a);
            if(periodic){//unwrap system by moving b
                bonds.end()->b->x = bonds.end()->a->x + bonds.end()->ba.getX();
                bonds.end()->b->y = bonds.end()->a->y + bonds.end()->ba.getY();
                bonds.end()->b->z = bonds.end()->a->z + bonds.end()->ba.getZ();
            }
            m.atoms.push_back(bonds.end()->b);
            bonds.pop_back();
            
            
            int oldMsize;
            do{
                oldMsize=m.atoms.size();
                //loop through bonds to see what binds to what is already in the molecule
                for (std::list<Bond>::iterator b = bonds.begin(); b != bonds.end(); ++b){
                    if(find (m.atoms.begin(), m.atoms.end(), b->a)!=m.atoms.end()) //molecule contains a
                        //double check if molecule contains b already (don't double count atoms in circular molecules)
                        if(find (m.atoms.begin(), m.atoms.end(), b->b)==m.atoms.end()){ // contains a, but not b
                            if(periodic){//unwrap system by moving b
                                b->b->x = bonds.end()->a->x + b->ba.getX();
                                b->b->y = bonds.end()->a->y + b->ba.getY();
                                b->b->z = bonds.end()->a->z + b->ba.getZ();
                            }
                            m.atoms.push_back(b->b); //add b to molecule
                            b=bonds.erase(b);
                            --b;
                        }
                    else if(find (m.atoms.begin(), m.atoms.end(), b->b)!=m.atoms.end()) //molecule contains b
                        //double check if molecule contains a already (don't double count atoms in circular molecules)
                        if(find (m.atoms.begin(), m.atoms.end(), b->a)==m.atoms.end()){ // contains b, but not a
                            if(periodic){//unwrap system by moving a
                                b->a->x = b->b->x - b->ba.getX();
                                b->a->y = b->b->y - b->ba.getY();
                                b->a->z = b->b->z - b->ba.getZ();
                            }
                            m.atoms.push_back(b->a); //add a to molecule
                            b=bonds.erase(b);
                            --b;
                        }

                }
            }while(oldMsize!=m.atoms.size());   //keep looping until nothing is added to m
            
            //we are done with this molecule, save it to mols
            mols.push_back(m);
            LOG(logDEBUG, *_log) << " BreakIntoMolecules(): put "<< m.atoms.size() 
                    << " atoms into molecule "<< mols.size() <<";\t" << bonds.size() <<" bonds left.\n" << flush;
        }
        
        if(periodic)
            LOG(logDEBUG, *_log) << " BreakIntoMolecules(): Molecules have been unwrapped.\n" << flush;
        
        return(mols);
    }
    
    
    
    ub::matrix<double> Bulkesp::BuildDenMat(Orbitals &_orb, std::string _state, std::string _spin, int _state_no)
    {
        ub::matrix<double> DMAT_tot;
        bool _do_transition=false;
        
        if(_state!="ground")
            LOG(logDEBUG, *_log) << " Bulkesp is not tested for any state other than the ground state. "
                    <<"It will likely not work correctly, as BSE Singlet and Triplet Coefficients are simply copied from the global orbital object.\n"  << flush; 
        
        if(_state=="transition"){
            _do_transition=true;
            if (_spin=="singlet"){
                DMAT_tot=_orb.TransitionDensityMatrix(_orb.MOCoefficients(), _orb.BSESingletCoefficients(), _state_no-1);
            }
            else if (_spin=="triplet"){
                DMAT_tot=_orb.TransitionDensityMatrix(_orb.MOCoefficients(), _orb.BSETripletCoefficients(), _state_no-1); 
            }
            else throw std::runtime_error("Spin entry not recognized");
        }
        else if (_state=="ground" || _state=="excited"){
            
        
            ub::matrix<double> &DMATGS=_orb.DensityMatrixGroundState(_orb.MOCoefficients());
            DMAT_tot=DMATGS;
            if ( _state_no > 0 && _state=="excited"){
                std::vector<ub::matrix<double> > DMAT;
                if (_spin=="singlet"){
                    DMAT = _orb.DensityMatrixExcitedState(_orb.MOCoefficients() , _orb.BSESingletCoefficients(), _state_no-1);
                }
                else if (_spin=="triplet"){
                    DMAT = _orb.DensityMatrixExcitedState( _orb.MOCoefficients() , _orb.BSETripletCoefficients(), _state_no-1);
                }
                else throw std::runtime_error("Spin entry not recognized");
                DMAT_tot=DMAT_tot-DMAT[0]+DMAT[1];
            }            
	   // Ground state + hole_contribution + electron contribution
	}
        else throw std::runtime_error("State entry not recognized");
        
        
        return(DMAT_tot);
    }
    
    
    
    
    ub::matrix<double> Bulkesp::BuildOverlapMat(Orbitals &_molOrb, Orbitals &_globalOrb){
        ub::matrix<double> ov;
        return (ov);
    }
    
    
    
    
    ub::vector<double> Bulkesp::ComputeESP(std::vector< QMAtom* > & _atomlist, ub::matrix<double> &_dmat,
            ub::matrix<double> &_ovmat, AOBasis &_basis,BasisSet &bs,string gridsize, Grid &_grid){
        // Calculating nuclear potential at gridpoints

        ub::vector<double> _ESPatGrid = ub::zero_vector<double>(_grid.getsize());
        // ub::vector<double> _NucPatGrid = ub::zero_vector<double>(_gridpoints.size());

        
#warning "TODO: Need to check if overlap already exists (overlap from CPMD is better because it takes PBC into account)"
        ub::vector<double> DMATasarray=_dmat.data();
        ub::vector<double> AOOasarray=_ovmat.data();
        double N_comp=0.0;
        #pragma omp parallel for reduction(+:N_comp) 
        for ( unsigned _i =0; _i < DMATasarray.size(); _i++ ){
                N_comp =N_comp+ DMATasarray(_i)*AOOasarray(_i);
            } 


        NumericalIntegration numway;

        numway.GridSetup(gridsize,&bs,_atomlist);
        LOG(logDEBUG, *_log) << TimeStamp() << " Calculate Densities at Numerical Grid with gridsize "<<gridsize  << flush; 
        double N=numway.IntegrateDensity_Atomblock(_dmat,&_basis);
        LOG(logDEBUG, *_log) << TimeStamp() << " Calculated Densities at Numerical Grid, Number of electrons is "<< N << flush; 

        if(std::abs(N-N_comp)>0.001){
            LOG(logDEBUG, *_log) <<"=======================" << flush; 
            LOG(logDEBUG, *_log) <<"WARNING: Calculated Densities at Numerical Grid, Number of electrons "
                    << N <<" is far away from the the real value "<< N_comp
                    <<", you should increase the accuracy of the integration grid."<< flush; 
            N=N_comp;
            LOG(logDEBUG, *_log) <<"WARNING: Electronnumber set to "<< N << flush; 
            LOG(logDEBUG, *_log) <<"=======================" << flush; 
        }

        double netcharge=getNetcharge( _atomlist,N );   

        LOG(logDEBUG, *_log) << TimeStamp() << " Calculating ESP at CHELPG grid points"  << flush;     
        //boost::progress_display show_progress( _grid.getsize() );
        #pragma omp parallel for
        for ( int i = 0 ; i < _grid.getsize(); i++){
            _ESPatGrid(i)=numway.IntegratePotential(_grid.getGrid()[i]*tools::conv::nm2bohr);
            //++show_progress;
        }

        LOG(logDEBUG, *_log) << TimeStamp() << " Electron contribution calculated"  << flush; 
        ub::vector<double> _NucPatGrid = EvalNuclearPotential(  _atomlist,  _grid );
        _ESPatGrid += _NucPatGrid;
        LOG(logDEBUG, *_log) << TimeStamp() << " Nuclear contribution calculated"  << flush; 

        
        return(_ESPatGrid);
    }
    
    
    
    void Bulkesp::FillElement2NBF(std::vector< QMAtom* >& _atomlist, BasisSet &bs){
        
        //Find list of elements
        _elements.clear();
        for (std::vector< QMAtom* >::iterator ait = _atomlist.begin(); ait < _atomlist.end(); ait++) {

            std::string element_name = (*ait)->type;
            list<std::string>::iterator ite;
            ite = find(_elements.begin(), _elements.end(), element_name);
            if (ite == _elements.end()) {            //this is the first atom of this element encountered
                _elements.push_back(element_name);
            }
        }
        
        //Find NBF for each element
        _element2NBF.clear();
        for(list<std::string>::iterator eit = _elements.begin(); eit!=_elements.end(); ++eit){
            vector<Shell*> shells = bs.getElement(*eit)->getShells();
            int nfunc=0;
            //loop over shells
            for(vector< Shell* >::iterator sit = shells.begin(); sit!=shells.end(); ++sit){
                nfunc += 1 + 2 * ((*sit)->getLmax());
            }
            //store number of functions
            _element2NBF[(*eit)]=nfunc;
        }
        
    }
    
    
    
    //Find indeces of MO coefficients of the first shell of each atom,
    //assuming that the order of _atomlist is the same as that of the _MO_Coefficients.
    std::map<QMAtom*,int> Bulkesp::MapAtom2MOCoefIndex(std::vector< QMAtom* >& _atomlist){
            std::map<QMAtom*,int> ret;
            int count=0;
            //loop over all atoms in system
            for(vector<QMAtom*>::iterator i = _atomlist.begin(); i != _atomlist.end(); ++i){
                ret[(*i)] = count;
                count += _element2NBF[(*i)->type];
            }
            return(ret);
    }
    
    
                
    void Bulkesp::Evaluate(std::vector< QMAtom* >& _atomlist, Orbitals& _globalOrb,
            ub::matrix<double> _global_MO_Coeffs, AOBasis &_basis, BasisSet &bs,
            string gridsize, double maxBondScale, std::string _state,
            std::string _spin, int _state_no){
        
        //find the individual molecules
        std::vector<Bulkesp::Molecule> mols = BreakIntoMolecules(_atomlist, maxBondScale);
        
        //map elements to number of basis functions per element
        FillElement2NBF(_atomlist, bs);
        //map QMAtom pointers to column indeces in _MO_Coefficients
        std::map<QMAtom*,int> atom2MOIndex = MapAtom2MOCoefIndex(_atomlist);
        
        //loop over molecules
        LOG(logDEBUG, *_log) << " Bulkesp::Evaluate(): found "<< mols.size() << "molecules.\n" << flush; 
        for (std::vector<Bulkesp::Molecule>::iterator m = mols.begin(); m != mols.end(); ++m){
            
            LOG(logDEBUG, *_log) << " Bulkesp::Evaluate(): "<< TimeStamp()<<" processing molecule "<< m-mols.begin() << endl << flush; 
            
            //Obtain AO data relevant to this molecule only.
            //extract basis
            AOBasis _m_basis;
            _m_basis.AOBasisFill(&bs, m->atoms);
            
            //extract MO coefficients
            int numorb=_global_MO_Coeffs.size1();
            ub::matrix<double> _m_coefs;
            _m_coefs.resize(numorb,_m_basis.AOBasisSize());
            int bfgi; //basis function global index
            int bfli=0; //basis function local index
            //loop over atoms in this molecule
            for(vector<QMAtom*>::iterator a = _atomlist.begin(); a != _atomlist.end(); ++a){
                bfgi=atom2MOIndex[(*a)];
                for(int i=0; i<_element2NBF[(*a)->type]; i++){//loop over basis functions of this atom
                    //copy MO coefficients to the molecule's matrix
                    ub::column(_m_coefs, bfli) = ub::column(_global_MO_Coeffs, bfgi);
                    bfli++;
                    bfgi++;
                }
            }
            if(bfli!=_m_basis.AOBasisSize())    //check number of basis functions
                throw std::runtime_error("Number of basis functions in molecule does not match that in its basis set.\n");
            
            //Set up the Orbital object for this molecule and calculate density matrix
            Orbitals _molOrb;
            _molOrb.MOCoefficients()=_m_coefs;
            _molOrb.setQMpackage("votca"); //MO coefs are already in votca's ordering
            _molOrb.setBasisSetSize(_m_basis.AOBasisSize());
            _molOrb.setBSEindices(_globalOrb.getBSEvmin(), _globalOrb.getBSEvmin(),
                    _globalOrb.getBSEcmin(), _globalOrb.getBSEcmax(), 0);
            _molOrb.BSESingletCoefficients()=_globalOrb.BSESingletCoefficients();
            _molOrb.BSETripletCoefficients()=_globalOrb.BSETripletCoefficients();
            
            //set up the density matrix
            ub::matrix<double> _m_dmat=BuildDenMat(_molOrb, _state, _spin, _state_no);
            
            //set up the overlap matrix
#warning "TODO: Need to check if overlap already exists (overlap from CPMD is better because it takes PBC into account)"
            ub::matrix<double> _m_ovmat;
            if(periodic)
                if(_globalOrb.hasAOOverlap())
                    _m_ovmat = _globalOrb.AOOverlap();
                else throw std::runtime_error("Periodic system, but periodic AO overlaps not present in Orbitals.");
            else{
                AOOverlap overlap;
                overlap.Initialize(_basis._AOBasisSize);
                overlap.Fill(&_basis);
                _m_ovmat = overlap._aomatrix;
            }
            
            //set up grid
            // setting up grid    
            Grid _grid;
            _grid.setAtomlist(&_atomlist);
            _grid.setupCHELPgrid();
            LOG(logDEBUG, *_log) << TimeStamp() <<  " Done setting up CHELPG grid with " << _grid.getsize() << " points " << endl;
            
            //calculate the ESP
            ub::vector<double> ESP=ComputeESP(m->atoms, _m_dmat, _m_ovmat, _m_basis, bs, gridsize, _grid);
            
            //output
            //CHELPG grids aren't periodic and equally spaced,
            //so can't output to .cube format.
            //Create own format.
            std::ostringstream fn;
            fn << "BulkEsp_" << m-mols.begin() << ".grid";
            _grid.writeIrregularGrid(fn.str(), ESP, m->atoms);
            
            //TODO: fit charges
            
        }
        LOG(logDEBUG, *_log) << " Bulkesp::Evaluate(): "<< TimeStamp()<<" All molecules processed." << endl << flush; 
        
    }
    
    
}}