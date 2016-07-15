
#include <votca/xtp/bulkesp.h>
#include <votca/xtp/aomatrix.h>
#include <votca/tools/linalg.h>
//#include <boost/progress.hpp>
#include <votca/xtp/numerical_integrations.h>
#include <math.h> 
#include <votca/tools/constants.h>

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
                double dif = ((*i)->x - (*j)->x);
                double dist = dif*dif;
                dif = ((*i)->y - (*j)->y);
                dist += dif*dif;
                dif = ((*i)->z - (*j)->z);
                dist += dif*dif;
                
                double aceptDist=_elements.getCovRad((*i)->type) + _elements.getCovRad((*j)->type);
                if(dist<=aceptDist*aceptDist*scale){
                    Bond nb;
                    nb.a=(*i);
                    nb.b=(*j);
                    bonds.push_back(nb);
                }
            }
        }
        LOG(logDEBUG, *_log) << " BreakIntoMolecules(): "<< bonds.size() <<" bonds found.\n" << flush;
  
        //assign bonds to molecules
        while(bonds.size()>0){
            //init the new molecule with the last bond we have
            Bulkesp::Molecule m;
            m.atoms.push_back(bonds.end()->a);
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
                            m.atoms.push_back(b->b); //add b to molecule
                            b=bonds.erase(b);
                            --b;
                        }
                    else if(find (m.atoms.begin(), m.atoms.end(), b->b)!=m.atoms.end()) //molecule contains b
                        //double check if molecule contains a already (don't double count atoms in circular molecules)
                        if(find (m.atoms.begin(), m.atoms.end(), b->a)==m.atoms.end()){ // contains b, but not a
                            m.atoms.push_back(b->a); //add a to molecule
                            b=bonds.erase(b);
                            --b;
                        }

                }
            }while(oldMsize!=m.atoms.size());   //keep looping until nothing is added to m
            
            //we are done with this molecule, save it to mols
            mols.push_back(m);
            LOG(logDEBUG, *_log) << " BreakIntoMolecules(): put "<< m.atoms.size() << " atoms into molecule "<< mols.size() <<";\t" << bonds.size() <<" bonds left.\n" << flush;
        }
        return(mols);
    }
    
    
    
    ub::vector<double> Bulkesp::ComputeESP(std::vector< QMAtom* > & _atomlist, ub::matrix<double> &_dmat, AOBasis &_basis,BasisSet &bs,string gridsize){
        // setting up grid    
        Grid _grid;
        _grid.setAtomlist(&_atomlist);
        _grid.setupCHELPgrid();
        //_grid.printGridtoxyzfile("grid.xyz");
        LOG(logDEBUG, *_log) << TimeStamp() <<  " Done setting up CHELPG grid with " << _grid.getsize() << " points " << endl;

        // Calculating nuclear potential at gridpoints

        ub::vector<double> _ESPatGrid = ub::zero_vector<double>(_grid.getsize());
        // ub::vector<double> _NucPatGrid = ub::zero_vector<double>(_gridpoints.size());

        
#warning "TODO: Need to check if overlap already exists (overlap from CPMD is better because it takes PBC into account)"
        AOOverlap overlap;
        overlap.Initialize(_basis._AOBasisSize);
        overlap.Fill(&_basis);
        ub::vector<double> DMATasarray=_dmat.data();
        ub::vector<double> AOOasarray=overlap._aomatrix.data();
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
            LOG(logDEBUG, *_log) <<"WARNING: Calculated Densities at Numerical Grid, Number of electrons "<< N <<" is far away from the the real value "<< N_comp<<", you should increase the accuracy of the integration grid."<< flush; 
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
    
    
                
    void Bulkesp::Evaluate(std::vector< QMAtom* >& _atomlist, ub::matrix<double> _MO_Coefficients, AOBasis &_basis,BasisSet &bs,string gridsize, double maxBondScale){
    
        //find the individual molecules
        std::vector<Bulkesp::Molecule> mols = BreakIntoMolecules(_atomlist, maxBondScale);
        
        //map elements to number of basis functions per element
        FillElement2NBF(_atomlist, bs);
        //map QMAtom pointers to column indeces in _MO_Coefficients
        std::map<QMAtom*,int> atom2MOIndex = MapAtom2MOCoefIndex(_atomlist);
        
        //loop over molecules
        for (std::vector<Bulkesp::Molecule>::iterator m = mols.begin(); m != mols.end(); ++m){
            //extract AO data relevant to this molecule only
                //basis
            AOBasis _m_basis;
            _m_basis.AOBasisFill(&bs, m->atoms);
            
                //MO coefficients
            int numorb=_MO_Coefficients.size1();
            ub::matrix<double> _m_coefs;
            _m_coefs.resize(numorb,_m_basis.AOBasisSize());
            int bfgi; //basis function global index
            int bfli=0; //basis function local index
            //loop over atoms in this molecule
            for(vector<QMAtom*>::iterator a = _atomlist.begin(); a != _atomlist.end(); ++a){
                bfgi=atom2MOIndex[(*a)];
                for(int i=0; i<_element2NBF[(*a)->type]; i++){//loop over basis functions of this atom
                    //copy MO coefficients to the molecule's matrix
                    ub::column(_m_coefs, bfli) = ub::column(_MO_Coefficients, bfgi);
                    bfli++;
                    bfgi++;
                }
            }
            if(bfli!=_m_basis.AOBasisSize())    //check number of basis functions
                throw std::runtime_error("Number of basis functions in molecule does not match that in its basis set.\n");
            
                //TODO: calculate density matrix
            ub::matrix<double> _m_dmat;
            
            
            //calculate the ESP
            ub::vector<double> ESP=ComputeESP(m->atoms, _m_dmat, _m_basis, bs, gridsize);
            
            //TODO: output the ESP to a cube file
        }
        
    }
    
}}