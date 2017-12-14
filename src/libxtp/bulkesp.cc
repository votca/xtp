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


#include <votca/xtp/bulkesp.h>
#include <votca/xtp/aomatrix.h>
#include <votca/tools/linalg.h>
//#include <boost/progress.hpp>
#include <votca/xtp/numerical_integrations_periodic.h>
#include <math.h> 
#include <votca/tools/constants.h>

#include "votca/xtp/orbitals.h"

#include <fstream>

using namespace votca::tools;


namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;
        

        std::vector<Bulkesp::Molecule> Bulkesp::BreakIntoMolecules(std::vector< ctp::QMAtom* > _atoms, double scale) {

            std::vector<Bulkesp::Molecule> mols;
            std::list<Bulkesp::Bond> bonds;

            Elements _elements;

            CTP_LOG(ctp::logDEBUG, *_log) << " BreakIntoMolecules(): locating bonds.\n" << flush;
            
            //find all the bonds;
            for (vector<ctp::QMAtom*>::iterator i = _atoms.begin(); i != _atoms.end(); ++i) {
                for (vector<ctp::QMAtom*>::iterator j = i + 1; j != _atoms.end(); ++j) {
                    tools::vec dif;
                    if(periodic){
//                        dif = WrapPoint((*i)->getPos(), boxLen) - WrapPoint((*j)->getPos(), boxLen); //Wrap points into periodic box
                        dif = (*i)->getPos() - (*j)->getPos();
                    }
                    else{
                        dif = (*i)->getPos() - (*j)->getPos();
                    }

                    //cout << "dif: " << dif[0] << " " << dif[1] << " " << dif[2] << "\t";
                    for (int k = 0; k < 3; k++) {
                        if (periodic && std::abs(dif[k]) > boxLen[k]*0.5) //correct for for bond crossing PBC, if it exists
                            if (dif[k] > 0) //i.x>j.x
                                dif[k] -= boxLen[k];
                            else //i.x<j.x
                                dif[k] += boxLen[k];
                    }
                    //cout << "PBC corrected: " << dif[0] << " " << dif[1] << " " << dif[2] << "\n";
                    vec v(dif);
                    double distSq = v*v;

                    double acceptDist = _elements.getCovRad((*i)->type) + _elements.getCovRad((*j)->type);
//                    if(fabs((*i)->x - 14.905)<0.002 && fabs((*j)->x - 14.266)<0.002){
//                    //if(fabs((*i)->x - 14.266)<0.002){
//                        cout<< (*i)->type << "\t" << (*i)->x << " " << (*i)->y << " " << (*i)->z <<endl;
//                        cout<< (*j)->type << "\t" << (*j)->x << " " << (*j)->y << " " << (*j)->z <<endl;
//                        cout<<"\tacceptDist:"<<acceptDist<<"\tscale:"<<scale<< "\tdist:"<<sqrt(distSq)<<endl;
//                        cout << "\tPBC corrected dist:" <<distSq <<endl<<flush;
//                        //exit(0);
//                    }
                    if (distSq <= acceptDist * acceptDist * scale * scale) {
                        Bond nb; //bond goes
                        nb.a = (*i); //from a
                        nb.b = (*j); //to b
                        nb.ba = v; //and has a (PBC corrected) vector ba (from b to a))
                        nb.a_indx = i - _atoms.begin();
                        nb.b_indx = j - _atoms.begin();
                        bonds.push_back(nb);
                    }
                }
            }
            CTP_LOG(ctp::logDEBUG, *_log) << " BreakIntoMolecules(): " << bonds.size() << " bonds found.\n" << flush;
            //        for (std::list<Bond>::iterator b = bonds.begin(); b != bonds.end(); ++b){
            //            CTP_LOG(ctp::logDEBUG, *_log) << b->a_indx <<" - "<<b->b_indx<<"\t"<<b->ba<<"\t"<<sqrt(b->ba*b->ba)<< flush;
            //        }


            //now add all the atoms that have no bonds as a separate molecule
            Bulkesp::Molecule leftover;
            for (vector<ctp::QMAtom*>::iterator i = _atoms.begin(); i != _atoms.end(); ++i) {
                bool found = false; //does this atom have bonds?
                for (std::list<Bond>::iterator b = bonds.begin(); b != bonds.end(); ++b) {
                    if (b->a == *i || b->b == *i) {
                        found = true;
                        break;
                    }
                }
                if (!found) //atom has no bonds
                    leftover.atoms.push_back(*i);
                leftover.atomIndeces.push_back(i - _atoms.begin());
            }
            if (leftover.atoms.size() > 0) {
                mols.push_back(leftover);
                CTP_LOG(ctp::logDEBUG, *_log) << " BreakIntoMolecules(): put " << leftover.atoms.size()
                        << " Unbonded atoms into molecule " << mols.size() - 1 << ";\t" << bonds.size() << " bonds left.\nThey were:" << flush;
                for(std::vector<ctp::QMAtom*>::iterator it = leftover.atoms.begin(); it != leftover.atoms.end(); ++it) {
                    /* std::cout << *it; ... */
                    CTP_LOG(ctp::logDEBUG, *_log) << (*it)->type << "\t" << (*it)->x << "\t" << (*it)->y << "\t" << (*it)->z << endl << flush;
                }
                exit(-1);
                
            }



            //assign bonds to molecules
            while (bonds.size() > 0) {
                //init the new molecule with the last bond we have
                Bulkesp::Molecule m;
                std::list<Bond>::iterator last = bonds.end();
                std::advance(last, -1);
                m.atoms.push_back(last->a);
                m.atomIndeces.push_back(last->a_indx);
                leftover.atomIndeces.push_back(bonds.end()->a_indx);
                if (periodic) {//unwrap system by moving b
                    last->b->x = last->a->x - last->ba.getX();
                    last->b->y = last->a->y - last->ba.getY();
                    last->b->z = last->a->z - last->ba.getZ();
                }
                m.atoms.push_back(last->b);
                m.atomIndeces.push_back(last->b_indx);
                bonds.pop_back();


                int oldMsize;
                do {
                    oldMsize = m.atoms.size();
                    //loop through bonds to see what binds to what is already in the molecule
                    for (std::list<Bond>::iterator b = bonds.begin(); b != bonds.end(); ++b) {
                        if (find(m.atoms.begin(), m.atoms.end(), b->a) != m.atoms.end()) { //molecule contains a
                            //double check if molecule contains b already (don't double count atoms in circular molecules)
                            if (find(m.atoms.begin(), m.atoms.end(), b->b) == m.atoms.end()) { // contains a, but not b
                                if (periodic) {//unwrap system by moving b
                                    b->b->x = b->a->x - b->ba.getX();
                                    b->b->y = b->a->y - b->ba.getY();
                                    b->b->z = b->a->z - b->ba.getZ();
                                }
                                m.atoms.push_back(b->b); //add b to molecule
                                m.atomIndeces.push_back(b->b_indx);
                                b = bonds.erase(b);
                                --b;
                            }
                        } else if (find(m.atoms.begin(), m.atoms.end(), b->b) != m.atoms.end()) { //molecule contains b
                            //double check if molecule contains a already (don't double count atoms in circular molecules)
                            if (find(m.atoms.begin(), m.atoms.end(), b->a) == m.atoms.end()) { // contains b, but not a
                                if (periodic) {//unwrap system by moving a
                                    b->a->x = b->b->x + b->ba.getX();
                                    b->a->y = b->b->y + b->ba.getY();
                                    b->a->z = b->b->z + b->ba.getZ();
                                }
                                m.atoms.push_back(b->a); //add a to molecule
                                m.atomIndeces.push_back(b->a_indx);
                                b = bonds.erase(b);
                                --b;
                            }
                        }

                    }
                } while (oldMsize != m.atoms.size()); //keep looping until nothing is added to m

                //we are done with this molecule, save it to mols
                mols.push_back(m);
                CTP_LOG(ctp::logDEBUG, *_log) << "BreakIntoMolecules(): put " << m.atoms.size()
                        << " atoms into molecule " << mols.size() - 1 << ";\t" << bonds.size() << " bonds left." << flush;
                
                if(m.atoms.size() != 3)
                {
                    cout<<"Unexpected number of atoms in molecule"<< mols.size() <<endl;
                }
            }

            if (periodic)
                CTP_LOG(ctp::logDEBUG, *_log) << " BreakIntoMolecules(): Molecules have been unwrapped." << flush;
            
            return (mols);
        }

        ub::matrix<double> Bulkesp::BuildDenMat(Orbitals &_orb, std::string _state, std::string _spin, int _state_no) {
            ub::matrix<double> DMAT_tot;
            bool _do_transition = false;

            if (_state != "ground")
                CTP_LOG(ctp::logDEBUG, *_log) << " Bulkesp is not tested for any state other than the ground state. "
                << "It will likely not work correctly, as BSE Singlet and Triplet Coefficients are simply copied from the global orbital object.\n" << flush;

            if (_state == "transition") {
                _do_transition = true;
                if (_spin == "singlet") {
                    DMAT_tot = _orb.TransitionDensityMatrix(_orb.MOCoefficients(), _orb.BSESingletCoefficients(), _state_no - 1);
                } else if (_spin == "triplet") {
                    DMAT_tot = _orb.TransitionDensityMatrix(_orb.MOCoefficients(), _orb.BSETripletCoefficients(), _state_no - 1);
                } else throw std::runtime_error("Spin entry not recognized");
            } else if (_state == "ground" || _state == "excited") {
                DMAT_tot = _orb.DensityMatrixGroundState(_orb.MOCoefficients());
                if (_state_no > 0 && _state == "excited") {
                    std::vector<ub::matrix<double> > DMAT;
                    if (_spin == "singlet") {
                        DMAT = _orb.DensityMatrixExcitedState(_orb.MOCoefficients(), _orb.BSESingletCoefficients(), _state_no - 1);
                    } else if (_spin == "triplet") {
                        DMAT = _orb.DensityMatrixExcitedState(_orb.MOCoefficients(), _orb.BSETripletCoefficients(), _state_no - 1);
                    } else throw std::runtime_error("Spin entry not recognized");
                    DMAT_tot = DMAT_tot - DMAT[0] + DMAT[1];
                }
                // Ground state + hole_contribution + electron contribution
            } else throw std::runtime_error("State entry not recognized");


            return (DMAT_tot);
        }

        ub::matrix<double> Bulkesp::BuildOverlapMat(Orbitals &_molOrb, Orbitals &_globalOrb) {
            ub::matrix<double> ov;
            return (ov);
        }

        /**
         * 
         * @param _global_atomlist
         * @param _local_atomlist
         * @param _local_atomIndeces Indeces of atoms in current molecule
         * @param _global_dmat
         * @param _global_basis
         * @param bs
         * @param gridsize
         * @param _grid
         * @return 
         */
        ub::vector<double> Bulkesp::ComputeESP(std::vector< ctp::QMAtom* > & _global_atomlist,
                std::vector< ctp::QMAtom* > & _local_atomlist, std::vector<unsigned> _local_atomIndeces,
                ub::matrix<double> &_global_dmat, AOBasis &_global_basis, BasisSet &bs, string gridsize, Grid &_grid, double &netcharge) {


            ub::vector<double> _ESPatGrid = ub::zero_vector<double>(_grid.getsize());

            NumericalIntegrationPeriodic numway;

            //numway.GridSetup(gridsize,&bs,_global_atomlist);
            numway.setBox(boxLen*tools::conv::ang2bohr);
            numway.SetRelevantAtomIds(_local_atomIndeces);
            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Entering Bulkesp::ComputeESP()"<< flush;
            numway.GridSetup(gridsize, &bs, _local_atomlist, &_global_basis);
            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculate Potentials at Numerical Grid with gridsize " << gridsize << flush;
            //As long as basis functions are well supported and molecules are smaller than 0.5*boxLen along any axis, then
            //density integration should be accurate enough without making it explicitly periodic
            //double N = numway.IntegrateDensity_Molecule(_global_dmat, &_global_basis, _local_atomIndeces);
            double N = numway.IntegrateDensity(_global_dmat);
            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculated Potentials at Numerical Grid, Number of electrons is " << N << flush;

            _do_round = false; //do not round net charge, we expect total charge to be non-int, as molecules can transfer some between themselves
            netcharge = getNetcharge(_local_atomlist, N, false); //don't round

            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating ESP at CHELPG grid points" << flush;
            //boost::progress_display show_progress( _grid.getsize() );





            double numK;

            //#define MADELUNG_TEST
#ifdef MADELUNG_TEST


            double exactMadelung = 1.74756459463318;

#define GROLIKE
#ifdef GROLIKE
            ofstream myfile("Energy_kmax8_gro.dat");
#else
            ofstream myfile("Energy_kmax8.dat");
#endif

            int natomsonside = 2;
            double numK = 5;
            double a = 5.6402 * tools::conv::ang2bohr;
            cout << "a = " << a << endl;
            cout << "nearest neighbour distance = " << a / 2 << endl;
            double BL[3];
            BL[0] = a;
            BL[1] = a;
            BL[2] = a;
            std::vector< ctp::QMAtom* > fake_atom_list;
            fake_atom_list.resize(0);
            std::vector< vec > points;
            points.push_back(vec(0, 0, 0) * tools::conv::bohr2nm); //0,0,0 in a
            points.push_back(vec(0.0, 0.0, a / 2) * tools::conv::bohr2nm); //0,0,0.5 in a
            Grid eval_grid(points); //in nm
            _ESPatGrid = ub::zero_vector<double>(eval_grid.getsize());


            for (double alpha = 0.01; alpha < 4; alpha += 0.01) {
                numway.FillMadelungGrid(BL, natomsonside);
#ifdef GROLIKE
                numway.PrepKspaceDensity_gromacs_like(alpha, fake_atom_list, _ECP, eval_grid, numK);
                numway.IntegratePotential_w_PBC_gromacs_like(eval_grid, _ESPatGrid);
#else                
                numway.PrepKspaceDensity(BL, alpha, fake_atom_list, _ECP, numK);
                for (int i = 0; i < eval_grid.getsize(); i++) {
                    _ESPatGrid(i) = numway.IntegratePotential_w_PBC(eval_grid.getGrid()[i] * tools::conv::nm2bohr, BL);
                }
                numway.FreeKspace();
#endif

                cout << "Madelung constant is: " << _ESPatGrid(0)*(a / 2) << "\n";
                myfile << natomsonside << " \t" << numway.numK[0] << " \t" << numway.alpha << " \t"
                        //<<std::abs(_ESPatGrid(0)*(a/natomsonside)) - exactMadelung<<" \t"
                        << _ESPatGrid(0)*(a / 2) << " \t"
                        << numway.E_rspace * (a / 2) << " \t" << numway.E_kspace * (a / 2) << " \t" << numway.E_erfc * (a / 2)
                        << endl;

            }
            exit(0);

#endif        











            if (periodic) {
                CTP_LOG(ctp::logDEBUG, *_log)  << "Bulkesp::ComputeESP(): " << ctp::TimeStamp() << " periodicity is on, including long range contributions." << flush;
                tools::vec BL = boxLen * tools::conv::ang2bohr; //bohr

                numK = 32;
                float alpha=1;

                numway.PrepKspaceDensity_gromacs_like(alpha, _local_atomlist, _ECP, _grid, numK);
                numway.IntegratePotential_w_PBC_gromacs_like(_grid, _ESPatGrid);
                
                
                
//                cout<<endl;
//                for(int i=0; i<10; i++){
//                    alpha=1.0 + (i*1.0);
//                    numway.PrepKspaceDensity_gromacs_like(alpha, _local_atomlist, _ECP, _grid, numK);
//                    numway.IntegratePotential_w_PBC_gromacs_like(_grid, _ESPatGrid);
//                    cout<<"numK= "<< numK<< "\talpha= "<< alpha<<"\t"<< _ESPatGrid(1)<<"\t"<< _ESPatGrid(50)<<"\t"<< _ESPatGrid(100)<< endl<<flush;
//                }
//                cout<<endl;
//                numK = 32;
//                for(int i=0; i<10; i++){
//                    alpha=1.0 + (i*1.0);
//                    numway.PrepKspaceDensity_gromacs_like(alpha, _local_atomlist, _ECP, _grid, numK);
//                    numway.IntegratePotential_w_PBC_gromacs_like(_grid, _ESPatGrid);
//                    cout<<"numK= "<< numK<< "\talpha= "<< alpha<<"\t"<< _ESPatGrid(1)<<"\t"<< _ESPatGrid(50)<<"\t"<< _ESPatGrid(100)<< endl<<flush;
//                }
//                cout<<endl;
//                numK = 64;
//                for(int i=0; i<10; i++){
//                    alpha=1.0 + (i*1.0);
//                    numway.PrepKspaceDensity_gromacs_like(alpha, _local_atomlist, _ECP, _grid, numK);
//                    numway.IntegratePotential_w_PBC_gromacs_like(_grid, _ESPatGrid);
//                    cout<<"numK= "<< numK<< "\talpha= "<< alpha<<"\t"<< _ESPatGrid(1)<<"\t"<< _ESPatGrid(50)<<"\t"<< _ESPatGrid(100)<< endl<<flush;
//                }
//                cout<<endl;
//                numK = 128;
//                for(int i=0; i<10; i++){
//                    alpha=1.0 + (i*1.0);
//                    numway.PrepKspaceDensity_gromacs_like(alpha, _local_atomlist, _ECP, _grid, numK);
//                    numway.IntegratePotential_w_PBC_gromacs_like(_grid, _ESPatGrid);
//                    cout<<"numK= "<< numK<< "\talpha= "<< alpha<<"\t"<< _ESPatGrid(1)<<"\t"<< _ESPatGrid(50)<<"\t"<< _ESPatGrid(100)<< endl<<flush;
//                }
//                exit(0);
                
                
                

                /*
                numway.PrepKspaceDensity(BL, 0.5, _local_atomlist, _ECP);
                CTP_LOG(ctp::logDEBUG, *_log) << " Bulkesp::ComputeESP(): Found density in Fourier space"<< endl;
                #pragma omp parallel for
                for ( int i = 0 ; i < _grid.getsize(); i++){
                    _ESPatGrid(i)=numway.IntegratePotential_w_PBC(_grid.getGrid()[i]*tools::conv::nm2bohr, BL);
                    //++show_progress;
                }
                 */


                CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Electron and Nuclear contributions calculated" << flush;
                //exit(0);

                //calculate and record the molecular dipole moments
                ub::vector<double> dipPos(3);
                ctp::QMAtom* atom = *(_local_atomlist.begin());
                dipPos(0) = atom->x * tools::conv::ang2bohr;
                dipPos(1) = atom->y * tools::conv::ang2bohr;
                dipPos(2) = atom->z * tools::conv::ang2bohr;
                double dipole = numway.CalcDipole_w_PBC(dipPos); //in bohr * e
                CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Molecular dipole: " << dipole / 0.393430307 << "Debye" << flush;
                *dipolesLog << dipole / 0.393430307 << endl;

                //numway.FreeKspace();

            } else {
                CTP_LOG(ctp::logDEBUG, *_log) << " Bulkesp::ComputeESP(): periodicity is off, no long range contributions." << endl;

                //numway.SetGridToCharges(_local_atomlist);
#pragma omp parallel for
                for (int i = 0; i < _grid.getsize(); i++) {
                    _ESPatGrid(i) = numway.IntegratePotential(_grid.getGrid()[i] * tools::conv::nm2bohr);
                    //++show_progress;
                }

                CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Electron contribution calculated" << flush;

                // Calculating nuclear potential at gridpoints
                ub::vector<double> _NucPatGrid = EvalNuclearPotential(_local_atomlist, _grid);
                _ESPatGrid += _NucPatGrid;
                CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Nuclear contribution calculated" << flush;
            }

            return (_ESPatGrid);
        }

        void Bulkesp::FillElement2NBF(std::vector< ctp::QMAtom* >& _atomlist, BasisSet &bs) {

            //Find list of elements
            _elementNames.clear();
            for (std::vector< ctp::QMAtom* >::iterator ait = _atomlist.begin(); ait < _atomlist.end(); ait++) {

                std::string element_name = (*ait)->type;
                list<std::string>::iterator ite;
                ite = find(_elementNames.begin(), _elementNames.end(), element_name);
                if (ite == _elementNames.end()) { //this is the first atom of this element encountered
                    _elementNames.push_back(element_name);
                }
            }

            //Find NBF for each element
            _element2NBF.clear();
            for (list<std::string>::iterator eit = _elementNames.begin(); eit != _elementNames.end(); ++eit) {
                vector<Shell*> shells = bs.getElement(*eit)->getShells();
                int nfunc = 0;
                //loop over shells
                for (vector< Shell* >::iterator sit = shells.begin(); sit != shells.end(); ++sit) {
                    nfunc += 1 + 2 * ((*sit)->getLmax());
                }
                //store number of functions
                _element2NBF[(*eit)] = nfunc;
            }

        }



        //Find indeces of MO coefficients of the first shell of each atom,
        //assuming that the order of _atomlist is the same as that of the _MO_Coefficients.

        std::map<ctp::QMAtom*, int> Bulkesp::MapAtom2MOCoefIndex(std::vector< ctp::QMAtom* >& _atomlist) {
            std::map<ctp::QMAtom*, int> ret;
            int count = 0;
            //loop over all atoms in system
            for (vector<ctp::QMAtom*>::iterator i = _atomlist.begin(); i != _atomlist.end(); ++i) {
                ret[(*i)] = count;
                count += _element2NBF[(*i)->type];
            }
            return (ret);
        }

        void Bulkesp::Evaluate(std::vector< ctp::QMAtom* >& _atomlist, ub::matrix<double> &_global_dmat, Orbitals& _globalOrb,
                ub::matrix<double> _global_MO_Coeffs, AOBasis &_basis, BasisSet &bs,
                string gridsize, double maxBondScale, std::string _state,
                std::string _spin, int _state_no) {

            //find the individual molecules
            std::vector<Bulkesp::Molecule> mols = BreakIntoMolecules(_atomlist, maxBondScale);

            //open/create dipolesLog
            dipolesLog->open("dipolesLog.dat", ios_base::trunc);
            
            double system_netcharge=0;

            //loop over molecules
            CTP_LOG(ctp::logDEBUG, *_log) << "Bulkesp::Evaluate(): found " << mols.size() << " molecules." << endl << flush;
            for (std::vector<Bulkesp::Molecule>::iterator m = mols.begin(); m != mols.end(); ++m) {

                CTP_LOG(ctp::logDEBUG, *_log) << "Bulkesp::Evaluate(): " << ctp::TimeStamp() << " processing molecule " << m - mols.begin() << flush;

                //verify atomic coordinates and units
                for (std::vector<ctp::QMAtom*>::iterator a = m->atoms.begin(); a != m->atoms.end(); ++a) {
                    ctp::QMAtom* ap = *a;
                    CTP_LOG(ctp::logDEBUG, *_log) << ap->type << '\t' << ap->x << '\t' << ap->y << '\t' << ap->z << "A" << flush;
                }
                CTP_LOG(ctp::logDEBUG, *_log) << "box: " << boxLen[0] << '\t' << boxLen[1] << '\t'<< boxLen[2] << "A" << flush;

                //set up grid
                Grid _grid(true, false, false); //create polarsites, so we can output grid to .cube file
                //            _grid.setAtomlist(&m->atoms);
                //            _grid.setupCHELPgrid();
                _grid.setCubegrid(true);
                if (periodic) {
                    //_grid.setPadding(0.0);
                    _grid.setPeriodicity(boxLen);
                }
                
//                //let's just consider the first Oxygen atom as a molecule
//                m->atoms= vector<ctp::QMAtom*>(m->atoms.begin()+0, m->atoms.begin()+2);
//                m->atomIndeces= vector<unsigned>(m->atomIndeces.begin()+0, m->atomIndeces.begin()+2);

                //test: set inner cutoff to 0 and calculate all potentials near nuclei
//                _grid.setCutoffs(20, 0.001); //between 1.5 and 3 A, as that is the region where water-water interactions take place
//                _grid.setSpacing(0.3); //defaults to 0.3 A
                _grid.setAtomlist(&m->atoms);
//                _grid.setAtomlist(&_atomlist);
                _grid.setupgrid();
                _grid.setupCHELPgrid();
                CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Done setting up CHELPG grid with " << _grid.getsize() << " points " << flush;

                //calculate the ESP
                //ub::vector<double> ESP=ComputeESP(m->atoms, _m_dmat, _m_ovmat, _m_basis, bs, gridsize, _grid);
                double netcharge = 0.0;
                ub::vector<double> ESP = ComputeESP(_atomlist, m->atoms, m->atomIndeces,
                        _global_dmat, _basis, bs, gridsize, _grid, netcharge);
                
//                std::vector< unsigned > allIndeces;
//                for (vector<ctp::QMAtom*>::iterator i = _atomlist.begin(); i != _atomlist.end(); ++i) {
//                    allIndeces.push_back(i - _atomlist.begin());
//                }
//                ub::vector<double> ESP = ComputeESP(_atomlist, _atomlist, allIndeces,
//                        _global_dmat, _basis, bs, gridsize, _grid, netcharge);

                //store the potential in apolarsites
                for (int i = 0; i < _grid.getsize(); i++) {
                    //                ub::vector<double> point = _grid.getGrid()[i];
                    //                APolarSite* site = _grid.Sites()[i];
                    //                site->setPhi(ESP(i), 0.0);
                    _grid.Sites()[i]->setPhi(ESP(i), 0.0);
                }
                
                system_netcharge+=netcharge;


                //and save it to a .cube file
                std::ostringstream fn;
                if (periodic && _grid.getCubegrid()) {
                    fn.str(std::string());
                    fn << "BulkEsp_" << m - mols.begin() << ".cube";
                    _grid.printgridtoCubefile(fn.str());
                }

                //output
                //CHELPG grids aren't periodic and equally spaced,
                //so can't output to .cube format.
                //Create own format.
                //note: just like printgridtoCubefile, this prints potential from apolar sites
                fn.clear();
                fn.str("");
                fn << "BulkEsp_" << m - mols.begin() << ".grid";
                _grid.writeIrregularGrid(fn.str(), m->atoms, _ECP);


                //TODO: fit charges
                std::vector< tools::vec > _fitcenters;

                for (unsigned j = 0; j < m->atoms.size(); j++) {
                    tools::vec _pos = m->atoms[j]->getPos() * tools::conv::ang2nm;
                    _fitcenters.push_back(_pos);
                }

                std::vector<double> _charges = FitPartialCharges(_fitcenters, _grid, ESP, netcharge);


                //Write charges to qmatoms
                for (unsigned _i = 0; _i < m->atoms.size(); _i++) {
                    m->atoms[_i]->charge = _charges[_i];
                }


                CTP_LOG(ctp::logDEBUG, *_log) << "Bulkesp::Evaluate(): " << ctp::TimeStamp() << " done with molecule " << m - mols.begin() << endl << flush;

            }
            CTP_LOG(ctp::logDEBUG, *_log) << "Bulkesp::Evaluate(): " << ctp::TimeStamp() << " All molecules processed." << endl << flush;
            dipolesLog->close();

            CTP_LOG(ctp::logDEBUG, *_log) << "Bulkesp::Evaluate(): " << ctp::TimeStamp() << " Net charge of the whole system is "<< system_netcharge << endl << flush;
            if(fabs(system_netcharge)>0.05){
                if(fmod(system_netcharge,1.0)<0.05){
                    CTP_LOG(ctp::logWARNING, *_log) << "Bulkesp::Evaluate(): " << ctp::TimeStamp() << " System is ionized. Is this intended?"<< endl << flush;
                }
                else{
                    CTP_LOG(ctp::logWARNING, *_log) << "Bulkesp::Evaluate(): " << ctp::TimeStamp() << " System has significant net partial charge. Projection was likely incomplete. "<< endl << flush;
                }
            }
        }

        
        
        
        
        /*
         * Modified version of espfit::FitPartialCharges() that allows for periodic boundary conditions
         * 
         */
        std::vector<double> Bulkesp::FitPartialCharges(std::vector< tools::vec >& _fitcenters, Grid& _grid, ub::vector<double>& _potential, double& _netcharge) {
            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Setting up Matrices for fitting of size " << _fitcenters.size() + 1 << " x " << _fitcenters.size() + 1 << flush;

            const std::vector< tools::vec >& _gridpoints=_grid.getGrid();   

            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Using " << _fitcenters.size() << " Fittingcenters and " << _gridpoints.size() << " Gridpoints." << flush;

            ub::matrix<double> _Amat = ub::zero_matrix<double>(_fitcenters.size() + 1, _fitcenters.size() + 1);
            ub::matrix<double> _Bvec = ub::zero_matrix<double>(_fitcenters.size() + 1, 1);
            //boost::progress_display show_progress( _fitcenters.size() );
            // setting up _Amat
            
            tools::vec L; //periodic box in nm
            if(periodic){
                L = boxLen * tools::conv::ang2nm;
            }
#pragma omp parallel for
            for (unsigned _i = 0; _i < _Amat.size1() - 1; _i++) {
                for (unsigned _j = _i; _j < _Amat.size2() - 1; _j++) {
                    for (unsigned _k = 0; _k < _gridpoints.size(); _k++) {
                        tools::vec dif_i = _fitcenters[_i] - _gridpoints[_k];
                        tools::vec dif_j = _fitcenters[_j] - _gridpoints[_k];
                        if(periodic){ //adjust for periodic boundary conditions
                            dif_i=WrapPoint(dif_i, L);
                            dif_j=WrapPoint(dif_j, L);
                        }
                        double dist_i = tools::abs(dif_i) * tools::conv::nm2bohr;
                        double dist_j = tools::abs(dif_j) * tools::conv::nm2bohr;

                        _Amat(_i, _j) += 1.0 / dist_i / dist_j;
                    }
                    _Amat(_j, _i) = _Amat(_i, _j);
                }
            }

            for (unsigned _i = 0; _i < _Amat.size1(); _i++) {
                _Amat(_i, _Amat.size1() - 1) = 1.0;
                _Amat(_Amat.size1() - 1, _i) = 1.0;
            }
            _Amat(_Amat.size1() - 1, _Amat.size1() - 1) = 0.0;

            // setting up Bvec
#pragma omp parallel for
            for (unsigned _i = 0; _i < _Bvec.size1() - 1; _i++) {
                for (unsigned _k = 0; _k < _gridpoints.size(); _k++) {
                    double dist_i = tools::abs(_fitcenters[_i] - _gridpoints[_k]) * tools::conv::nm2bohr;
                    _Bvec(_i, 0) += _potential(_k) / dist_i;
                }
            }

            _Bvec(_Bvec.size1() - 1, 0) = _netcharge; //netcharge!!!!
            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Inverting Matrices " << flush;
            // invert _Amat
            ub::matrix<double> _Amat_inverse = ub::zero_matrix<double>(_fitcenters.size() + 1, _fitcenters.size() + 1);



            if (_do_svd) {
                int notfittedatoms = linalg_invert_svd(_Amat, _Amat_inverse, _conditionnumber);
                CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << "SVD Done. " << notfittedatoms << " could not be fitted and are set to zero." << flush;
            } else {
                linalg_invert(_Amat, _Amat_inverse);
            }

            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Inverting Matrices done." << flush;
            //_Amat.resize(0,0);



            ub::matrix<double> _charges = ub::prod(_Amat_inverse, _Bvec);

            std::vector<double> _result;
            for (unsigned _i = 0; _i < _charges.size1(); _i++) {
                _result.push_back(_charges(_i, 0));
            }

            double _sumcrg = 0.0;
            for (unsigned _i = 0; _i < _fitcenters.size(); _i++) {

                //CTP_LOG(logDEBUG, *_log) << " Center " << _i << " FitCharge: " << _result[_i] << flush;
                _sumcrg += _result[_i];
            }

            CTP_LOG(ctp::logDEBUG, *_log) << " Sum of fitted charges: " << _sumcrg << flush;

            // get RMSE
            double _rmse = 0.0;
            double _totalPotSq = 0.0;
            for (unsigned _k = 0; _k < _gridpoints.size(); _k++) {
                double temp = 0.0;
                for (unsigned _i = 0; _i < _fitcenters.size(); _i++) {
                    double dist = tools::abs(_gridpoints[_k] - _fitcenters[_i]) * tools::conv::nm2bohr;
                    temp += _result[_i] / dist;
                }
                _rmse += (_potential(_k) - temp)*(_potential(_k) - temp);
                _totalPotSq += _potential(_k) * _potential(_k);
            }
            CTP_LOG(ctp::logDEBUG, *_log) << " RMSE of fit:  " << sqrt(_rmse / _gridpoints.size()) << flush;
            CTP_LOG(ctp::logDEBUG, *_log) << " RRMSE of fit: " << sqrt(_rmse / _totalPotSq) << flush;

            return _result;
        }
        
        
        
        double Bulkesp::getNetcharge(std::vector< ctp::QMAtom* >& _atoms, double N, bool doround) {
            double netcharge = 0.0;
            if (std::abs(N) < 0.05) {
                //CTP_LOG(ctp::logDEBUG, *_log) << "Number of Electrons is "<<N<< " transitiondensity is used for fit"  << flush;
                _do_Transition = true;
            } else {
                double Znuc_ECP = 0.0;
                double Znuc = 0.0;
                for (unsigned j = 0; j < _atoms.size(); j++) {
                    Znuc_ECP += _elements.getNucCrgECP(_atoms[j]->type);
                    Znuc += _elements.getNucCrg(_atoms[j]->type);
                }

                if (_ECP) {
                    if (std::abs(Znuc_ECP - N) < 4) {
                        CTP_LOG(ctp::logDEBUG, *_log) << "Number of Electrons minus ECP_Nucleus charge is " << Znuc_ECP - N << " you use ECPs, sounds good" << flush;
                    } else if (std::abs(Znuc - N) < 4) {
                        CTP_LOG(ctp::logDEBUG, *_log) << "Number of Electrons minus real Nucleus charge is " << Znuc - N << " you are sure you want ECPs?" << flush;
                    } else {
                        CTP_LOG(ctp::logDEBUG, *_log) << "Warning: Your molecule is highly ionized and you want ECPs, sounds interesting" << flush;
                    }
                    netcharge = Znuc_ECP - N;
                } else {
                    if (std::abs(Znuc - N) < 4) {
                        CTP_LOG(ctp::logDEBUG, *_log) << "Number of Electrons minus Nucleus charge is " << Znuc - N << " you probably do not use ECPs, if you do use ECPs please use the option. Otherwise you are fine" << flush;
                    } else if (std::abs(Znuc_ECP - N) < 4) {
                        CTP_LOG(ctp::logDEBUG, *_log) << "Number of Electrons minus ECP_Nucleus charge is " << Znuc_ECP - N << " you probably use ECPs, if you do use ECPs please use the option to switch on" << flush;
                    } else {
                        CTP_LOG(ctp::logDEBUG, *_log) << "Warning: Your molecule is highly ionized and you use real core potentials, sounds interesting" << flush;
                    }
                }
                _do_Transition = false;
            }
            if(doround){
                CTP_LOG(ctp::logDEBUG, *_log) << "Rounding Netcharge from " << netcharge<< flush;
            }
            CTP_LOG(ctp::logDEBUG, *_log) << "Netcharge constrained to " << netcharge << flush;

            return netcharge;
        }    


    }//xtp
}//votca
