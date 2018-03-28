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

#include "cpmd.h"
#include <votca/ctp/segment.h>

#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/constants.h>
#include <stdio.h>
#include <iomanip>
#include <sys/stat.h>
#include <vector>

#ifdef DEBUG
#include <votca/xtp/aomatrix.h>
#endif



namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        void Cpmd::Initialize(Property *options) {

            // NWChem file names
            string fileName = "system";

            _xyz_file_name = fileName + ".xyz";
            _input_file_name = fileName + ".inp";
            _log_file_name = fileName + ".log";
            _orb_file_name = "WFNCOEF" ;

            string key = "package";
            string _name = options->get(key+".name").as<string> ();

            if ( _name != "cpmd" ) {
                cerr << "Tried to use " << _name << " package. ";
                throw std::runtime_error( "Wrong options file");
            }

            _executable =       options->get(key + ".executable").as<string> ();
            _charge =           options->get(key + ".charge").as<int> ();
            _spin =             options->get(key + ".spin").as<int> ();
            _options =          options->get(key + ".options").as<string> ();
            _cleanup =          options->get(key + ".cleanup").as<string> ();
            _threads =          options->get(key + ".threads").as<int> ();
            if(_threads!=1){
                CTP_LOG(ctp::logWARNING, *_pLog) << "CPMD: CPMD only supports "
                        "MPI paralellization.\n"
                        "If you want a parallel run of CPMD, set "
                        " the executable name to: mpirun -np <# of threads> cpmd-mpi.x"
                        "Setting number of threads to 1 and continuing." << flush;
                _threads=1;
            }



            //restart?
            if (options->exists(key + ".restart")) {
                _rsrt=true;
                _rsrt_kwds = options->get(key + ".restart").as<std::string> ();
            }
            else _rsrt=false;

            //optimize wavefunction?
            if (options->exists(key + ".optimizewf")) {
                _optWF=true;
                _convCutoff = options->get(key + ".optimizewf").as<double> ();
            }
            else _optWF=false;
            
            //use VDW correction?
            if (options->exists(key + ".useGrimmeVDW")) {
                _useGrimmeVDW=true;
            }
            else _useGrimmeVDW=false;
            
            //custom CPMD controlls
            if (options->exists(key + ".customCPMDcontrolls")) {
                _custom_CPMD_controlls=options->get(key + ".customCPMDcontrolls").as<std::string> ();
            }
            else _custom_CPMD_controlls="";


            //functional and pseudopotentials
            if (options->exists(key + ".functional")) {
                if (!(options->exists(key + ".functional"))) throw std::runtime_error("Functional name missing");
                _functional = options->get(key + ".functional.name").as<std::string> ();
                
                list<Property*>::iterator pit;
                std::list<Property *> props = options->Select(key + ".functional.pseudopotentials.*");
                for (pit = props.begin(); pit != props.end(); ++pit)
                {
                    Property* p= *pit;
                    _ppFileNames[p->name()]=p->as<std::string> ();
                }
                props = options->Select(key + ".functional.l.*");
                for (pit = props.begin(); pit != props.end(); ++pit)
                {
                    Property* p= *pit;
                    _ppLData[p->name()]=p->as<std::string> ();
                }
                
            }
            else throw std::runtime_error("No functional and pseudopotentials specified");
            
            //symmetry
            _symmetry=0;
            if (options->exists(key + ".symmetry")) {
                _symmetry = options->get(key + ".symmetry").as<int> ();
            }
            else
                CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: no symmetry provided, assuming simple cubic." << flush;
            
            //cell
            _cell="20.0   1.0   1.0  0.0  0.0  0.0";
            if (options->exists(key + ".cell")) {
                _cell = options->get(key + ".cell").as<std::string> ();
            }
            else
                CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: no cell provided, assuming cube with side length of 20 Bohr." << flush;
            
            //plane wave cutoff
            _pwCutoff=80.0;
            if (options->exists(key + ".pwcutoff")) {
                _pwCutoff = options->get(key + ".pwcutoff").as<double> ();
            }
            else
                CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: no plane wave cutoff provided, assuming "<< _pwCutoff <<" Ry." << flush;

            //output electrostatic potential?
            _elpot=false;
            if (options->exists(key + ".elpot")) {
                _elpot=options->get(key + ".elpot").as<bool> ();
            }

            //project wavefunction onto atomic orbitals?
            _projectWF=false;
            if (options->exists(key + ".projectwf")) {
                _projectWF=options->get(key + ".projectwf").as<bool> ();
            }
            
            //basis set file name
            if (options->exists(key + ".basisset")) {
            _basisset_name=options->get(key + ".basisset").as<string> ();
            }

            //do population analysis? required to have access to WF coefs in atom-centric basis
            //requires _projectWF
            _popAnalysis=false;
            if (options->exists(key + ".popanalysis")) {
                _popAnalysis=options->get(key + ".popanalysis").as<bool> ();
            }

            //get density and overlap matrices during post processing?
            //requires _popAnalysis and _projectWF
            _getMat=false;
            if (options->exists(key + ".getmatrices")) {
                _getMat=options->get(key + ".getmatrices").as<bool> ();
            }

            if(_getMat) _popAnalysis=true;
            if(_popAnalysis) _projectWF=true;

            if(_projectWF && _optWF){
                cerr << "Warning: Wavefunction optimization and projection onto atom-centric orbitals can not be done together.\nCPMD would crash.\n";
                cerr << "Splitting into optimization on first run and projection/population analysis on second.\n";
                CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: Splitting run into two steps." << flush;
            }

        }

        bool Cpmd::WriteInputFile(std::vector<ctp::Segment* > segments, Orbitals* orbitals_guess, std::vector<ctp::PolarSeg*> PolarSegments) {

            std::vector< ctp::Atom* > _atoms;
            std::vector< ctp::Atom* > ::iterator ait;
            std::vector< ctp::Segment* >::iterator sit;
            std::string temp_suffix = "/id";

            ofstream _com_file;

            std::string _com_file_name_full = _run_dir + "/" + _input_file_name;
            //if need to run CPMD twice, create new names for input and log files of the first run
            if(_projectWF && _optWF){
                //input
                size_t pointpos = _input_file_name.find('.');
                if (pointpos!=std::string::npos){
                    _wfOpt_input_file_name=_input_file_name.substr(0,pointpos) + "_wfOpt" + _input_file_name.substr(pointpos);
                }else{
                _wfOpt_input_file_name=_input_file_name + "_wfOpt";
                }
                _com_file_name_full= _run_dir + "/" + _wfOpt_input_file_name;

                //log
                pointpos = _log_file_name.find('.');
                if (pointpos!=std::string::npos){
                    _wfOpt_log_file_name=_log_file_name.substr(0,pointpos) + "_wfOpt" + _log_file_name.substr(pointpos);
                }else{
                _wfOpt_log_file_name=_log_file_name + "_wfOpt";
                }
            }

            _com_file.open(_com_file_name_full.c_str());

            // header
            _com_file << "&INFO\nGenerated by VOTCA\n&END" << endl;


            //control
            _com_file << "\n&CPMD\n";
            if(_rsrt) _com_file << "  RESTART " << _rsrt_kwds << endl;  //restart
            if(_optWF){                                                 //optimize WF
                _com_file << "  OPTIMIZE WAVEFUNCTION" << endl;
                _com_file << "  CONVERGENCE ORBITALS" << endl;
                _com_file << "  " << FortranFormat(_convCutoff) << endl;
                _com_file << "  PCG MINIMIZE" << endl;              //use the more stable optimizer
                _com_file << "  TIMESTEP" << endl;
                _com_file << "   20" << endl;
                _com_file << "  STORE WAVEFUNCTIONS" <<endl;
                _com_file << "   20 SC=20" <<endl;


            }
            if(_elpot){                                                 //output electrostatic potential
                _com_file << "  ELECTROSTATIC POTENTIAL" << endl;
                _com_file << "  RHOOUT" << endl;
            }
            if(_projectWF && !_optWF){
                _com_file << "  PROPERTIES" << endl;
            }
            
            _com_file << endl <<"  " << _custom_CPMD_controlls << endl;        //custom CPMD controls from the .xml file
            if(_useGrimmeVDW){                                          //VDW
                _com_file << "  VDW CORRECTION" << endl;
            }
            _com_file << "&END" << endl;

            //functional
            _com_file << "\n&DFT\n";
            _com_file << "  FUNCTIONAL " << _functional << endl;
            _com_file << "  GC-CUTOFF" << endl;
            _com_file << "   1.0d-06" << endl;
            _com_file << "&END" << endl;
            
            //VDW
            if(_useGrimmeVDW){  
                _com_file << "\n&VDW" << endl;
                _com_file << " VDW CORRECTION" << endl;
                _com_file << " ALL GRIMME" << endl;
                _com_file << " VDW-CELL" << endl;
                _com_file << " 1 1 1" << endl;
                _com_file << "&END" << endl;
            }
            
            //cell
            _com_file << "\n&SYSTEM\n";
            _com_file << "  SYMMETRY" << endl;
            _com_file << "   " << _symmetry <<endl;
            _com_file << "  CELL" << endl;
            _com_file << "   " << _cell <<endl;
            _com_file << "  CUTOFF" << endl;
            _com_file << "   " << FortranFormat(_pwCutoff) <<endl;
            _com_file << "&END" << endl;
            
            //properties
            if(_projectWF){
                _com_file << "\n&PROP\n";
                _com_file << "  PROJECT WAVEFUNCTION" << endl;
                if(_popAnalysis){
                    _com_file << "  CHARGES" << endl;
                    _com_file << "  POPULATION ANALYSIS MULLIKEN" << endl;
                }
                _com_file << "&END" << endl;
            }
            
            //basis
            if(_projectWF && !_optWF){ //prevent initializing _optWF runs with this basis
                _com_file << "\n&BASIS\n";
                WriteBasisSet(segments, _com_file);
                _com_file << "&END" << endl;
            }
            
            //atoms
            _com_file << "\n&ATOMS\n";
                //find how many atoms of each element there are
            //list<std::string> elements;
            int numAtoms=0;
            if(_elements.empty()){ //we have not tabulated # atoms of each element yet
                //(this is the first time this function is called)
                for (sit = segments.begin(); sit != segments.end(); ++sit) {

                    _atoms = (*sit)-> Atoms();

                    for (ait = _atoms.begin(); ait < _atoms.end(); ait++) {

                        std::string element_name = (*ait)->getElement();
                        list<std::string>::iterator ite;
                        ite = find(_elements.begin(), _elements.end(), element_name);
                        if (ite == _elements.end()) {            //this is the first atom of this element encountered
                            _elements.push_back(element_name);
                            _nAtomsOfElement[element_name]=1;
                        }
                        else
                        {
                            _nAtomsOfElement[element_name]++;
                        }
                        numAtoms++;
                    }
                }
            }
                //now loop over elements and store all atoms of that element
            bool atomOrderMapSet = (VOTCA2CPMD_map!=NULL && CPMD2VOTCA_map!=NULL);
            int Vind=0, Cind=0; //atom indexes in VOTCA and CPMD
            if(!atomOrderMapSet){
                VOTCA2CPMD_map=new int[numAtoms];
                CPMD2VOTCA_map=new int[numAtoms];
            }
            list<std::string>::iterator ite;
            for (ite = _elements.begin(); ite != _elements.end(); ite++) {
                if(_ppFileNames.find(*ite)==_ppFileNames.end()) {
                    cerr << "Error: Element "<<(*ite)<<" has no pseudopotential specified in CPMD options file.\n" << flush;
                    throw std::runtime_error("Encountered element with no pseudopotential.\n");
                }
                else{
                    _com_file << "*" << _ppFileNames[(*ite)] << endl; //store name of the pseudopotential file and it's read options
                    if(_ppLData.find(*ite)==_ppLData.end()) {
                        cerr << "Warning: Element "<<(*ite)<<" has no angular momentum data (<l></l>) specified in CPMD options file.\n\tAttempting to read it from basis set. This may produce errors.\n" << flush;
                        if(_basisset_name.empty()){
                            cerr << "Error: Basis set file not specified.\n" << flush;
                            throw std::runtime_error("Encountered element with no angular momentum data.\n");
                        }
                        else{
                            BasisSet _bs;
                            _bs.LoadBasisSet(_basisset_name);
                            int Lmax = 0;
                            //find Lmax by checking all shells of the element
                            Element* el=_bs.getElement(*ite);
                            for (Element::ShellIterator its = el->firstShell(); its != el->lastShell(); its++) {
                                        int Ls=(*its)->getLmax();
                                        if(Lmax<Ls) Lmax=Ls;
                            }
                            _com_file << "   "<<Lmax<<" "<<Lmax<<" "<<Lmax<< endl; //LMAX LOC SKIP
                        }
                    }
                    else{
                        _com_file << "   "<<_ppLData[(*ite)]<< endl; //LMAX LOC SKIP
                    }
                    _com_file << "   "<< _nAtomsOfElement[(*ite)] <<endl;  //# atoms of element
                    
                    //store atomic positions
                    for (sit = segments.begin(); sit != segments.end(); ++sit) {
                        Vind=0;
                        _atoms = (*sit)-> Atoms();
                        for (ait = _atoms.begin(); ait < _atoms.end(); ait++) {
                            if((*ait)->getElement().compare(*ite)==0){     //this element
                                vec pos = (*ait)->getQMPos();
                                _com_file << "   ";
                                _com_file << setw(12) << setiosflags(ios::fixed) << setprecision(5) << conv::nm2bohr*pos.getX() << "   ";
                                _com_file << setw(12) << setiosflags(ios::fixed) << setprecision(5) << conv::nm2bohr*pos.getY() << "   ";
                                _com_file << setw(12) << setiosflags(ios::fixed) << setprecision(5) << conv::nm2bohr*pos.getZ() << "   ";
                                _com_file << endl;
                                
                                //cache the mapping between VOTCA and CPMD atomic ordering
                                if(!atomOrderMapSet){
                                    VOTCA2CPMD_map[Vind]=Cind;
                                    CPMD2VOTCA_map[Cind]=Vind;
                                    CPMD2TYPE_map[Cind]=(*ite);
                                }
                                Cind++;
                            }
                            Vind++;
                        }
                    }
                }
            }
            
            //#warning "TODO: copy pseudopotentials to the _run_dir"
            _com_file << "&END" << endl;
            
            
            
            _com_file << endl;
            _com_file.close();
            
            
            //now write the output file for the second run, if necessary
            if(_projectWF && _optWF){
                _optWF=false;
                //force reading of restart file on second run
                bool old_rsrt=_rsrt;
                std::string old_rsrt_kwds = _rsrt_kwds;
                _rsrt = true;
                _rsrt_kwds = "WAVEFUNCTION COORDINATES LATEST";
                WriteInputFile(segments, orbitals_guess);
                //reset variables
                _optWF     = true;
                _rsrt      = old_rsrt;
                _rsrt_kwds = old_rsrt_kwds;
            }

            return true;
        }
        
        
        /**
         * Writes the basis set files to disk in a format that CPMD can understand
         */
        void Cpmd::WriteBasisSet(std::vector<ctp::Segment* > segments, ofstream &_com_file) {
            
            std::vector< ctp::Atom* > _atoms;
            std::vector< ctp::Atom* > ::iterator ait;
            std::vector< ctp::Segment* >::iterator sit;
            list<std::string> elements;
            
            BasisSet _bs;
            _bs.LoadBasisSet(_basisset_name);
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;

            for (sit = segments.begin(); sit != segments.end(); ++sit) {

                std::vector< ctp::Atom* > atoms = (*sit)-> Atoms();
                std::vector< ctp::Atom* >::iterator it;

                for (it = atoms.begin(); it < atoms.end(); it++) {

                    std::string element_name = (*it)->getElement();

                    list<std::string>::iterator ite;
                    ite = find(elements.begin(), elements.end(), element_name);

                    if (ite == elements.end()) {
                        elements.push_back(element_name);
                        CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: writing Gaussian basis for element "<< element_name << flush;

                        Element* element = _bs.getElement(element_name);
                        
                        std::string _short_el_file_name = element_name + "_" + _basisset_name + ".basis";
                        std::string _el_file_name = _run_dir + "/" + _short_el_file_name;
                        
                        
                        //write the element to the input file
                        _com_file << "*" << _short_el_file_name << " " << std::distance(element->firstShell(), element->lastShell()) << " GAUSSIAN"<<endl;
                        _com_file << "   ";
                        
                        //create the .basis file
                        ofstream _el_file;
                        _el_file.open(_el_file_name.c_str());
                        
                        //comment
                        _el_file << element_name << " with the "<< _basisset_name << " basis." << endl;
                        
                        //Lmax
                        int Lmax = 0;
                        //find Lmax by checking all shells of the element
                        for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                    int Ls=(*its)->getLmax();
                                    if(Lmax<Ls) Lmax=Ls;
                        }
                        _el_file << Lmax+1 << endl; //number of L-values in this .basis file
                        
                        
                        //sort shells by L
                        for (int L=0; L <= Lmax; L++)
                        {
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: sorting shells by L; L="<< L<<"; Lmax="<< Lmax << flush;
                            std::vector<Shell*> Lshells;
                            
                            int ndecays=0;
                            for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                Shell* shell = (*its);
                                int Ls=shell->getLmax();
                                if(shell->getType().size()>1){
                                    cerr << "CPMD does not support " << shell->getType() << " basis functions." << endl;
                                    cerr << "Please break the basis set into basis functions with only one L-value each." << endl << flush;
                                    CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: multi-L basis functions not supported." << flush;
                                    throw std::runtime_error("Unsupported basis function");
                                }

                                //For now assume all shells have only one L-value.
                                //Can decompose the basis set into such shells later, in another place.
                                //Any subsequent analysis will have to use the decomposed basis set too.
                                if (Ls==L) //this shell has the correct L
                                {
                                    ndecays+=shell->getSize();
                                    Lshells.push_back(shell);
                                    
                                    //write the shell's L-value to the input file
                                    _com_file << L << " ";
                                }
                            }
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: writing decays; L="<< L<<"; Lmax="<< Lmax << flush;
                            if(!Lshells.empty()){   //only write things if there are shells with this L
                                _el_file << "  Functions for l="<<L<<endl;
                                _el_file << "  " << Lshells.size()<< " " << ndecays << endl;
                                _el_file << endl;

                                //decays
                                ios::fmtflags old_settings = _el_file.flags();
                                _el_file << std::scientific << std::setprecision(6);
                                _el_file << "  ";
                                for (Element::ShellIterator its = Lshells.begin(); its !=  Lshells.end(); its++)
                                {
                                    Shell* shell = (*its);
                                    for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                        GaussianPrimitive* gaussian = *itg;
                                        _el_file << gaussian->decay << "\t";
                                    }
                                }
                                _el_file << endl;

                                //coefficients (scale*contraction)
                                int gs=0; //number of decays already handled
                                for (Element::ShellIterator its = Lshells.begin(); its !=  Lshells.end(); its++)
                                {
                                    Shell* shell = (*its);
                                    if(shell->getSize()!=0) //there are gaussians in this shell
                                    {
                                        int gi=0; //index of the current decay
                                        _el_file << "  ";
                                        //output zeros for all decays already handled
                                        for(gi=0; gi<gs; gi++)
                                        {
                                            _el_file << 0.0 << "\t";
                                        }
                                        //output coefficients for this shell's gaussians
                                        for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                            GaussianPrimitive* gaussian = *itg;
                                            _el_file << shell->getScale() * gaussian->contraction[L] << "\t";
                                            gi++;
                                        }
                                        gs+=shell->getSize();
                                        //output zeros till the end of decays
                                        for(;gi<ndecays; gi++)
                                        {
                                            _el_file << 0.0 << "\t";
                                        }
                                        _el_file << endl;
                                    }
                                }
                                _el_file.flags(old_settings);
                            }
                        }

                        _el_file << endl;
                        _el_file.close();
                        
                        _com_file<<endl;
                    }
                }
            }
        }
        

        /**
         * Runs the CPMD job.
         */
        bool Cpmd::Run(Orbitals* _orbitals) {
            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: Run()" << flush;
            if(_optWF && _projectWF){ //CPMD needs to run twice, once for _optWF and once for _projectWF
                //_optWF run:
                 CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: running [" << _executable << " " << _wfOpt_input_file_name << "]" << flush;
                 if (std::system(NULL)) {
                    std::string _command;
                    _command = "cd " + _run_dir + "; rm -f LocalError*.log; " + _executable + " " + _wfOpt_input_file_name + " | tee " + _wfOpt_log_file_name;
                    int check=std::system(_command.c_str());
	            if (check==-1){
    	                CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
    	                return false;
    	            }
                    
                    check = std::system("mv LocalError*.log LocalError_wfOpt*.log");
                    if (check==0){
    	                CTP_LOG(ctp::logWARNING, *_pLog) << "CPMD produced an error log. Moving it to LocalError_wfOpt*.log" << flush;
    	            }

                    if (CheckLogFile()) {
                        CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: finished wavefunction optimization job. Continuing to projection onto AOs." << flush;
                    } else {
                        CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: wavefunction optimization job failed" << flush;
                        return false;
                    }
                    
                } else {
                    CTP_LOG(ctp::logERROR, *_pLog) << _wfOpt_input_file_name << " failed to start. No shell accessible." << flush;
                    return false;
                }
                 
                _optWF = false;
                //continue as usual
            }
            
            //CPMD only needs to run once, or _optWF just finished running
            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: running [" << _executable << " " << _input_file_name << "]" << flush;

            if (std::system(NULL)) {
                std::string _command;
                _command = "cd " + _run_dir + "; rm -f LocalError*.log; " + _executable + " " + _input_file_name + ">" + _log_file_name;
                int check=std::system(_command.c_str());
                if (check==-1){
                    CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                    return false;
                }
                if (CheckLogFile()) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: finished job" << flush;
                    return true;
                } else {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: job failed" << flush;
                }
            } else {
                CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                return false;
            }

            return true;

        }
               
        
        

        /**
         * Cleans up after the CPMD job.
         */
        void Cpmd::CleanUp() {

            // cleaning up the generated files
            if (_cleanup.size() != 0) {

                CTP_LOG(ctp::logDEBUG, *_pLog) << "Removing " << _cleanup << " files" << flush;
                Tokenizer tok_cleanup(_cleanup, ", ");
                std::vector <std::string> _cleanup_info;
                tok_cleanup.ToVector(_cleanup_info);

                std::vector<std::string> ::iterator it;

                for (it = _cleanup_info.begin(); it != _cleanup_info.end(); ++it) {

                    if (*it == "inp") {
                        std::string file_name = _run_dir + "/" + _input_file_name;
                        remove(file_name.c_str());
                        if(_projectWF && _optWF)
                        {
                            remove(_wfOpt_input_file_name.c_str()); //also clean up the WF opt input
                        }
                    }

                    if (*it == "log") {
                        std::string file_name = _run_dir + "/" + _log_file_name;
                        remove(file_name.c_str());
                        if(_projectWF && _optWF)
                        {
                            remove(_wfOpt_log_file_name.c_str()); //also clean up the WF opt input
                        }
                    }

                    if (*it == "chk") {
                        std::string file_name = _run_dir + "/" + "LATEST";
                        remove(file_name.c_str());
                        //if user sets custom execution options for CPMD,
                        //they may override this name and its not easy to track.
                        file_name = _run_dir + "/" + "RESTART.1";
                        remove(file_name.c_str());
                    }

                    if (*it == "fort.7") {
                        std::string file_name = _run_dir + "/" + *it;
                        remove(file_name.c_str());
                        file_name = _run_dir + "/" + "OVERLAP";
                        remove(file_name.c_str());
                        file_name = _run_dir + "/" + "SPINDEN";
                        remove(file_name.c_str());
                        file_name = _run_dir + "/" + "WFNCOEF";
                        remove(file_name.c_str());
                        file_name = _run_dir + "/" + "CHOUT";
                        remove(file_name.c_str());
                    }

                    if (*it == "basis") {
                        std::vector<std::string> fileswithfileending;
                        boost::filesystem::recursive_directory_iterator fit(_run_dir);
                        boost::filesystem::recursive_directory_iterator endit;

                        while (fit != endit) {
                            if (boost::filesystem::is_regular_file(* fit) &&
                                    fit->path().extension() == *it)
                            {
                                fileswithfileending.push_back(
                                    fit->path().filename().string());
                            }
                            ++fit;
                        }
                        for (const auto filename : fileswithfileending) {
                            std::string file_name = _run_dir + "/" + filename;
                            remove(file_name.c_str());
                        }
                    }
                    
                    if (*it == "density") {
                        std::string file_name = _run_dir + "/" + "DENSITY";
                        remove(file_name.c_str());
                    }
                    
                    if (*it == "elpot") {
                        std::string file_name = _run_dir + "/" + "ELPOT";
                        remove(file_name.c_str());
                    }

                }
            }
            return;

        }



        bool Cpmd::CheckLogFile() {

            // check if the log file exists
            boost::filesystem::path arg_path;

            std::string _full_name = (arg_path / _run_dir / _log_file_name).c_str();
            if(_optWF && _projectWF){ //CPMD needs to run twice; this is the _optWF run
                _full_name = (arg_path / _run_dir / _wfOpt_log_file_name).c_str();
            }
            ifstream _input_file(_full_name.c_str());

            if (_input_file.fail()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "CPMD: " << _full_name << " is not found." << endl << flush;
                return false;
            };

            //Use brute force. Search every line for the termination string.
            //It doesn't appear at the very end, like in gaussian
            std::string::size_type self_energy_pos=std::string::npos;
            std::string _line;
            do {
                getline(_input_file, _line);
                self_energy_pos=_line.find("PROGRAM CPMD ENDED AT");
            } while (self_energy_pos==std::string::npos && !(_input_file.eof()));

            _input_file.close();

            if (self_energy_pos == std::string::npos) {
                CTP_LOG(ctp::logERROR, *_pLog) << "CPMD: " << _full_name << " is incomplete."<< endl << flush;
                return false;
            } else {
                CTP_LOG(ctp::logDEBUG,*_pLog) << "CPMD LOG is complete." <<endl << flush;
                return true;
            }
        }

        /**
         * Parses the CPMD Log file and stores data in the Orbitals object
         */
        bool Cpmd::ParseLogFile(Orbitals * _orbitals) {
            std::string _line;
            std::vector<std::string> results;
            std::vector<tools::vec> positions;
            
            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: parsing " << _log_file_name << flush;
            
            std::string _log_file_name_full = _log_file_name;
            if (_run_dir != "") _log_file_name_full = _run_dir + "/" + _log_file_name;
            
            // check if LOG file is complete
            if (!CheckLogFile()) return false;

            // save qmpackage name
            _orbitals->setQMpackage("cpmd");
            
            ifstream _input_file(_log_file_name_full.c_str());
            while (_input_file) {
                
                getline(_input_file, _line);
                boost::trim(_line);
                
                
                /*
                 * number of electrons
                 */
                std::string::size_type electrons_pos = _line.find("alpha electrons");
                if (electrons_pos != std::string::npos) {
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    int _number_of_electrons = (int) boost::lexical_cast<double>(results.back());
                    _orbitals->setNumberOfElectrons(_number_of_electrons);
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Alpha electrons: " << _number_of_electrons << flush;
                }
                
                /*
                 * atomic positions
                 */
                if (_line.find("*** ATOMS ***") != std::string::npos) {
                    getline(_input_file, _line);
                    do{
                        getline(_input_file, _line);
                        boost::trim(_line);
                        if(_line.find("******") != std::string::npos)
                            break;

                        boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                        vec v;
                        v.setX(boost::lexical_cast<double>(results[2]));
                        v.setY(boost::lexical_cast<double>(results[3]));
                        v.setZ(boost::lexical_cast<double>(results[4]));
                        v=v*tools::conv::bohr2ang;
                        positions.push_back(v); //store positions and core charges (later))
                    }while(true);
                }
                
                
                /*
                 * Occupied/unoccupied states
                 */
                if (_line.find("NUMBER OF STATES:") != std::string::npos) {
                    boost::trim(_line);
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    int nstates = (boost::lexical_cast<int>(results[3]));
                    getline(_input_file, _line);
                    boost::trim(_line);
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    int nelectrons = (int) (boost::lexical_cast<double>(results[3]));
                    _orbitals->setNumberOfLevels(nelectrons/2, nstates-(nelectrons/2));
                }
                

                
                
            }
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Done parsing" << flush;
            _input_file.close();
            
            //check that CPMD2TYPE_map is available
            if(CPMD2VOTCA_map==NULL){
                CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: Can not convert atom order from CPMD to VOTCA." << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: Please rerun with writing CPMD input (<tasks>input, parse</tasks>)." << flush;
                throw std::runtime_error("CPMD2TYPE_map unavailable, rerun with <tasks>input, parse</tasks>");
                exit(-1);
            }
            
            //store atoms to Orbitals in VOTCA's order
            for(unsigned int v=0; v<positions.size(); v++){
                int c = ConvAtomIndex_VOTCA2CPMD(v);
                _orbitals->AddAtom(CPMD2TYPE_map[c], positions[c].getX(), positions[c].getY(), positions[c].getZ(), 0); //core charges don't matter, VOTCA computes them on its own
            }
            
            if(_projectWF){
                //MO coefficient and overlap matrices
                if(!loadMatrices(_orbitals)) return false;

                //atom info
                if(!(_orbitals->hasQMAtoms())){ //no atoms defined for the orbitals
                    //lets fill them in, in CPMD's order

                    //iterate over elements
                    list<std::string>::iterator ite;
                    int i=0;
                    for (ite = _elements.begin(); ite != _elements.end(); ite++, i++) {
                        for(int a=0; a<_NA[i]; a++){
                            _orbitals->AddAtom(*ite, 0, 0, 0, _ZV[i]); //store core charge in the atomic charge field
                        }
                    }
                }
//                else
//                {
//                    cerr << "CPMD: _orbitals already has some atoms. Need to implement atom reordering for this case." << flush;
//                    CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: _orbitals already has some atoms. Need to implement atom reordering for this case." << flush;
//                    throw std::runtime_error("Unimplemented case");
//                    return false;
//                }
            }
            
            
            //basis set
            if(_orbitals->hasDFTbasis()){
                if(_orbitals->getDFTbasis().compare(_basisset_name)!=0){
                    cerr << "CPMD: _orbitals already has a basis set and it does not match the basis set CPMD was initialized with." << flush;
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: _orbitals already has a basis set and it does not match the basis set CPMD was initialized with." << flush;
                    throw std::runtime_error("Basis set mismatch");
                    return false;
                }
            }
            else{
                _orbitals->setDFTbasis(_basisset_name);
            }
            
#ifdef DEBUG
            /*
            {
                //double check completeness of projection
                cout << "\nMO linear independence diagonal components using non-reordered data from CPMD:"<<endl;
                ub::matrix<double> MO = _orbitals->MOCoefficients();
                const ub::matrix<double>& AO=_orbitals->AOOverlap();
                ub::range all_basis_funcs = ub::range(0, MO.size2());
                for (unsigned i = 0; i < MO.size1(); i++) {
                    ub::range ri = ub::range(i, i+1);
                    ub::matrix<double> Ci = ub::project(MO, ri, all_basis_funcs);
                    unsigned j = i;
                    ub::range rj = ub::range(j, j+1);
                    ub::matrix<double> Cj = ub::trans(ub::project(MO, rj, all_basis_funcs));
                    ub::matrix<double> SCj = ub::prod(AO,Cj);
                    cout << ub::prod(Ci,SCj)(0,0) <<"\t";
                }
                cout<<endl<<flush;
            }
            */
#endif
            
            
            //fix order for version 5 of .orb files
            ReorderOutput(_orbitals);
            
#ifdef DEBUG
            

            {   
                ios::fmtflags f( cout.flags() );
                
                std::vector< ctp::QMAtom* > atoms = _orbitals->QMAtoms();
                
//                cout << "\n Is wrapping done right? Atom positions inside the basis (Angstroms):" <<endl;
//                for (const auto a : atoms) {
//                    cout <<a->type<<"\t" << std::fixed << setw(6) << setprecision(3) << a->x <<"\t"<< a->y <<"\t"<< a->z <<"\t"<<endl;
//                }
//                cout<<endl<<flush;
//                cout.flags( f );
                
                
                tools::vec box(10.0*tools::conv::bohr2ang);
                for (const auto a : atoms) {
                        tools::vec r(a->x, a->y, a->z);
                        for(int k=0; k<3; k++)
                        {
                            r[k] = fmod(r[k], box[k]);
                            if(r[k]<0) r[k] = box[k]+r[k];
                        }
                        a->x=r[0];
                        a->y=r[1];
                        a->z=r[2];
                }
                
//                cout << "\n Wrapped positions inside the basis (Angstroms):" <<endl;
//                for (const auto a : atoms) {
//                    cout <<a->type<<"\t" << std::fixed << setw(6) << setprecision(3) << a->x <<"\t"<< a->y <<"\t"<< a->z <<"\t"<<endl;
//                }
//                cout<<endl<<flush;
//                cout.flags( f );
                

                const ub::matrix<double>& CPMD_AO=_orbitals->AOOverlap();
                
                BasisSet bs;
                bs.LoadBasisSet(_orbitals->getDFTbasis());
                AOBasis basis;
                basis.AOBasisFill(&bs, _orbitals->QMAtoms());
                AOOverlap overlap;
                overlap.Fill(basis);
                ub::matrix<double>& VOTCA_AO=overlap.Matrix();

                AOBasis p_basis;
                cout << "\n Box (Angstroms):"<< box[0]<<"\t"<< box[1] <<"\t"<< box[2] <<endl;
                p_basis.AOBasisFill(&bs, _orbitals->QMAtoms());
                AOOverlapPeriodic p_overlap;
                p_overlap.setBox(box*tools::conv::ang2bohr); //in Angstrom
                p_overlap.Fill(p_basis);
                const ub::matrix<double>& p_AO=p_overlap.Matrix();


//                cout << std::fixed << setw(10) << setprecision(7);
//                cout << "AO overlaps:"<<endl;
//                cout << "   \t    \t\tO S  \t\t        \t|\t    \t\tO Pz \t\t        \t|\t    \t\tO Py \t\t        \t|\t    \t\tO Px \t\t        "<<endl;
//                cout << "   \tCPMD\t\tVOTCA\t\tperiodic\t|\tCPMD\t\tVOTCA\t\tperiodic"<<endl;
//                const std::vector<AOShell*> shells=basis.getShells();
//                for (const auto s : shells) {
//                    for (unsigned j = s->getStartIndex(); j < s->getStartIndex() + s->getNumFunc(); j++) {
//                            cout << s->getName() <<" "<< s->getType()<< "\t" << CPMD_AO(0,j) <<"\t" << VOTCA_AO(0,j) <<"\t"<< p_AO(0,j) << "\t|\t"
//                                    << CPMD_AO(4,j) <<"\t" << VOTCA_AO(4,j) <<"\t"<< p_AO(4,j) << "\t|\t"
//                                    << CPMD_AO(5,j) <<"\t" << VOTCA_AO(5,j) <<"\t"<< p_AO(5,j) << "\t|\t"
//                                    << CPMD_AO(6,j) <<"\t" << VOTCA_AO(6,j) <<"\t"<< p_AO(6,j) <<endl;
//                    }
//                }
                
                cout << std::fixed << setw(10) << setprecision(7);
                cout << "MO linear independence:"<<endl;
                cout << "   \tCPMD\t\tVOTCA\t\tperiodic"<<endl;
                ub::matrix<double> MO = _orbitals->MOCoefficients();
                ub::range all_basis_funcs = ub::range(0, MO.size2());
                for (unsigned i = 0; i < MO.size1(); i++) {
                    ub::range ri = ub::range(i, i+1);
                    ub::matrix<double> Ci = ub::project(MO, ri, all_basis_funcs);
                    ub::matrix<double> Cj = ub::trans(Ci);
                    ub::matrix<double> CPMD_SCj = ub::prod(CPMD_AO,Cj);
                    ub::matrix<double> VOTCA_SCj = ub::prod(VOTCA_AO,Cj);
                    ub::matrix<double> p_SCj = ub::prod(p_AO,Cj);
                    cout << "   \t" << ub::prod(Ci,CPMD_SCj)(0,0) <<"\t" << ub::prod(Ci,VOTCA_SCj)(0,0) <<"\t"<< ub::prod(Ci,p_SCj)(0,0) << endl;
                }
                
                
                cout<<endl<<flush;
                cout.flags( f );
            }
            
            
            
            /*
            {
                ios::fmtflags f( cout.flags() );
                //check if reordering is correct
                cout << "\nMO linear independence diagonal components using reordered overlap and MO coefficients:"<<endl;
                ub::matrix<double> MO = _orbitals->MOCoefficients();
                const ub::matrix<double>& AO=_orbitals->AOOverlap();
                ub::range all_basis_funcs = ub::range(0, MO.size2());
                for (unsigned i = 0; i < MO.size1(); i++) {
                    ub::range ri = ub::range(i, i+1);
                    ub::matrix<double> Ci = ub::project(MO, ri, all_basis_funcs);
                    ub::matrix<double> Cj = ub::trans(Ci);
                    ub::matrix<double> SCj = ub::prod(AO,Cj);
                    cout << ub::prod(Ci,SCj)(0,0) <<"\t";
                }
                cout<<endl<<flush;
                
                cout << "\nReordered Overlap ("<<AO.size1()<<","<<AO.size2()<<")" <<endl;
                //for (unsigned i = 0; i < AO.size1(); i++) {
                for (unsigned i = 0; i < 1; i++) {
                    for (unsigned j = 0; j < AO.size2(); j++) {
                        cout << std::fixed << setw(6) << setprecision(3) << AO(i,j) <<" ";
                    }
                    cout<<endl<<endl;
                }
                cout<<endl<<flush;
                cout.flags( f );
            }
            {
                ios::fmtflags f( cout.flags() );
                cout << "\nMO linear independence diagonal components using reordered MO coefficients and non-periodic overlap matrix produced by VOTCA:"<<endl;
                ub::matrix<double> MO = _orbitals->MOCoefficients();
                
                BasisSet bs;
                bs.LoadBasisSet(_orbitals->getDFTbasis());
                AOBasis basis;
                basis.AOBasisFill(&bs, _orbitals->QMAtoms());
                
                AOOverlap overlap;
                //overlap.setBox(tools::vec(10.0,10.0,10.0)); //in Bohr
                overlap.Fill(basis);   //AOOverlapPeriodic will build an overlap matrix taking periodicity into account here
                ub::matrix<double>& AO=overlap.Matrix();
                ub::range all_basis_funcs = ub::range(0, MO.size2());
                for (unsigned i = 0; i < MO.size1(); i++) {
                    ub::range ri = ub::range(i, i+1);
                    ub::matrix<double> Ci = ub::project(MO, ri, all_basis_funcs);
                    ub::matrix<double> Cj = ub::trans(Ci);
                    ub::matrix<double> SCj = ub::prod(AO,Cj);
                    cout << ub::prod(Ci,SCj)(0,0) <<"\t";
                }
                cout<<endl<<flush;
                
                cout << "\nGenerated Overlap ("<<AO.size1()<<","<<AO.size2()<<")" <<endl;
                //for (unsigned i = 0; i < AO.size1(); i++) {
                for (unsigned i = 0; i < 1; i++) {
                    for (unsigned j = 0; j < AO.size2(); j++) {
                        cout << std::fixed << setw(6) << setprecision(3) << AO(i,j) <<" ";
                    }
                    cout<<endl<<endl;
                }
                cout<<endl<<flush;
                cout.flags( f );
                
                
                cout << "\n Is wrapping done right? Atom positions inside the basis (Angstroms):" <<endl;
                std::vector< ctp::QMAtom* > atoms = _orbitals->QMAtoms();
                for (const auto a : atoms) {
                    cout <<a->type<<"\t" << std::fixed << setw(6) << setprecision(3) << a->x <<"\t"<< a->y <<"\t"<< a->z <<"\t"<<endl;
                }
                cout<<endl<<flush;
                cout.flags( f );
                
                
                tools::vec box(10.0*tools::conv::bohr2ang);
                cout << "\n Box (Angstroms):"<< box[0]<<"\t"<< box[1] <<"\t"<< box[2] <<endl;
                cout << "\n Wrapped atom positions:" <<endl;
                for (const auto a : atoms) {
                    tools::vec r(a->x, a->y, a->z);
                    for(int k=0; k<3; k++)
                    {
                        r[k] = fmod(r[k], box[k]);
                        if(r[k]<0) r[k] = box[k]+r[k];
                    }
                    a->x=r[0];
                    a->y=r[1];
                    a->z=r[2];

                    cout <<a->type<<"\t" << std::fixed << setw(6) << setprecision(3) << a->x <<"\t"<< a->y <<"\t"<< a->z <<"\t"<<endl;
                }
                cout<<endl<<flush;
                cout.flags( f );
                
                
                
                cout << "\nGenerated Periodic Overlap from wrapped positions ("<<AO.size1()<<","<<AO.size2()<<")" <<endl;
                AOBasis p_basis;
                p_basis.AOBasisFill(&bs, _orbitals->QMAtoms());
                AOOverlapPeriodic p_overlap;
                p_overlap.setBox(box*tools::conv::ang2bohr); //in Angstrom
                p_overlap.Fill(p_basis);
                AO=p_overlap.Matrix();
                for (unsigned i = 0; i < 1; i++) {
                    for (unsigned j = 0; j < AO.size2(); j++) {
                        cout << std::fixed << setw(6) << setprecision(3) << AO(i,j) <<" ";
                    }
                    cout<<endl<<endl;
                }
                cout<<endl<<flush;
                cout.flags( f );
                


                
            }
            */
            
            
#endif
            
            return true;

        }
        
        bool Cpmd::loadMatrices(Orbitals * _orbitals)
        {
            //sizes of int and real in whichever Fortran compiler CPMD was compiled with
            //these should be read from the xml file
            int F_int_size=4;
            int F_real_size=8;
            
            int totAtoms=0;
            
            //check if WFNCOEF exists
            boost::filesystem::path arg_path;
            std::string _full_name = (arg_path / _run_dir / "WFNCOEF").c_str();
            ifstream wf_file(_full_name.c_str());
            if(wf_file.fail())
            {
                CTP_LOG(ctp::logERROR, *_pLog) << "CPMD: " << _full_name << " is not found." << endl << flush;
            }
            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: parsing " << _full_name << flush;
            
            //read WFNCOEF
            //variable names in all CAPS are variables from CPMD source code
            int count=0, endcount=0;
            wf_file.read((char*)&count, F_int_size); //bytes in this record
            int bl=count;

            int NATTOT=0;
            wf_file.read((char*)&NATTOT, F_int_size);    //number of basis functions
            bl-=F_int_size;
            _orbitals->setBasisSetSize(NATTOT);
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Basis functions: " << NATTOT << flush;

            _NSP=0;                  //number of atom types
            wf_file.read((char*)&_NSP, F_int_size);
            bl-=F_int_size;

            //double ZV[NSP];             //core charge
            //int NA[NSP], NUMAOR[NSP];
            if(_ZV!=NULL)
            {
                delete[] _ZV;
                delete[] _NA;
                delete[] _NUMAOR;
            }
            _ZV=new double[_NSP];
            _NA=new int[_NSP];
            _NUMAOR=new int[_NSP];
            if(_ZV==NULL || _NA==NULL || _NUMAOR==NULL)
                throw std::runtime_error("Memory allocation failed");
            for(int i=0; i<_NSP; i++)
            {
                wf_file.read((char*)&_ZV[i], F_real_size);    //core charge of atom type i
                wf_file.read((char*)&_NA[i], F_int_size);     //number of atoms of type i
                wf_file.read((char*)&_NUMAOR[i], F_int_size); //number of atomic orbitals of atom type i
                bl-=F_real_size+F_int_size*2;
                
                totAtoms+=_NA[i];
            }
            
            //check footer
            wf_file.read((char*)&endcount, F_int_size);
            if(bl!=0 || endcount!=count){ //number of bytes read was wrong
                cerr << "CPMD: " << "could not parse record in "<< _full_name << endl << flush;
                CTP_LOG(ctp::logERROR, *_pLog) << "CPMD: " << "could not parse record in "<< _full_name << flush;
                throw std::runtime_error("IO error");
                return false;
            }
                
            //NEW RECORD
            wf_file.read((char*)&count, F_int_size); //bytes in this record
            bl=count;

            int NUMORB=count/F_real_size/NATTOT;          //number of MOs (energy levels))
            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: number of energy levels: " << NUMORB << flush;
                //resize the coefficient matrix
            
            
            
            
            //map atomic orbitals to (CPMD) atom indeces so we can reorder the MO and Overlap matrices to VOTCA's atom order
            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: Reordering orbitals."<< flush;
            int* AO_CPMD2VOTCA_map= new int[NATTOT];
            int* VOTCA2numAOs= new int[totAtoms];   //number of AOs for each atom in VOTCA atomic order 
            int* VOTCA2firstAO= new int[totAtoms];  //index of first AO of each atom in VOTCA atomic order 
            if(AO_CPMD2VOTCA_map==NULL || VOTCA2numAOs==NULL || VOTCA2firstAO==NULL){
                throw std::runtime_error("Memory allocation failed");
            }
            int c=0;
            for(int i=0; i<_NSP; i++)
            {
                for(int a=0; a<_NA[i]; a++){
                    VOTCA2numAOs[ConvAtomIndex_CPMD2VOTCA(c)]=_NUMAOR[i];
                    c++; //increment CPMD atom index
                }
            }
                //find the beginning of each atom in the VOTCA order
            int orb=0;
            for(int v=0; v<totAtoms; v++){
                VOTCA2firstAO[v]=orb;
                orb+=VOTCA2numAOs[v];
            }
                //map each CMPD AO to a VOTCA AO
            c=0;    //CPMD atom index
            int co=0;
            int vo=0;
            for(int i=0; i<_NSP; i++)
            {
                for(int a=0; a<_NA[i]; a++){
                    for(int ofst=0; ofst<_NUMAOR[i]; ofst++){
                        int v = ConvAtomIndex_CPMD2VOTCA(c); //VOTCA atom index
                        vo = VOTCA2firstAO[v] + ofst; //VOTCA orbital index
                        AO_CPMD2VOTCA_map[co] = vo;
                        
                        co++;  //increment CPMD orbital index
                    }
                    c++; //increment CPMD atom index
                }
            }

            
            
            
            
            
            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: Reading MO coefficients."<< flush;
            ub::matrix<double> &mo_coefficients = _orbitals->MOCoefficients();
            mo_coefficients.resize(NUMORB, NATTOT);
            
            //mo_coefficients need to be in VOTCA's atomic order
            //AO reordering comes later
            double XXMAT;
            for(int i=0; i<NUMORB; i++){
                for(int j=0; j<NATTOT; j++){  
                    wf_file.read((char*)&XXMAT, F_real_size);
                    mo_coefficients(AO_CPMD2VOTCA_map[i],j)=XXMAT;
                    bl-=F_real_size;
                }
            }
            
            
            //check footer
            wf_file.read((char*)&endcount, F_int_size);
            if(bl!=0 || endcount!=count){ //number of bytes read was wrong
                cerr << "CPMD: " << "could not parse record in "<< _full_name << endl << flush;
                CTP_LOG(ctp::logERROR, *_pLog) << "CPMD: " << "could not parse record in "<< _full_name << endl << flush;
                throw std::runtime_error("IO error");
                return false;
            }
            
            wf_file.close();
            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: Done parsing" << flush;
            
            
            
            
            //check if OVERLAP exists
            _full_name = (arg_path / _run_dir / "OVERLAP").c_str();
            ifstream ov_file(_full_name.c_str());
            if(ov_file.fail())
            {
                CTP_LOG(ctp::logERROR, *_pLog) << "CPMD: " << _full_name << " is not found." << endl << flush;
            }
            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: parsing " << _full_name << flush;
            
            //read OVERLAP
            count=0, endcount=0;
            ov_file.read((char*)&count, 4); //bytes in this record
            bl=count;
            
            if(NATTOT*NATTOT!=count/8)
            {
                cerr << "CPMD: " << "Number of basis functions in the overlap and coefficient matrices do not match."<< endl << flush;
                CTP_LOG(ctp::logERROR, *_pLog) << "CPMD: " << "Number of basis functions in the overlap and coefficient matrices do not match."<< endl << flush;
                throw std::runtime_error("IO error");
                return false;
            }
            
            
                //resize the overlap matrix
            ub::symmetric_matrix<double> &overlap = _orbitals->AOOverlap();
            overlap.resize(NATTOT);
            
            //read
            //Overlap need to be in VOTCA's atomic order
            CTP_LOG(ctp::logDEBUG, *_pLog) << "CPMD: Reading Overlap matrix."<< flush;
            double XSMAT;
            for(int i=0; i<NATTOT; i++){
                for(int j=0; j<NATTOT; j++){  
                    ov_file.read((char*)&XSMAT, F_real_size);
                    overlap(AO_CPMD2VOTCA_map[i],AO_CPMD2VOTCA_map[j])=XSMAT;
                    bl-=F_real_size;
                }
            }
            
            //check footer
            ov_file.read((char*)&endcount, F_int_size);
            if(bl!=0 || endcount!=count){ //number of bytes read was wrong
                cerr << "CPMD: " << "could not parse record in "<< _full_name << endl << flush;
                CTP_LOG(ctp::logERROR, *_pLog) << "CPMD: " << "could not parse record in "<< _full_name << flush;
                throw std::runtime_error("IO error");
                return false;
            }
            
            ov_file.close();
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Done parsing" << flush;
            
            delete[] VOTCA2numAOs;
            delete[] VOTCA2firstAO;
            delete[] AO_CPMD2VOTCA_map;
            
            return true;
        }


        
        std::string Cpmd::FortranFormat(const double &number) {
            std::stringstream _ssnumber;
            if (number >= 0) _ssnumber << " ";
            _ssnumber << setiosflags(ios::fixed) << setprecision(8) << std::scientific << number;
            std::string _snumber = _ssnumber.str();
            boost::replace_first(_snumber, "e", "D");
            return _snumber;
        }



    }
}
