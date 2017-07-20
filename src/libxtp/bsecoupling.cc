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

#include <votca/xtp/bsecoupling.h>
#include <votca/tools/constants.h>
#include <boost/format.hpp>



namespace votca { namespace xtp {

namespace ub = boost::numeric::ublas;
using boost::format;

void BSECoupling::Initialize(Property* options){
    
    #if (GWBSE_DOUBLE)
        CTP_LOG(ctp::logDEBUG, *_pLog) <<  " Compiled with full double support" << flush;   
    #else
        CTP_LOG(ctp::logDEBUG, *_pLog) <<  " Compiled with float/double mixture (standard)" << flush;   
    #endif
    
    std::string key = Identify(); 
    _doSinglets=false;
    _doTriplets=false;
   _output_perturbation=false;
    
    
    _openmp_threads = 0;
    
    if ( options->exists( key + ".openmp") ) {
                 _openmp_threads = options->get(key + ".openmp").as<int> ();
            }
    
    string spintype   = options->get(key + ".spin").as<string> ();
        if(spintype=="all"){
            _doSinglets=true;
            _doTriplets=true;
        }
        else if(spintype=="triplet"){
            _doTriplets=true;
        }
        else if(spintype=="singlet"){
            _doSinglets=true;
        }
        else{
            throw std::runtime_error((boost::format("Choice % for type not known. Available singlet,triplet,all") % spintype).str());
        }
    _degeneracy = options->get(key + ".degeneracy").as<double> ();
     
     
   if ( options->exists( key + ".algorithm") ) {
                string algorithm = options->get(key + ".algorithm").as<string> ();
                 if(algorithm=="perturbation"){
                    _output_perturbation=true;
                 }
   }
    
   
    
        _levA  = options->get(key + ".moleculeA.states").as<int> ();
        _levB  = options->get(key + ".moleculeB.states").as<int> ();
        _occA  = options->get(key + ".moleculeA.occLevels").as<int> ();
        _occB  = options->get(key + ".moleculeB.occLevels").as<int> ();
        _unoccA  = options->get(key + ".moleculeA.unoccLevels").as<int> ();
        _unoccB  = options->get(key + ".moleculeB.unoccLevels").as<int> ();
        
   
        
        
}

void BSECoupling::addoutput(Property *_type_summary,Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB){
   
    string algorithm="full_diag";
    int methodindex=1;
    if (_output_perturbation){
        algorithm="perturbation";
        methodindex=0;
    }
    _type_summary->setAttribute("algorithm",algorithm);
    if (_doSinglets){
        Property *_singlet_summary = &_type_summary->add("singlets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB <_levB ; ++stateB ) {
               double JAB = getSingletCouplingElement( stateA , stateB, methodindex);
              
               Property *_coupling_summary = &_singlet_summary->add("coupling", (format("%1$1.6e") % JAB).str()); 
               double energyA = _orbitalsA->BSESingletEnergies()(stateA)*conv::hrt2ev;
               double energyB = _orbitalsB->BSESingletEnergies()(stateB)*conv::hrt2ev;
               _coupling_summary->setAttribute("excitonA", stateA);
               _coupling_summary->setAttribute("excitonB", stateB);
               _coupling_summary->setAttribute("energyA", (format("%1$1.6e") % energyA).str());
               _coupling_summary->setAttribute("energyB", (format("%1$1.6e") % energyB).str());
               _coupling_summary->setAttribute("pert", (format("%1$1.6e") % getSingletCouplingElement( stateA , stateB, 0)).str());
               _coupling_summary->setAttribute("diag", (format("%1$1.6e") % getSingletCouplingElement( stateA , stateB, 1)).str());
               
           } 
        }
    }
    
    //cout << JAB_triplet<<endl;
    //cout << JAB_triplet*conv::hrt2ev<<endl;
    if ( _doTriplets){
        Property *_triplet_summary = &_type_summary->add("triplets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB < _levA ; ++stateB ) {
               double JAB = getTripletCouplingElement( stateA , stateB,methodindex );
               //real_gwbse energyAD = getTripletDimerEnergy( stateA  );
               //real_gwbse energyBD = getTripletDimerEnergy( stateB  );
               Property *_coupling_summary = &_triplet_summary->add("coupling", (format("%1$1.6e") % JAB).str()); 
               double energyA = _orbitalsA->BSETripletEnergies()(stateA)*conv::hrt2ev;
               double energyB = _orbitalsB->BSETripletEnergies()(stateB)*conv::hrt2ev;
               _coupling_summary->setAttribute("excitonA", stateA);
               _coupling_summary->setAttribute("excitonB", stateB);
               _coupling_summary->setAttribute("energyA", (format("%1$1.6e") % energyA).str());
               _coupling_summary->setAttribute("energyB", (format("%1$1.6e") % energyB).str());
               _coupling_summary->setAttribute("pert", (format("%1$1.6e") % getTripletCouplingElement( stateA , stateB, 0)).str());
               _coupling_summary->setAttribute("diag", (format("%1$1.6e") % getTripletCouplingElement( stateA , stateB, 1)).str());
              
           } 
        }
    }       
}


double BSECoupling::getSingletCouplingElement( int levelA, int levelB, int methodindex) {
    return JAB_singlet[methodindex]( levelA  , levelB +  _levA ) * votca::tools::conv::hrt2ev;
}



double BSECoupling::getTripletCouplingElement( int levelA, int levelB, int methodindex) {

    return JAB_triplet[methodindex]( levelA  , levelB + _levA ) * votca::tools::conv::hrt2ev;
}


/**
 * \brief evaluates electronic couplings  
 *   
 * @param _orbitalsA molecular orbitals of molecule A
 * @param _orbitalsB molecular orbitals of molecule B
 * @param _orbitalsAB molecular orbitals of the dimer AB
 * @param _JAB matrix with electronic couplings
 * @return false if failed
 */
bool BSECoupling::CalculateCouplings(Orbitals* _orbitalsA, Orbitals* _orbitalsB, Orbitals* _orbitalsAB) {
       CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Calculating exciton couplings" << flush;
     // set the parallelization 
    #ifdef _OPENMP
    
    if ( _openmp_threads > 0 ) omp_set_num_threads(_openmp_threads);      
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << " Using "<< omp_get_max_threads()<<" threads" << flush;
    #endif
    
    
    
    _orbitalsAB->setCoupledExcitonsA(_levA);
    _orbitalsAB->setCoupledExcitonsB(_levB);
    //check to see if ordering of atoms agrees
    const std::vector<ctp::QMAtom*> atomsA=_orbitalsA->QMAtoms();
    const std::vector<ctp::QMAtom*> atomsB=_orbitalsB->QMAtoms();
    const std::vector<ctp::QMAtom*> atomsAB=_orbitalsAB->QMAtoms();
    
    for (unsigned i=0;i<atomsAB.size();i++){
        ctp::QMAtom* dimer=atomsAB[i];
        ctp::QMAtom* monomer=NULL;
        if (i<atomsA.size()){
            monomer=atomsA[i];
        }
        else if (i<atomsB.size()+atomsA.size() ){
            monomer=atomsB[i-atomsA.size()];
        }
        else{
            throw runtime_error((boost::format("Number of Atoms in dimer %3i and the two monomers A:%3i B:%3i does not agree") %atomsAB.size() %atomsA.size() %atomsB.size()).str());
        }
        
        if(monomer->type != dimer->type){
            throw runtime_error("\nERROR: Atom types do not agree in dimer and monomers\n");
        }
        if(std::abs(monomer->x-dimer->x)>0.001 || std::abs(monomer->y-dimer->y)>0.001 || std::abs(monomer->z-dimer->z)>0.001){
            CTP_LOG(ctp::logERROR,*_pLog) << "======WARNING=======\n Coordinates of monomers and dimer atoms do not agree, do you know what you are doing?\n " << flush;
            break;
        }
        
    }
    
    
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
    
    if ( ( _basisA == 0 ) || ( _basisB == 0 ) ) {
        CTP_LOG(ctp::logERROR,*_pLog) << "Basis set size is not stored in monomers" << flush;
        return false;
    }

    // number of levels stored in monomers
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
        
    // get exciton information of molecule A
    int _bseA_cmax        = _orbitalsA->getBSEcmax();
    int _bseA_cmin        = _orbitalsA->getBSEcmin();
    int _bseA_vmax        = _orbitalsA->getBSEvmax();
    int _bseA_vmin        = _orbitalsA->getBSEvmin();
    int _bseA_vtotal      = _bseA_vmax - _bseA_vmin +1 ;
    int _bseA_ctotal      = _bseA_cmax - _bseA_cmin +1 ;
    int _bseA_size        = _bseA_vtotal * _bseA_ctotal;
    int _bseA_singlet_exc = _orbitalsA->BSESingletCoefficients().size2();
    int _bseA_triplet_exc = _orbitalsA->BSETripletCoefficients().size2();

    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   molecule A has " << _bseA_singlet_exc << " singlet excitons with dimension " << _bseA_size << flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   molecule A has " << _bseA_triplet_exc << " triplet excitons with dimension " << _bseA_size << flush;
    
    // now, two storage assignment matrices for two-particle functions
    ub::matrix<int> _combA;
    _combA.resize(_bseA_size,2);
    int _cnt = 0;
    for ( int _v = 0; _v < _bseA_vtotal; _v++){
        for ( int _c = 0; _c < _bseA_ctotal; _c++){
            _combA(_cnt,0) = _v;
            _combA(_cnt,1) = _bseA_vtotal + _c;
            _cnt++;
        }
    }
    
    // get exciton information of molecule B
    int _bseB_cmax        = _orbitalsB->getBSEcmax();
    int _bseB_cmin        = _orbitalsB->getBSEcmin();
    int _bseB_vmax        = _orbitalsB->getBSEvmax();
    int _bseB_vmin        = _orbitalsB->getBSEvmin();
    int _bseB_vtotal      = _bseB_vmax - _bseB_vmin +1 ;
    int _bseB_ctotal      = _bseB_cmax - _bseB_cmin +1 ;
    int _bseB_size        = _bseB_vtotal * _bseB_ctotal;
    int _bseB_singlet_exc = _orbitalsB->BSESingletCoefficients().size2();
    int _bseB_triplet_exc = _orbitalsB->BSETripletCoefficients().size2();
    

    
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   molecule B has " << _bseB_singlet_exc << " singlet excitons with dimension " << _bseB_size << flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   molecule B has " << _bseB_triplet_exc << " triplet excitons with dimension " << _bseB_size << flush;
    
    // now, two storage assignment matrices for two-particle functions
    ub::matrix<int> _combB;
    _combB.resize(_bseB_size,2);
    _cnt = 0;
    for ( int _v = 0; _v < _bseB_vtotal; _v++){
        for ( int _c = 0; _c < _bseB_ctotal; _c++){
            _combB(_cnt,0) = _bseA_vtotal + _bseA_ctotal + _v;
            _combB(_cnt,1) = _bseA_vtotal + _bseA_ctotal + _bseB_vtotal + _c;
            _cnt++;
        }
    }
    
    if(_levA>_bseA_singlet_exc){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of excitons you want is greater than stored for molecule A. Setting to max number available" << flush; 
        _levA=_bseA_singlet_exc;
    }
    if(_levB>_bseB_singlet_exc){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of excitons you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _levB=_bseB_singlet_exc;
    }
    
    
    if(_levA>_bseA_singlet_exc){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of Frenkel states you want is greater than stored for molecule A. Setting to max number available" << flush; 
        _levA=_bseA_singlet_exc;
    }
    if(_levB>_bseB_singlet_exc){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of Frenkel states you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _levB=_bseB_singlet_exc;
    }
    
    if(_unoccA>_bseA_ctotal){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of occupied orbitals in molecule A for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _unoccA=_bseA_ctotal;
    }
    else if (_unoccA<0){
        _unoccA=_bseA_ctotal;
    }
    if(_unoccB>_bseB_ctotal){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of occupied orbitals in molecule B for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _unoccB=_bseB_ctotal;
    }
    else if (_unoccB<0){
        _unoccB=_bseB_ctotal;
    }
    
    if(_occA>_bseA_vtotal){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of unoccupied orbitals in molecule A for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _occA=_bseA_vtotal;
    }
    else if (_occA<0){
        _occA=_bseA_vtotal;
    }
    if(_occB>_bseB_vtotal){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of unoccupied orbitals in molecule B for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _occB=_bseB_vtotal;
    }else if (_occB<0){
        _occB=_bseB_vtotal;
    }
    
    
    if(_orbitalsA->getQMpackage()!=_orbitalsB->getQMpackage() || _orbitalsA->getQMpackage()!=_orbitalsAB->getQMpackage()){
        throw runtime_error("Qmpackages in Orbfiles were made using different QMPackages, that doe snot work at the moment due to not reorddeering the MOCOefficients/Overlapmatrix.\n Redo the calculation please. ");
    }
    
    // get exciton information of pair AB
    int _bseAB_cmax = _orbitalsAB->getBSEcmax();
    int _bseAB_cmin = _orbitalsAB->getBSEcmin();
    int _bseAB_vmax = _orbitalsAB->getBSEvmax();
    int _bseAB_vmin = _orbitalsAB->getBSEvmin();
    int _bseAB_vtotal = _bseAB_vmax - _bseAB_vmin +1 ;
    int _bseAB_ctotal = _bseAB_cmax - _bseAB_cmin +1 ;
    int _bseAB_size   = _bseAB_vtotal * _bseAB_ctotal;
    // check if electron-hole interaction matrices are stored
    if ( ! _orbitalsAB->hasEHinteraction() ){
        CTP_LOG(ctp::logERROR,*_pLog) << "BSE EH int not stored in dimer " << flush;
        return false;
    }
    const ub::matrix<real_gwbse>&    _eh_d = _orbitalsAB->eh_d(); 
    const ub::matrix<real_gwbse>&    _eh_x = _orbitalsAB->eh_x(); 
    
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   dimer AB has BSE EH interaction (direct)   with dimension " << _eh_d.size1() << " x " <<  _eh_d.size2() << flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   dimer AB has BSE EH interaction (exchange) with dimension " << _eh_x.size1() << " x " <<  _eh_x.size2() << flush;
    // now, two storage assignment matrices for two-particle functions
    ub::matrix<int> _combAB;
    _combAB.resize(_bseAB_size,2);
    _cnt = 0;
    for ( int _v = 0; _v < _bseAB_vtotal; _v++){
        for ( int _c = 0; _c < _bseAB_ctotal; _c++){
            //_combAB(_cnt,0) = _v;
            //_combAB(_cnt,1) = _bseAB_vtotal + _c;
            
            _combAB(_cnt,0) = _bseAB_vmin + _v;
            _combAB(_cnt,1) = _bseAB_vmin + _bseAB_vtotal + _c;
            
            _cnt++;
        }
    }
    
    

    
    // DFT levels of monomers can be reduced to those used in BSE
    _levelsA = _bseA_vtotal + _bseA_ctotal;
    _levelsB = _bseB_vtotal + _bseB_ctotal;
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   levels used in BSE of molA: " << _bseA_vmin << " to " << _bseA_cmax << " total: " << _bseA_vtotal + _bseA_ctotal <<  flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   levels used in BSE of molB: " << _bseB_vmin << " to " << _bseB_cmax << " total: " << _bseB_vtotal + _bseB_ctotal <<  flush;
    
    
    if ( ( _levelsA == 0 ) || (_levelsB == 0) ) {
        CTP_LOG(ctp::logERROR,*_pLog) << "No information about number of occupied/unoccupied levels is stored" << flush;
        return false;
    } 
    
    //       | Orbitals_A          0 |      | Overlap_A |     
    //       | 0          Orbitals_B |  X   | Overlap_B |  X  Transpose( Orbitals_AB )
    ub::zero_matrix<double> zeroB( _levelsA, _basisB ) ;
    ub::zero_matrix<double> zeroA( _levelsB, _basisA ) ;
    ub::matrix<double> _psi_AxB ( _levelsA + _levelsB, _basisA + _basisB  );
    
    // constructing merged orbitals
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = ub::project( _orbitalsA->MOCoefficients() , ub::range(_bseA_vmin, _bseA_cmax+1) , ub::range ( 0, _basisA ));
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = ub::project( _orbitalsB->MOCoefficients(), ub::range(_bseB_vmin, _bseB_cmax+1) , ub::range ( 0, _basisB )); 
    
    // psi_AxB * S_AB * psi_AB
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   projecting monomer onto dimer orbitals" << flush; 
    
    if ( !_orbitalsAB->hasAOOverlap() ) {
            CTP_LOG(ctp::logERROR,*_pLog) << "Overlap matrix is not stored"; 
            return false;
    }
    //convert to full matrix from symmetric
    ub::matrix<double> _overlapAB =_orbitalsAB->AOOverlap();  
    ub::matrix<double> _psi_AB = ub::prod( _overlapAB,ub::trans(_orbitalsAB->MOCoefficients()) );  
    ub::matrix<double> _psi_AxB_dimer_basis = ub::prod( _psi_AxB, _psi_AB );  
    _psi_AB.resize(0,0);
    _overlapAB.resize(0,0);
    //cout<< "_psi_AxB_dimer"<<endl;
    unsigned int LevelsA = _levelsA;
    for (unsigned i=0;i<_psi_AxB_dimer_basis.size1();i++){
        double mag=0.0;
        for (unsigned j=0;j<_psi_AxB_dimer_basis.size2();j++){
            mag+=_psi_AxB_dimer_basis(i,j)*_psi_AxB_dimer_basis(i,j);
            
    }
         if (mag<0.95){
            int monomer = 0;
            int level = 0;
            if ( i < LevelsA ) {
                monomer = 1;
                level   = _bseA_vmin + i;
            } else {
                monomer = 2;
                level   = _bseB_vmin + i -_levelsA;
                
            }
            CTP_LOG(ctp::logERROR,*_pLog) << "\nERROR: " << i << " Projection of orbital " << level << " of monomer " << monomer << " on dimer is insufficient,mag="<<mag<<" maybe the orbital order is screwed up, otherwise increase dimer basis.\n"<<flush;
        }
    }
   
    
    //notation AB is CT states with A+B-, BA is the counterpart
    //Setting up CT-states:
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   Setting up CT-states" << flush; 
    //Number of A+B- states
    int noAB=_occA*_unoccB;
    //Number of A-B+ states
    int noBA=_unoccA*_occB;
    
    
    ub::matrix<int> comb_CTAB;
    comb_CTAB.resize(noAB,2);
    _cnt = 0;
    

    
    
    // iterate A over occupied, B over unoccupied
    int v_start=_bseA_vtotal-_occA;
    for ( int _v = v_start; _v < _bseA_vtotal; _v++){
        for ( int _c = 0; _c <_unoccB; _c++){            
            comb_CTAB(_cnt,0) =_v;
            comb_CTAB(_cnt,1) = _bseA_vtotal+_bseA_ctotal+_bseB_vtotal + _c;
           
            _cnt++;
        }
    }
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  <<"  "<<noAB <<" CT states A+B- created" << flush;
 
    ub::matrix<int> comb_CTBA;
    comb_CTBA.resize(noBA,2);
    _cnt = 0;
    // iterate A over unoccupied, B over occupied
    v_start=_bseB_vtotal-_occB;
    for ( int _v = v_start; _v < _bseB_vtotal; _v++){
        for ( int _c = 0; _c <_unoccA; _c++){            
            comb_CTBA(_cnt,0) =_bseA_vtotal+_bseA_ctotal+_v;
            comb_CTBA(_cnt,1) = _bseA_vtotal+ _c;
            
            _cnt++;
        }
    }
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  <<"  "<<noBA <<" CT states B+A- created" << flush;
    
    
    
    // these 4 matrixes, matrix(i,j) contains the j-th dimer MO component of the i-th excitation
   
    ctAB.resize(noAB,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_CT = 0 ; _i_CT < noAB ; _i_CT++){
    for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
        ctAB(_i_CT,_i_bseAB)=_psi_AxB_dimer_basis( comb_CTAB(_i_CT,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( comb_CTAB(_i_CT,1), _combAB( _i_bseAB,1) );
        }
    }
    
    
    ctBA.resize(noBA,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_CT = 0 ; _i_CT < noBA ; _i_CT++){
    for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
        ctBA(_i_CT,_i_bseAB)=_psi_AxB_dimer_basis( comb_CTBA(_i_CT,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( comb_CTBA(_i_CT,1), _combAB( _i_bseAB,1) );
        }
    }
    
      
    // some more convenient storage
    
    
    _kap.resize(_bseA_size,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) );
            
        }
    }

    
    

    _kbp.resize(_bseB_size,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) );
        }
    }
    
    // Same routines but also take <v|c'> <c|v'> projections into account 
    /*
 
    _kap.resize(_bseA_size,_bseAB_size);
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) )+
              _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,1) )     ;
            
        }
    }

    

    
   
    _kbp.resize(_bseB_size,_bseAB_size);
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) )+
                    _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,1) );
        }
    }
    */ 
  

    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   construct projection of product functions " << flush; 

   
    
        
    //     cout << "Size of _kap " << _kap.size1() << " : " <<  _kap.size2() << "\n" << flush; 
    //     cout << "Size of _kbp " << _kbp.size1() << " : " <<  _kbp.size2() << "\n" << flush; 
    _psi_AxB_dimer_basis.resize(0,0);
    _combAB.resize(0,0);
    _combA.resize(0,0);
    _combB.resize(0,0);
    // now the different spin types
            if (_doSinglets) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Evaluating singlets" << flush;
                // get singlet BSE Hamiltonian from _orbitalsAB
                ub::matrix<double> _Hamiltonian_AB = _eh_d + 2.0 * _eh_x;
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Setup Hamiltonian" << flush;
                ub::matrix<real_gwbse> temp = ub::project(_orbitalsA->BSESingletCoefficients(),
                        ub::range(0, _orbitalsA->BSESingletCoefficients().size1()), ub::range(0, _levA));

                const ub::matrix<double> _bseA_T = ub::trans(temp);
                temp = ub::project(_orbitalsB->BSESingletCoefficients(),
                        ub::range(0, _orbitalsB->BSESingletCoefficients().size1()), ub::range(0, _levB));
                const ub::matrix<double> _bseB_T = ub::trans(temp);
                temp.resize(0, 0);
                
                JAB_singlet = ProjectExcitons(_bseA_T, _bseB_T, _Hamiltonian_AB);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   calculated singlet couplings " << flush;
            }



            if (_doTriplets) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Evaluating triplets" << flush;
                // get triplet BSE Hamiltonian from _orbitalsAB
                ub::matrix<double> _Hamiltonian_AB = _eh_d;
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "  Converted Hamiltonian to double" << flush;
                ub::matrix<real_gwbse> temp = ub::project(_orbitalsA->BSETripletCoefficients(),
                        ub::range(0, _orbitalsA->BSETripletCoefficients().size1()), ub::range(0, _levA));
                const ub::matrix<double> _bseA_T = ub::trans(temp);
                temp = ub::project(_orbitalsB->BSETripletCoefficients(), ub::range(0, _orbitalsB->BSETripletCoefficients().size1()), ub::range(0, _levB));
                const ub::matrix<double> _bseB_T = ub::trans(temp);
                temp.resize(0, 0);

                JAB_triplet = ProjectExcitons(_bseA_T, _bseB_T, _Hamiltonian_AB);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   calculated triplet couplings " << flush;
            }
    
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Done with exciton couplings" << flush;
    return true;   
};


std::vector< ub::matrix<double> > BSECoupling::ProjectExcitons(const ub::matrix<double>& _bseA_T, const ub::matrix<double>& _bseB_T, 
                                  ub::matrix<double>& _H){
    
    
    
    
     // get projection of monomer excitons on dimer product functions
     ub::matrix<double> _proj_excA = ub::prod( _bseA_T, _kap);
     ub::matrix<double> _proj_excB = ub::prod( _bseB_T, _kbp);
     
     _bseA_exc = _proj_excA.size1();
     _bseB_exc = _proj_excB.size1();
     _bse_exc=_bseA_exc+_bseB_exc;
     
     
     
     unsigned _ctAB=ctAB.size1();
     
     unsigned _ctBA=ctBA.size1();
     _ct=_ctAB+_ctBA;
     unsigned nobasisfunc=_H.size1();
     
     
     
     ub::matrix<double> fe_states=ub::matrix<double>(_bse_exc,nobasisfunc);
     ub::project(fe_states, ub::range ( 0, _bseA_exc ),ub::range (0, nobasisfunc )  )=_proj_excA;
     ub::project(fe_states, ub::range ( _bseA_exc, _bse_exc ),ub::range (0, nobasisfunc )  )=_proj_excB;
      
     ub::matrix<double> ct_states=ub::matrix<double>(_ct,nobasisfunc);
     
     //cout<< _ct<< "ct states"<<endl;
    if(_ct>0){ 
     //orthogonalize ct-states with respect to the FE states. 
       CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << " Orthogonalizing CT-states with respect to FE-states" << flush;
   
     if(_ctAB>0){
     ub::project(ct_states, ub::range ( 0 , _ctAB ) ,ub::range (0,nobasisfunc ) )=ctAB;
    }
    if(_ctBA>0){
     ub::project(ct_states, ub::range ( _ctAB, _ct ),ub::range (0, nobasisfunc )  )=ctBA;
     }
       
       
        //orthogonalize ct-states with respect to FE states
     ub::matrix<double> overlaps=ub::prod(fe_states,ub::trans(ct_states));
     
     //cout << "overlap"<< overlaps.size1()<< "x"<<overlaps.size2()<<endl;
     //cout << overlaps<<endl;
     ub::matrix<double> correction=ub::prod(ub::trans(overlaps),fe_states);
     //cout << "correction"<< correction.size1()<< "x"<<correction.size2()<<endl;
    // cout << "ct_states"<< ct_states.size1()<< "x"<<ct_states.size2()<<endl;
  
   

    ct_states=ct_states-correction;    

     overlaps.resize(0,0);
     correction.resize(0,0);
     //normalize
    
     for (unsigned i=0;i<_ct;i++){
         double norm=0.0;
         for (unsigned j=0;j<nobasisfunc;j++){
         norm+=ct_states(i,j)*ct_states(i,j);    
         }
         //cout << "norm ["<<i<<"]:" <<norm<<endl;
         norm=1/std::sqrt(norm);
         if(norm<0.95){
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << " WARNING: CT-state "<< i<< " norm is only"<< norm << flush; 
         }
         for (unsigned j=0;j<nobasisfunc;j++){
            ct_states(i,j)=norm*ct_states(i,j);    
         }
         
     }
    //cout <<ub::prod(fe_states,ub::trans(ct_states))<<endl; 
      
    } 

     
    
    
     ub::matrix<double> projection =ub::zero_matrix<double>(_bse_exc+_ct,nobasisfunc);
     
     
     CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << " merging projections into one vector  " << flush;
    
  ub::project(projection, ub::range (0 , _bse_exc) ,ub::range (0,nobasisfunc ) )=fe_states;
   
     if(_ct>0){
    ub::project(projection, ub::range ( _bse_exc , _bse_exc+_ct ) ,ub::range (0,nobasisfunc ) )=ct_states;
     }
      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   Setting up coupling matrix size "<< _bse_exc +_ct<<"x"<<_bse_exc +_ct << flush;
     // matrix _J
     
    //  E_A         J_AB        J_A_ABCT        J_A_BACT
    //  J_BA        E_B         J_B_ABCT        J_B_BACT
    //  J_ABCT_A    J_ABCT_B    E_ABCT          J_ABCT_BACT
    //  J_BACT_A   J_BACT_B    J_BACT_ABCT     E_BACT
     
     // I think this only works for hermitian/symmetric H so only in TDA
     // setup J
     
    
     
     ub::matrix<double> _temp=ub::prod(_H,ub::trans(projection));
     _H.resize(0,0);
     ub::matrix<double> _J_dimer=ub::prod(projection,_temp);
     _temp.resize(0,0);
     

    
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   Setting up overlap matrix size "<< _bse_exc +_ct<<"x"<<_bse_exc +_ct << flush;
     // setup S
    
    ub::matrix<double> _S_dimer=ub::prod(projection,ub::trans(projection));
    
    projection.resize(0,0);
    if(tools::globals::verbose &&  _bse_exc+_ct<100){
         CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "_J_dimer[Ryd]"<<flush;
     
     CTP_LOG(ctp::logDEBUG, *_pLog) << _J_dimer<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "_S_dimer"<<flush;
     
     CTP_LOG(ctp::logDEBUG, *_pLog) << _S_dimer<<flush;
      CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    }
   
    
    double small=linalg_loewdin(_J_dimer,_S_dimer);
    
    if(tools::globals::verbose && _bse_exc+_ct<100){
         CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << "_J_ortho[Ryd]"<<flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << _J_dimer<<flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << "_S-1/2"<<flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << _S_dimer<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    }
     CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   Smallest value of dimer overlapmatrix is "<< small<< flush;
     
    std::vector< ub::matrix<double> >_J;
     
     CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   Running Perturbation algorithm"<< flush;
    _J.push_back( Perturbation(_J_dimer));
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "    Running Projection algorithm"<< flush;
    _J.push_back( Fulldiag(_J_dimer));
    
    
       if(tools::globals::verbose){
     CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "Jeff_pert[Hrt]"<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << _J[0]<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "Jeff_diag[Hrt]"<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << _J[1]<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     }
      
     return _J;
}

ub::matrix<double> BSECoupling::Perturbation(const ub::matrix<double>& _J_dimer){
    
    ub::matrix<double> _J = ub::zero_matrix<double>(_bse_exc, _bse_exc);
    bool _diag_ct = true;
    ub::matrix<double> _J_result=_J_dimer;
    if (_ct > 0 && _diag_ct) {

        ub::matrix<double> transformation = ub::identity_matrix<double>(_bse_exc + _ct, _bse_exc + _ct);
        ub::vector<double> eigenvalues_ct;



        ub::matrix<double> Ct = ub::project(_J_dimer, ub::range(_bse_exc, _bse_exc + _ct), ub::range(_bse_exc, _bse_exc + _ct));
        linalg_eigenvalues(eigenvalues_ct, Ct);
        ub::project(transformation, ub::range(_bse_exc, _bse_exc + _ct), ub::range(_bse_exc, _bse_exc + _ct)) = Ct;

        Ct.resize(0, 0);

        if (tools::globals::verbose) {

            CTP_LOG(ctp::logDEBUG, *_pLog) << "FE state hamiltonian" << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << ub::project(_J_dimer, ub::range(0, _bse_exc), ub::range(0, _bse_exc)) << flush;
            if (_ct > 0) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "eigenvalues of CT states" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << eigenvalues_ct << flush;
            }

        }

        ub::matrix<double> _temp = ub::prod(_J_dimer, transformation);
        _J_result = ub::prod(ub::trans(transformation), _temp);
        if (tools::globals::verbose && _bse_exc + _ct < 100) {
            CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "_J_ortho[Hrt] CT-state diag" << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << _J_result << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
        }
    }
    for (int stateA = 0; stateA < _levA; stateA++) {
        double Ea = _J_result(stateA, stateA);
        for (int stateB = 0; stateB < _levB; stateB++) {
            int stateBd = stateB + _bseA_exc;
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Calculating coupling between exciton A" << stateA + 1 << " and exciton B" << stateB + 1 << flush;
            double J = _J_result(stateA, stateBd);

            double Eb = _J_result(stateBd, stateBd);
            for (unsigned k = _bse_exc; k < (_bse_exc + _ct); k++) {
                double Eab = _J_result(k, k);
                if (std::abs(Eab - Ea) < 0.001) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "Energydifference between state A " << stateA + 1 << "and CT state " << k + 1 << " is " << Eab - Ea << "[Hrt]" << flush;
                }
                if (std::abs(Eab - Eb) < 0.001) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "Energydifference between state B " << stateB + 1 << "and CT state " << k + 1 << " is " << Eab - Eb << "[Hrt]" << flush;

                }
                J += 0.5 * _J_result(k, stateA) * _J_result(k, stateBd)*(1 / (Ea - Eab) + 1 / (Eb - Eab)); // Have no clue why 0.5
            }
            _J(stateA, stateBd) = J;
            _J(stateBd, stateA) = J;


        }
    }

            
    return _J;
}


ub::matrix<double> BSECoupling::Fulldiag(const ub::matrix<double>& _J_dimer){
    ub::matrix<double> _J = ub::zero_matrix<double>(_bse_exc, _bse_exc);
   

    ub::vector<double> _J_eigenvalues;
    ub::matrix<double> J_eigenvectors;

    linalg_eigenvalues(_J_dimer,_J_eigenvalues, J_eigenvectors);
    if (tools::globals::verbose && _bse_exc + _ct < 10) {
        CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << "Eigenvectors of J" << flush;

        CTP_LOG(ctp::logDEBUG, *_pLog) << J_eigenvectors << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << "J_eigenvalues[Hrt]" << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << _J_eigenvalues << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
    }
    //Calculate projection on subspace for every pair of excitons separately
    for (int stateA = 0; stateA < _levA; stateA++) {
        for (int stateB = 0; stateB < _levB; stateB++) {
            int stateBd = stateB + _bseA_exc;
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Calculating coupling between exciton A" << stateA + 1 << " and exciton B" << stateB + 1 << flush;
            std::vector<unsigned> index;
            std::vector<int> signvec;
            for (unsigned i = 0; i < _bse_exc + _ct; i++) {
                if (i == unsigned(stateA) || i == unsigned(stateBd)) {

                    double close = 0.0;
                    unsigned ind = 0;
                    int sign = 0;
                    //row
                    for (unsigned j = 0; j < _bse_exc + _ct; j++) {
                        bool check = true;
                        // if index i is already in index
                        // should not happen but if one vector was similar to two others.
                        for (unsigned l = 0; l < index.size(); l++) {
                            if (j == index[l]) {
                                check = false;
                                break;
                            }
                        }

                        if (check && std::abs(J_eigenvectors(i, j)) > close) {
                            ind = j;
                            close = std::abs(J_eigenvectors(i, j));
                            if (J_eigenvectors(i, j) >= 0) {
                                sign = 1;
                            } else {
                                sign = -1;
                            }
                        }
                    }
                    index.push_back(ind);
                    signvec.push_back(sign);
                }
            }

            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Order is: [Initial state n->nth eigenvalue]" << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "    A" << stateA + 1 << ":" << stateA + 1 << "->" << index[0] + 1 << " ";
            CTP_LOG(ctp::logDEBUG, *_pLog) << "    B" << stateB + 1 << ":" << stateBd + 1 << "->" << index[1] + 1 << " " << flush;

            //setting up transformation matrix _T and diagonal matrix _E for the eigenvalues;

            ub::matrix<double> _E = ub::zero_matrix<double>(2, 2);
            ub::matrix<double> _T = ub::zero_matrix<double>(2, 2);
            //find the eigenvectors which are most similar to the initial states

            //row 
            for (unsigned i = 0; i < 2; i++) {
                unsigned k = index[i];
                double sign = signvec[i];
                double normr = 1 / std::sqrt(J_eigenvectors(stateA, k) * J_eigenvectors(stateA, k) + J_eigenvectors(stateBd, k) * J_eigenvectors(stateBd, k));
                _T(0, i) = sign * J_eigenvectors(stateA, k) * normr;
                _T(1, i) = sign * J_eigenvectors(stateBd, k) * normr;
                _E(i, i) = _J_eigenvalues(k);
            }


            if ((_T(1, 1) * _T(0, 0) - _T(1, 0) * _T(0, 1)) < 0) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << " Reduced state matrix is not in a right handed basis, multiplying second eigenvector by -1 " << flush;
                _T(0, 1) = -_T(0, 1);
                _T(1, 1) = -_T(1, 1);
            }

            if (tools::globals::verbose) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << "_T" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << _T << flush;

            }

            ub::matrix<double> S_small = ub::prod(_T, ub::trans(_T));
            if (tools::globals::verbose) {

                CTP_LOG(ctp::logDEBUG, *_pLog) << "S_small" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << S_small << flush;

            }
            //orthogonalize that matrix
            double small = linalg_loewdin(_E, S_small);
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Smallest value of dimer overlapmatrix is " << small << flush;
            if (tools::globals::verbose) {

                CTP_LOG(ctp::logDEBUG, *_pLog) << "S-1/2" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << S_small << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << "E_ortho" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << _E << flush;
            }
            _T = ub::prod(_T, S_small);
            //cout <<  ub::prod(_T,ub::trans(_T))<<endl;
            if (tools::globals::verbose) {

                CTP_LOG(ctp::logDEBUG, *_pLog) << "T_ortho" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << _T << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
            }

            ub::matrix<double> temp = ub::prod(_T, _E);

            ub::matrix<double> _J_small = ub::prod(temp, ub::trans(_T));
            if (tools::globals::verbose) {

                CTP_LOG(ctp::logDEBUG, *_pLog) << "T_ortho*E_ortho" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << temp << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << "T_ortho*E_ortho*T_ortho^T" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << _J_small << flush;
            }

            _J(stateA, stateBd) = _J_small(0, 1);
            _J(stateBd, stateA) = _J_small(1, 0);

        }
    }
       
    return _J;
}


    
}}
