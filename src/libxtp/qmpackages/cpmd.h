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

#ifndef __VOTCA_XTP_CPMD_H
#define	__VOTCA_XTP_CPMD_H


#include <votca/ctp/apolarsite.h>
#include <votca/xtp/qmpackage.h>

#include <string> 



namespace votca { namespace xtp {
/**
    \brief Wrapper for the CPMD program
 
    The Cpmd class executes the CPMD package 
    and extracts information from its log and io files
    
*/
class Cpmd : public QMPackage
{
public:   

   std::string getPackageName() { return "cpmd"; }

   void Initialize( Property *options );

   /* Writes CPMD input file with coordinates of all atoms
    */
   bool WriteInputFile( std::vector< ctp::Segment* > segments, Orbitals* orbitals_guess = NULL, std::vector<ctp::PolarSeg*> PolarSegments = {});
   
   bool WriteShellScript();

   bool Run(Orbitals* _orbitals = NULL );

   void CleanUp();
   
   bool CheckLogFile();

   bool ParseLogFile( Orbitals* _orbitals );

   bool ParseOrbitalsFile( Orbitals* _orbitals ){return true;};
   
   bool setMultipoleBackground( std::vector<ctp::PolarSeg*> multipoles){ return true; };
   
      
   std::string getScratchDir( ) { return _scratch_dir; }
   
   bool loadMatrices(Orbitals * _orbitals);
   
   Cpmd(){
       _ZV=NULL;
       _NA=NULL;
       _NUMAOR=NULL;
       VOTCA2CPMD_map=NULL;
       CPMD2VOTCA_map=NULL;
   }
   
   ~Cpmd(){
       //free the memory
       if(_ZV!=NULL)
       {
           delete[] _ZV;
           delete[] _NA;
           delete[] _NUMAOR;
       }
       if(VOTCA2CPMD_map!=NULL){ delete[] VOTCA2CPMD_map;}
       if(CPMD2VOTCA_map!=NULL){ delete[] CPMD2VOTCA_map;}
   };
   
private:  

    std::string                              _shell_file_name;
    std::string                              _chk_file_name;
    std::string                              _scratch_dir;
    std::string                              _input_vxc_file_name;    
    std::string                              _cleanup;
    std::string                              _vdWfooter;
    
    
    bool _rsrt;             //have data from previous run of CPMD we want to reuse
    bool _optWF;            //optimize wavefunction
    bool _elpot;            //calculate and output electrostatic potential (needs conversion with cpmd2cube.x)
    bool _projectWF;        //project wavefunction onto localized atomic orbitals
    bool _popAnalysis;      //do population analysis, required to extract WF coefficients in atomic basis
    bool _getMat;           //get density and overlap matrices, requires _popAnalysis and _projectWF
    bool _useGrimmeVDW;     //whether to use the Grimme correction to get VDW effects
    double _pwCutoff;       //plane wave cutoff (in Ry)
    double _convCutoff;     //cutoff for MO convergence
    int _symmetry;          //symmetry number (0=isolated, 1=simple cubic, 2=FCC, 3=BCC, ... check CPMD manual under "SYMETRY")
    std::string _cell;      //cell dimensions, check CPMD manual under "CELL"
    std::string _functional;//BLYP, B3LYP, HCTH, LDE, etc.
    std::string _rsrt_kwds; //what parts to reuse from previous run
    std::string _pplib_path;//full path to the pseudopotential library of CPMD
    std::string _wfOpt_input_file_name; //name of the input file for WF optimization, only used when running CPMD twice
    std::string _wfOpt_log_file_name;
    std::string _custom_CPMD_controlls; //string to be added into the &CPMD section of the input file(s)
    
    
    std::map<std::string,std::string> _ppFileNames;   //pseudopotential file names indexed by element name
    std::map<std::string,std::string> _ppLData;       //LMAX, LOC and SKIP data for pseudopotential file
    std::map<std::string,int> _nAtomsOfElement;       //number of atoms of element indexed by element name
    list<std::string> _elements;                      //keeps track of element names and their order in CPMD
    
    int _NSP;                    
    double *_ZV;             //core charge
    int *_NA, *_NUMAOR;
    int *VOTCA2CPMD_map, *CPMD2VOTCA_map;
    std::map<int,std::string> CPMD2TYPE_map;
    
    
    ub::symmetric_matrix<double>            _overlap; //overlap matrix, from OVERLAP file
    ub::symmetric_matrix<double>            _density; //density matrix, calculated here
    ub::matrix<double>                      _coefs;   //coefficients of MOs expressed in basis set, from WFNCOEF
    

    int NumberOfElectrons( std::string _line ); 
    int BasisSetSize( std::string _line ); 
    int EnergiesFromLog( std::string _line, std::ifstream inputfile ); 
    std::string FortranFormat( const double &number );
    int NumbfQC( std::string _shell_type);
    int NumbfGW( std::string _shell_type);
    int NumbfQC_cart( std::string _shell_type);
    void WriteBasisSet(std::vector<ctp::Segment* > segments, ofstream &_com_file);
    
    int ConvAtomIndex_CPMD2VOTCA(int indx){
        return(CPMD2VOTCA_map[indx]);
    }
    
    int ConvAtomIndex_VOTCA2CPMD(int indx){
        return(VOTCA2CPMD_map[indx]);
    }

    
    
};


}}

#endif	/* __VOTCA_XTP_CPMD_H */

