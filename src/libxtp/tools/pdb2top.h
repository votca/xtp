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

#ifndef _PDB2Top__H
#define _PDB2Top__H


#include <votca/xtp/topology.h>
#include <votca/xtp/atom.h>
#include <votca/xtp/logger.h>
#include <boost/algorithm/string.hpp>
#include <votca/tools/vec.h>
#include <boost/format.hpp>
#include <vector>
#include <sstream> 


using boost::format;
using boost::io::group;

namespace votca { namespace xtp {
    namespace ba = boost::algorithm;
    using namespace std;
    
class PDB2Top : public QMTool
{
public:

    PDB2Top() { };
   ~PDB2Top() { };
   
    string Identify() { return "pdb2top"; }
    // run sequence
    void   Initialize(tools::Property *options);
    bool   Evaluate();

    // helpful guys
    void readPDB();
    void readGRO();
    
    // make top file
    void top2txt();
    
    // error formating function
    void error1(string line){ cout << endl; throw runtime_error(line); }
    
private:
    string      _input_pdb;
    string      _input_gro;
    string      _out_top;
    
    bool        _has_pdb;
    bool        _has_gro;
    
    int         _numMol;
    
    Topology    _top;
};

void PDB2Top::Initialize(tools::Property* options) {
    // read options    
    string key = "options.pdb2top.";

    // set boolean constants to false
    _has_pdb = false;
    _has_gro = false;
    
    // find PDB, then GRO, then error
    if ( options->exists(key+"pdb") ){
        _input_pdb      = options->get(key+"pdb").as<string> ();
        _has_pdb        = true;
        cout << endl << "... ... PDB file: \t" << _input_pdb;
    }
    else if ( options->exists(key+"gro") ){
        _input_gro      = options->get(key+"gro").as<string> ();
        _has_gro        = true;
        cout << endl << "... ... GRO file: \t" << _input_gro;
    }
    else{
        error1( "... ... Error. Unsupported input file format. \n"
                "... ... Supported formats: pdb, gro \n");     
    }

    if ( options->exists(key+"top") ){
        _out_top = options->get(key+"top").as<string> ();
        cout << "\n" "... ... User defined top file: \t" << _numMol;
    }
    else{
        _out_top = "topol.top"; 
        cout <<"\n" "... ... Default top file: topol.top \n";
    }

    
    if ( options->exists(key+"num") ){
        _numMol = options->get(key+"num").as<int> ();
        cout << "\n" "... ... User defined num of mols: \t" << _numMol;
    }
    else{_numMol = 1; cout <<"\n" "... ... Default num of mols: 1 \n";}
}

bool PDB2Top::Evaluate() {
    if (_has_pdb){
        readPDB();
    }
    else if (_has_gro){
        readGRO();
    }
    else {error1("... ... I have no file to read from.\n"
                 "... ... tags: gro, pdb \n"); }
    
    // make actual topol.top
    top2txt();
    return true;
}

void PDB2Top::top2txt(){

    // generating part
   Topology * _topPtr = &_top;

    // preps
    stringstream ss, sbody, stype;
    
    string weight("1.000"), charge("0.000"), resref("-1"), 
           sigma ("0.00000e+00"), eps ("0.00000e+00");
    
    vector <char> atTypesLst;
    
    // iterations
    vector <Fragment*> _fragsPtr = _topPtr->Fragments();
    vector <Fragment*>::iterator _fragIt;
    for (_fragIt = _fragsPtr.begin(); _fragIt!=_fragsPtr.end(); _fragIt++){
        
        // read atoms from fragments
        // iterate
        Fragment* _frag = *(_fragIt);
        vector <Atom*> _atsPtr = _frag->Atoms();
        vector <Atom*>::iterator _atIt;
        
        // print name as comment
        sbody << "; res " << _frag->getName() << endl;
        
        for (_atIt = _atsPtr.begin(); _atIt!=_atsPtr.end(); _atIt++){
            
            Atom* _at = *(_atIt);

            // print  name, mass, charge, ptype, sigma, eps
            // if not found add, else dont
            if( find(atTypesLst.begin(), atTypesLst.end(),
                    _at->getName()[0]) == atTypesLst.end()){
            
                atTypesLst.push_back(_at->getName()[0]);
                
                stype << format("%1$6s%2$11s%3$12s%4$9s%5$13s%6$13s\n") 
                        % _at->getName()[0]     % weight
                        % charge                % "A" 
                        % sigma                 % eps           ;
            }

            // print atid, elem, resid,resname, atname, id, charge, mass
            sbody << format(    "%1$10s%2$10s%3$10s%4$10s"
                                "%5$10s%6$10s%7$10s%8$10s\n"  )
                        % _at->getId()       % _at->getName()[0]
                        % _at->getResnr()    % _at->getResname()
                        % _at->getName()     % _at->getId()
                        % charge             % weight           ;
        }
    }
    // output part
    ss      << "; "                             << '\n' 
            << ";  generated by votca top2map " << '\n'
            << ";  for xtp_map ONLY"            << '\n'
            << ";  you CAN NOT use this top "   << '\n'
            << ";  for simulations! "           << '\n'
            << "; "                             << '\n'
            << ""                               << '\n'
            << "[ defaults ]"                   << '\n'
            << format(";%1$7s%2$17s%3$16s%4$14s%5$8s\n") 
                        % "nbfunc"      % "comb-rule" 
                        % "gen-pairs"   % "fudgeLJ" 
                        % "fudgeLJ"
            << format("%1$8s%2$17s%3$16s%4$14s%5$8s\n") 
                        % '1'           % '3'
                        % "yes"         % "0.5"
                        % "0.5"
            << ""                               << '\n'
            << "[ atomtypes ]"                  << '\n'
            << format(";%1$5s%2$11s%3$12s%4$9s%5$13s%6$13s\n") 
                        % "name"     % "mass"
                        % "charge"   % "ptype" 
                        % "sigma"    % "eps" ;
    // types
    ss << stype.str();
            
    // 
    ss      << ""                               << '\n'
            << "[ moleculetype ]"               << '\n'
            << "; Name            nrexcl"       << '\n'
            << "Other               3"          << '\n'
            << ""                               << '\n'
            << "[ atoms ]"                      << '\n'
            << format(    ";%1$9s%2$10s%3$10s%4$10s"
                          "%5$10s%6$10s%7$10s%8$10s\n"  )
                        % "nr"       % "type"
                        % "resnr"    % "residue" 
                        % "atom"     % "cgnr" 
                        % "charge"   % "mass" ;
    // body
    ss << sbody.str();
    
    ss << ""                                    << '\n'
       << "[ system ]"                          << '\n'
       << "; Name"                              << '\n'
       << "Protein"                             << '\n'
       << ""                                    << '\n'
       << "[ molecules ]"                       << '\n'
       << "; Compound        #mols"             << '\n'
       << format("Other               %s") 
                                     % _numMol  << '\n' ;

    
    // print entire stream to file;
    std::ofstream _outfile( _out_top.c_str());
    _outfile << ss.str();
    _outfile.close();
}

void PDB2Top::readPDB(){
   cout << endl << "... ... Assuming: PDB";
    
    // set molecule >> segment >> fragment
    // reconnect them all
    Topology * _topPtr = 0;
    _topPtr = &_top;
    
    Molecule * _molPtr = 0;
    // direct
    _molPtr = _topPtr->AddMolecule("M1");
                // inverse
                _molPtr->setTopology(_topPtr);
    
    Segment  * _segPtr  = 0;
    // direct
    _segPtr = _topPtr->AddSegment("S1");
               _molPtr->AddSegment(_segPtr);
               // inverse
                _segPtr->setTopology(_topPtr);
                _segPtr->setMolecule(_molPtr);

    // try: read PDB file
    std::ifstream _file( _input_pdb.c_str());
    if (_file.fail()) {
        error1( "... ... Can not open: " + _input_pdb + "\n"
                "... ... Does it exist? Is it correct file name?\n");
    }
    else{
        cout << endl << 
                ("... ... File " + _input_pdb + ""
                 " was opened successfully.\n");
    }

    // read PDB line by line
    string _line;
    
    // counters for loops
//    int _atom_id = 0;
    int _newResNum = 0;

    while ( std::getline(_file, _line,'\n') ){
        if(     boost::find_first(_line, "ATOM"  )   || 
                boost::find_first(_line, "HETATM")      
                ){
            
            //      according to PDB format
            string _recType    (_line,( 1-1),6); // str,  "ATOM", "HETATM"
            string _atNum      (_line,( 7-1),6); // int,  Atom serial number
            string _atName     (_line,(13-1),4); // str,  Atom name
            string _atAltLoc   (_line,(17-1),1); // char, Alternate location indicator
            string _resName    (_line,(18-1),4); // str,  Residue name
            string _chainID    (_line,(22-1),1); // char, Chain identifier
            string _resNum     (_line,(23-1),4); // int,  Residue sequence number
            string _atICode    (_line,(27-1),1); // char, Code for insertion of res
            string _x          (_line,(31-1),8); // float 8.3 ,x
            string _y          (_line,(39-1),8); // float 8.3 ,y
            string _z          (_line,(47-1),8); // float 8.3 ,z
            string _atOccup    (_line,(55-1),6); // float  6.2, Occupancy
            string _atTFactor  (_line,(61-1),6); // float  6.2, Temperature factor
            string _segID      (_line,(73-1),4); // str, Segment identifier
            string _atElement  (_line,(77-1),2); // str, Element symbol
            string _atCharge   (_line,(79-1),2); // str, Charge on the atom

            ba::trim(_recType);
            ba::trim(_atNum);
            ba::trim(_atName);
            ba::trim(_atAltLoc);
            ba::trim(_resName);
            ba::trim(_chainID);
            ba::trim(_resNum);
            ba::trim(_atICode);
            ba::trim(_x);
            ba::trim(_y);
            ba::trim(_z);
            ba::trim(_atOccup);
            ba::trim(_atTFactor);
            ba::trim(_segID);
            ba::trim(_atElement);
            ba::trim(_atCharge);
            
            double _xd(0),_yd(0),_zd(0);
            int _resNumInt(0); 
            
            try
            {
            _xd = stod(_x);
            _yd = stod(_y);
            _zd = stod(_z);
            _resNumInt = boost::lexical_cast<int>(_resNum);
            }
            catch(boost::bad_lexical_cast &)
            {
                error1( "... ... Can not convert PDB coord line!\n"
                        "... ... Atom number: " + _atNum + "\n"
                        "... ... Make sure this line is PDB style\n");
            }
            
            tools::vec r(_xd , _yd , _zd);

            // set fragment
            // reconnect to topology, molecule, segment
            Fragment * _fragPtr = 0;
            // make new frag for new res number
            // otherwise use last created
            if ( _newResNum != _resNumInt ){

                _newResNum = _resNumInt;
                string _newResName = _resName+'_'+_resNum;
                
                // direct
                _fragPtr = _topPtr->AddFragment(_newResName);
                           _molPtr->AddFragment(_fragPtr);
                           _segPtr->AddFragment(_fragPtr);
                          // inverse
                          _fragPtr->setTopology(_topPtr);
                          _fragPtr->setMolecule(_molPtr);
                          _fragPtr->setSegment(_segPtr);        
            }
            else{
                _fragPtr = _topPtr->Fragments().back();
            }
            if (_fragPtr==0) {error1("Zero pointer in GRO reader. Why?");}
                        
            // set atom
            // reconnect to topology, molecule, segment, fragment
            Atom * _atmPtr = 0;
            // direct
            _atmPtr = _topPtr->AddAtom(_atName);
                      _molPtr->AddAtom(_atmPtr);
                      _segPtr->AddAtom(_atmPtr);
                     _fragPtr->AddAtom(_atmPtr);
                      // inverse
                      _atmPtr->setTopology(_topPtr);
                      _atmPtr->setMolecule(_molPtr);        
                      _atmPtr->setSegment(_segPtr);
                      _atmPtr->setFragment(_fragPtr);
                      
            _atmPtr->setResnr        (_resNumInt);
            _atmPtr->setResname      (_resName);
            _atmPtr->setPos          (r);
        }
    }
    
    return;
}

void PDB2Top::readGRO(){
    cout << endl << "... ... Assuming: GRO";

    // set molecule >> segment >> fragment
    // reconnect them all
    Topology * _topPtr = 0;
    _topPtr = &_top;
    
    Molecule * _molPtr = 0;
    // direct
    _molPtr = _topPtr->AddMolecule("M1");
                // inverse
                _molPtr->setTopology(_topPtr);
    
    Segment  * _segPtr  = 0;
    // direct
    _segPtr = _topPtr->AddSegment("S1");
               _molPtr->AddSegment(_segPtr);
               // inverse
                _segPtr->setTopology(_topPtr);
                _segPtr->setMolecule(_molPtr);

    // try: read GRO file
    std::ifstream _file( _input_gro.c_str());
    if (_file.fail()) {
        error1( "... ... Can not open: " + _input_gro + "\n"
                "... ... Does it exist? Is it correct file name?\n");
    }
    else{
        cout << endl << 
                ("... ... File " + _input_gro + ""
                 " was opened successfully.\n");
    }

    // read GRO line by line
    string _line;
    
    // counters for loops
    int _newResNum = -1; // res reference
    int _atTotl = 0;  // total num of atoms in GRO
    int _atCntr = 0;  // atom number counter
    
    // GRO: first two lines are tech specs -> ignore them
    // ignore first line, it's a comment
    std::getline(_file, _line,'\n');

    // GRO check: if second line can cast to int, then ok

    try
    {   
        // first line, number of atoms in XYZ
        std::getline(_file, _line,'\n');
        ba::trim(_line);
        _atTotl = boost::lexical_cast<int>(_line);
    }
    catch(boost::bad_lexical_cast &)
    {
        error1( "... ... Bad GRO file format!\n"
                "... ... First two line must contain technical specs.\n");
    }

    // actual loop
    while ( std::getline(_file, _line,'\n') ){
        if (_atCntr < _atTotl){
            
            string _resNum     (_line, 0,5); // int,  Residue number
            string _resName    (_line, 5,5); // str,  Residue name
            string _atName     (_line,10,5); // str,  Atom name
            string _atNum      (_line,15,5); // int,  Atom number
            string _x          (_line,20,8); // float 8.3 ,x
            string _y          (_line,28,8); // float 8.3 ,y
            string _z          (_line,36,8); // float 8.3 ,z
            
            ba::trim(_atNum);
            ba::trim(_atName);
            ba::trim(_resNum);
            ba::trim(_resName);
            ba::trim(_x);
            ba::trim(_y);
            ba::trim(_z);
            
            // try cast
            int _resNumInt(0);//,_atNumInt(0);
            double _xd(0),_yd(0),_zd(0);
            try
            {
                _resNumInt = boost::lexical_cast<int>(_resNum);
                //_atNumInt  = boost::lexical_cast<int>(_atNum);

                _xd = stod(_x);
                _yd = stod(_y);
                _zd = stod(_z);
            }
            catch (boost::bad_lexical_cast &)
            {
                error1( "... ... Can not convert GRO coord line!\n"
                        "... ... Atom number: " + _atNum + "\n"
                        "... ... Make sure this line is GRO style\n");
            }
            
            tools::vec r(_xd , _yd , _zd);
                
            // set fragment
            // reconnect to topology, molecule, segment
            Fragment * _fragPtr = 0;
            // make new frag for new res number
            // otherwise use last created
            if ( _newResNum != _resNumInt ){

                _newResNum = _resNumInt;
                string _newResName = _resName+'_'+_resNum;
                
                // direct
                _fragPtr = _topPtr->AddFragment(_newResName);
                           _molPtr->AddFragment(_fragPtr);
                           _segPtr->AddFragment(_fragPtr);
                          // inverse
                          _fragPtr->setTopology(_topPtr);
                          _fragPtr->setMolecule(_molPtr);
                          _fragPtr->setSegment(_segPtr);        
            }
            else{
                _fragPtr = _topPtr->Fragments().back();
            }
            if (_fragPtr==0) {error1("Zero pointer in GRO reader. Why?");}
                        
            // set atom
            // reconnect to topology, molecule, segment, fragment
            Atom * _atmPtr = 0;
            // direct
            _atmPtr = _topPtr->AddAtom(_atName);
                      _molPtr->AddAtom(_atmPtr);
                      _segPtr->AddAtom(_atmPtr);
                     _fragPtr->AddAtom(_atmPtr);
                      // inverse
                      _atmPtr->setTopology(_topPtr);
                      _atmPtr->setMolecule(_molPtr);        
                      _atmPtr->setSegment(_segPtr);
                      _atmPtr->setFragment(_fragPtr);
        
            _atmPtr->setResnr        (_resNumInt);
            _atmPtr->setResname      (_resName);
            _atmPtr->setPos          (r);
        
        }
        _atCntr++;
    }
    
    return;
}

}}


#endif
