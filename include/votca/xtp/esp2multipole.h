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

#ifndef _VOTCA_XTP_ESP2MULTIPOLE_H
#define _VOTCA_XTP_ESP2MULTIPOLE_H

#include <stdio.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/mulliken.h>
#include <votca/xtp/lowdin.h>
#include <votca/xtp/nbo.h>
#include <votca/ctp/logger.h>
#include <votca/xtp/qmmachine.h>
#include <boost/filesystem.hpp>

namespace votca { namespace xtp {
    using namespace std;
    
class Esp2multipole 
{
public:

    Esp2multipole (ctp::Logger* log) {_log=log; }
   ~Esp2multipole () {};

    string Identify() { return "esp2multipole"; }

    void   Initialize(Property *options);
    
   
    void Extractingcharges( Orbitals& _orbitals );
    void WritetoFile(string _output_file,  string identifier="esp2multipole");
    string GetIdentifier();

private:
    
    int         _state_no;  
    int         _openmp_threads;
    string      _state;
    string      _method;
    string      _spin;
    string      _integrationmethod;
    string      _gridsize;
    bool        _use_mulliken;
    bool        _use_lowdin;
    bool        _use_CHELPG;
    bool        _use_bulkESP;
    bool        _use_GDMA;
    bool        _use_CHELPG_SVD;
    bool        _use_NBO;
    bool        _use_ecp;
    bool        _do_svd;
    double      _conditionnumber;
    vector< QMAtom* > _Atomlist;
    
    ctp::Logger*      _log;
    
    

};



}}


#endif