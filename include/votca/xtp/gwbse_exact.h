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

#ifndef _VOTCA_XTP_GWBSE_EXACT_H
#define _VOTCA_XTP_GWBSE_EXACT_H

#include <votca/ctp/logger.h>
#include <votca/xtp/eigen.h>

namespace votca
{
namespace xtp
{

class GWBSE_Exact {
    
private:
    
    // ****** Misc. ******
    
    ctp::Logger* _pLog;
    
    // ****** Orbitals ******
    
    Orbitals& _orbitals;

    int _homo; // HOMO index
    
    // RPA
    int _rpamin;
    int _rpamax;
    
    // QP
    int _qpmin;
    int _qpmax;
    int _qptotal;

    // BSE
    int _bse_vmin;
    int _bse_vmax;
    int _bse_cmin;
    int _bse_cmax;
    int _bse_maxeigenvectors;
    
    // ****** Basis sets ******
    
    std::string _auxbasis_name;
    std::string _dftbasis_name;
    
public:

    void setLogger(ctp::Logger* pLog) {
        _pLog = pLog;
    }

    GWBSE_Exact(Orbitals& orbitals)
        : _orbitals(orbitals) {
    };

    bool Evaluate();

};

}
}

#endif /* _VOTCA_XTP_GWBSE_EXACT_H */