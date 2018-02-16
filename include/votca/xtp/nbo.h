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

#ifndef __XTP_NBO__H
#define	__XTP_NBO__H


#include <votca/xtp/elements.h>
#include <votca/xtp/aobasis.h>
#include <votca/ctp/logger.h>
#include <votca/xtp/qmatom.h>


/**
* \brief Takes a list of atoms, and the corresponding density and overlap matrices and puts out a table of partial charges
*
* 
* 
*/
using namespace votca::tools;


namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
class NBO{
public:
    
    NBO(ctp::Logger *log){_log = log;}
   ~NBO(){};
       
   void EvaluateNBO(std::vector< QMAtom* >& _atomlist,const ub::matrix<double> &_dmat,const AOBasis &_basis, BasisSet &bs);
  
private:
    
     ctp::Logger *_log;
     Elements _elements; 
    
    ub::matrix<double> IntercenterOrthogonalisation(ub::matrix<double> &P,ub::matrix<double> &Overlap,vector< QMAtom* >& _atomlist, BasisSet &bs);
    void TransformMatrixtoNewBasis(ub::matrix<double>& Matrix,const ub::matrix<double>& transformation);
};
}}

#endif /* NBO_H */


