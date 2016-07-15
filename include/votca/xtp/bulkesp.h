/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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

#ifndef BULKESP_H
#define BULKESP_H

#include <votca/xtp/espfit.h>

using namespace votca::tools;

namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
class Bulkesp: public Espfit{

public:

    struct Bond{
        QMAtom* a;
        QMAtom* b;
    };
    
    struct Molecule{
        std::vector< QMAtom* > atoms;
        std::vector< Bond > bonds;
    };
    
    
public:
    Bulkesp(Logger *log):Espfit(log){
    }
    
    std::vector<Bulkesp::Molecule> BreakIntoMolecules(std::vector< QMAtom* > a, double scale);
    
    ub::vector<double> ComputeESP(std::vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat, AOBasis &_basis,BasisSet &bs,string gridsize);
    
    void Evaluate(std::vector< QMAtom* >& _atomlist, ub::matrix<double> _MO_Coefficients, AOBasis &_basis,BasisSet &bs,string gridsize, double maxBondScale);
    
    void FillElement2NBF(std::vector< QMAtom* >& _atomlist, BasisSet &bs);
    
private:

    std::map<std::string,int> _element2NBF; //Number of Basis Functions for each element in the basis set
    list<std::string> _elements;             //list of all elements in the QMatoms vector
    
    std::map<QMAtom*,int> MapAtom2MOCoefIndex(std::vector< QMAtom* >& _atomlist);
};    
  
    
}}




#endif /* BULKESP_H */

