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
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aoshell.h"
#include "votca/xtp/qmatom.h"
#include "votca/xtp/elements.h"
#include <votca/tools/constants.h>



namespace votca { namespace xtp {

 AOBasis::~AOBasis() {
        for (vector< AOShell* >::iterator it = _aoshells.begin(); it != _aoshells.end() ; it++ ) delete (*it);
        _aoshells.clear();
         }

AOShell* AOBasis::addShell( string shellType,int Lmax,int Lmin, double shellScale, int shellFunc, int startIndex, int offset, vec pos, string name, int index )
    {
        AOShell* aoshell = new AOShell( shellType,Lmax,Lmin, shellScale, shellFunc, startIndex, offset, pos, name, index, this );
        _aoshells.push_back(aoshell);
        return aoshell;
        }

void AOBasis::ReorderMOs(ub::matrix<double> &v, const string& start, const string& target) {

    // cout << " Reordering MOs from " << start << " to " << target << endl;

    if (start == target) {
        return;
    }
    
    if(target=="orca" || target=="nwchem"){
        vector<int> multiplier = getMultiplierVector(target,start);
        // and reorder rows of _orbitals->_mo_coefficients() accordingly
        MultiplyMOs(v, multiplier);
    }

    // get reordering vector _start -> target

    vector<int> order = getReorderVector(start, target);
    
    // Sanity check
    if (v.size2() != order.size()) {
        cerr << "Size mismatch in ReorderMOs" << v.size2() << ":" << order.size() << endl;
        throw std::runtime_error("Abort!");
    }

    // actual swapping of coefficients
    for (unsigned _i_orbital = 0; _i_orbital < v.size1(); _i_orbital++) {
        for (unsigned s = 1, d; s < order.size(); ++s) {
            for (d = order[s]; d < s; d = order[d]) {
                ;
            }
            if (d == s) while (d = order[d], d != s) swap(v(_i_orbital, s), v(_i_orbital, d));
        }
    }

    // NWChem has some strange minus in d-functions
    if (start == "nwchem" || start == "orca") {
        
        vector<int> multiplier = getMultiplierVector(start, target);
        // and reorder rows of _orbitals->_mo_coefficients() accordingly
        MultiplyMOs(v, multiplier);

    }


    return;
}

void AOBasis::ReorderMatrix(ub::symmetric_matrix<double> &v,const string& start,const string& target ){
    if (start==target){
        return;
    }
    vector<int> order = getReorderVector(start, target);
    
     if (v.size2() != order.size()) {
        cerr << "Size mismatch in ReorderMatrix" << v.size2() << ":" << order.size() << endl;
        throw std::runtime_error("Abort!");
    }

    ub::symmetric_matrix<double> temp=v;
    for(unsigned i=0;i<temp.size1();i++){
        int i_index=order[i];
        for(unsigned j=0;j<temp.size1();j++){
            int j_index=order[j];
            v(i_index,j_index)=temp(i,j);
        }
    }
    
    
    return;
}


void AOBasis::MultiplyMOs(ub::matrix<double> &v, vector<int> const &multiplier )  {
          // Sanity check
          if ( v.size2() != multiplier.size() ) {
              cerr << "Size mismatch in MultiplyMOs" << v.size2() << ":" << multiplier.size() << endl;
              throw std::runtime_error( "Abort!");
          }

          for ( unsigned _i_orbital = 0; _i_orbital < v.size1(); _i_orbital++ ){

               for ( unsigned _i_basis = 0; _i_basis < v.size2(); _i_basis++ ){

                   v(_i_orbital, _i_basis ) = multiplier[_i_basis] * v(_i_orbital, _i_basis );

               }


           }
       }
//this is for gaussian only to transform from gaussian ordering cartesian to spherical in gaussian ordering not more
ub::matrix<double> AOBasis::getTransformationCartToSpherical(const string& package){
    ub::matrix<double>_trafomatrix;
    if ( package != "gaussian" ){
        cout << " I should not have been called, will do nothing! " << endl;
    } else {
        // go through basisset, determine function sizes
        int _dim_sph = 0;
        int _dim_cart = 0;
        for (AOShellIterator _is = this->firstShell(); _is != this->lastShell() ; _is++ ) {
            const AOShell* _shell = getShell( _is );
            const string& _type =  _shell->getType();

            _dim_sph  += NumFuncShell( _type );
            _dim_cart +=NumFuncShell_cartesian( _type );

        }

     _trafomatrix= ub::zero_matrix<double>( _dim_sph , _dim_cart );
        // initialize _trafomatrix
        

        // now fill it
        int _row_start = 0;
        int _col_start = 0;
        for (AOShellIterator _is = this->firstShell(); _is != this->lastShell() ; _is++ ) {
            const AOShell* _shell = getShell( _is );
            string _type =  _shell->getType();
            int _row_end = _row_start +NumFuncShell( _type );
            int _col_end = _col_start +NumFuncShell_cartesian( _type );

            ub::matrix_range< ub::matrix<double> > _submatrix = ub::subrange( _trafomatrix, _row_start, _row_end, _col_start, _col_end);
            addTrafoCartShell(  _shell, _submatrix  );

            _row_start = _row_end;
            _col_start = _col_end;

        }
    }
    return _trafomatrix;
}

//only for gaussian package
void AOBasis::addTrafoCartShell(const AOShell* shell , ub::matrix_range< ub::matrix<double> >& _trafo ){


    // fill _local according to _lmax;
    int _lmax = shell->getLmax();
    string _type = shell->getType();

    int _sph_size =NumFuncShell( _type ) + OffsetFuncShell( _type );
    int _cart_size = NumFuncShell_cartesian( _type ) + OffsetFuncShell_cartesian( _type )  ;

    // cout << "    local size : " << _sph_size << " : " << _cart_size << endl;

    ub::matrix<double> _local =  ub::zero_matrix<double>(_sph_size,_cart_size);

    // s-functions
    _local(0,0) = 1.0; // s

    // p-functions
    if ( _lmax > 0 ){
        _local(1,1) = 1.0;
        _local(2,2) = 1.0;
        _local(3,3) = 1.0;
    }

    // d-functions
    if ( _lmax > 1 ){
        _local(4,4) = -0.5;             // d3z2-r2 (dxx)
        _local(4,5) = -0.5;             // d3z2-r2 (dyy)
        _local(4,6) =  1.0;             // d3z2-r2 (dzz)
        _local(5,8) =  1.0;             // dxz
        _local(6,9) =  1.0;             // dyz
        _local(7,4) = 0.5*sqrt(3.0);    // dx2-y2 (dxx)
        _local(7,5) = -_local(7,4);      // dx2-y2 (dyy)
        _local(8,7) = 1.0;              // dxy
     }

    if ( _lmax > 2 ){
        cerr << " Gaussian input with f- functions or higher not yet supported!" << endl;
        exit(1);
    }

    // now copy to _trafo
    for ( int _i_sph = 0 ; _i_sph < NumFuncShell( _type ) ; _i_sph++ ){
        for  ( int _i_cart = 0 ; _i_cart < NumFuncShell_cartesian( _type ) ; _i_cart++ ){


            _trafo( _i_sph , _i_cart ) = _local( _i_sph + OffsetFuncShell( _type ) , _i_cart +  OffsetFuncShell_cartesian( _type ) );

        }
    }
    return;
}


int AOBasis::getMaxFunctions () {

    int _maxfunc = 0;

    // go through basisset
    for (AOShellIterator _is = firstShell(); _is != lastShell() ; _is++ ) {
        const AOShell* _shell = this->getShell( _is );
        int _func_here = _shell->getNumFunc();
        if ( _func_here > _maxfunc ) _maxfunc = _func_here;
    }
    return _maxfunc;
}


vector<int> AOBasis::getMultiplierVector( const string& start, const string& target){
    vector<int> multiplier;
    multiplier.reserve(_AOBasisSize);
    // go through basisset
    for (AOShellIterator _is = firstShell(); _is != lastShell() ; _is++ ) {
        const AOShell* _shell = this->getShell( _is );
        addMultiplierShell(  start, target, _shell->getType(), multiplier );
    }
    return multiplier;
    }

void AOBasis::addMultiplierShell(const string& start, const string& target, const string& shell_type, vector<int>& multiplier) {


    if (target == "xtp" || start=="xtp") {
        // current length of vector
        //int _cur_pos = multiplier.size() - 1;

        // single type shells defined here
        if (shell_type.length() == 1) {
            if (shell_type == "S") {
                multiplier.push_back(1);
            }
            else if (shell_type == "P") {
                multiplier.push_back(1);
                multiplier.push_back(1);
                multiplier.push_back(1);
            }
            else if (shell_type == "D") {
                if (start == "nwchem") {
                    multiplier.push_back(-1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                } else if (start == "orca"){
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                } else {
                    cerr << "Tried to get multipliers d-functions from package " << start << ".";
                    throw std::runtime_error("Multiplication not implemented yet!");
                }
            }
            else if (shell_type == "F") {
                if ( start == "orca" ){
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(-1);
                    multiplier.push_back(-1);

                } else {
                cerr << "Tried to get multipliers d-functions from package " << start << ".";
                throw std::runtime_error("Multiplication not implemented yet!");
                }
            }
            else if (shell_type == "G") {
                cerr << "Tried to get multipliers g-functions . ";
                throw std::runtime_error("Multiplication not implemented yet!");
            }
        } else {
            // for combined shells, iterate over all contributions
            //_nbf = 0;
            for (unsigned i = 0; i < shell_type.length(); ++i) {
                string local_shell = string(shell_type, i, 1);
                addMultiplierShell(start, target, local_shell, multiplier);
            }
        }
    } else {

        cerr << "Tried to reorder functions (multiplier)  from " << start << " to " << target << endl;
        throw std::runtime_error("Reordering not implemented yet!");


    }
    return;
}


vector<int>  AOBasis::getReorderVector(const string& start,const string& target){
    vector<int> neworder;
    neworder.reserve(_AOBasisSize);
    // go through basisset
    for (AOShellIterator _is = firstShell(); _is != lastShell() ; _is++ ) {
        const AOShell* _shell = getShell( _is );
        addReorderShell( start, target, _shell->getType(), neworder );
    }
    return neworder;
}

vector<int> AOBasis::invertOrder(const vector<int>& order ){
    vector<int>neworder=vector<int>(order.size());
    for(unsigned i=0;i<order.size();i++){
        neworder[order[i]]=int(i);
    }
    
    return neworder;
}



void AOBasis::addReorderShell(const string& start,const string& target,const string& shell_type, vector<int>& neworder ) {
    vector<int> order;
    // current length of vector
    int _cur_pos = neworder.size() -1 ;

    if ( target == "xtp" || start=="xtp"){

    // single type shells defined here
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){
           order.push_back( _cur_pos + 1 );
       }//for S


       //votca order is z,y,x e.g. Y1,0 Y1,-1 Y1,1
       else if (shell_type == "P") {
                if (start == "orca") {
                    order.push_back(_cur_pos + 1);
                    order.push_back(_cur_pos + 3);
                    order.push_back(_cur_pos + 2);
                } else if (start == "gaussian" || start == "nwchem") {
                    order.push_back(_cur_pos + 3);
                    order.push_back(_cur_pos + 2);
                    order.push_back(_cur_pos + 1);
                } else if (start == "votca") {//for usage with old orb files
                    order.push_back(_cur_pos + 3);
                    order.push_back(_cur_pos + 2);
                    order.push_back(_cur_pos + 1);
                } else if (start == "xtp") {
                    order.push_back(_cur_pos + 1);
                    order.push_back(_cur_pos + 2);
                    order.push_back(_cur_pos + 3);
                }else {
               cerr << "Tried to reorder p-functions from package " << start << ".";
               throw std::runtime_error( "Reordering not implemented yet!");
           }
        } //for P
       //votca order is d3z2-r2 dyz dxz dxy dx2-y2 e.g. Y2,0 Y2,-1 Y2,1 Y2,-2 Y2,2
       else if ( shell_type == "D" ){
           //orca order is d3z2-r2 dxz dyz dx2-y2 dxy
           if ( start == "gaussian"|| start=="orca"){
               order.push_back( _cur_pos + 1 );
               order.push_back( _cur_pos + 3 );
               order.push_back( _cur_pos + 2 );
               order.push_back( _cur_pos + 5 );
               order.push_back( _cur_pos + 4 );
           } else if ( start == "nwchem") {
               // nwchem order is dxy dyz d3z2-r2 -dxz dx2-y2
               order.push_back( _cur_pos + 4  );
               order.push_back( _cur_pos + 2 );
               order.push_back( _cur_pos + 1 );
               //neworder.push_back( -(_cur_pos + 1) ); // bloody inverted sign // BUG!!!!!!!
               order.push_back( _cur_pos + 3 );
               order.push_back( _cur_pos + 5 );

           }else if ( start == "votca") { //for usage with old orb files

               order.push_back( _cur_pos + 3 );
               order.push_back( _cur_pos + 2 );
               order.push_back( _cur_pos + 4 );
               order.push_back( _cur_pos + 1 );
               order.push_back( _cur_pos + 5 );
            }else if ( start == "xtp") {

               order.push_back( _cur_pos + 1 );
               order.push_back( _cur_pos + 2 );
               order.push_back( _cur_pos + 3 );
               order.push_back( _cur_pos + 4 );
               order.push_back( _cur_pos + 5 );
            }else {
               cerr << "Tried to reorder d-functions from package " << start << ".";
               throw std::runtime_error( "Reordering not implemented yet!");
           }
       }
       else if ( shell_type == "F" ){
           if ( start == "gaussian" || start == "orca" ){
               order.push_back( _cur_pos + 1 );
               order.push_back( _cur_pos + 3 );
               order.push_back( _cur_pos + 2 );
               order.push_back( _cur_pos + 5 );
               order.push_back( _cur_pos + 4 );
               order.push_back( _cur_pos + 7 );
               order.push_back( _cur_pos + 6 );
           } else if ( start == "xtp" ){
               order.push_back( _cur_pos + 1 );
               order.push_back( _cur_pos + 2 );
               order.push_back( _cur_pos + 3 );
               order.push_back( _cur_pos + 4 );
               order.push_back( _cur_pos + 5 );
               order.push_back( _cur_pos + 6 );
               order.push_back( _cur_pos + 7 );
           } else {
               cerr << "Tried to reorder f-functions from package " << start << ".";
               throw std::runtime_error( "Reordering not implemented yet!");
           }
       }
       else{
           cerr << "Tried to reorder functions  of shell type "<<shell_type<<endl;
           throw std::runtime_error( "Reordering not implemented yet!");
       }
    } else {
        // for combined shells, iterate over all contributions
        //_nbf = 0;
        for( unsigned i = 0; i < shell_type.length(); ++i) {
           string local_shell =    string( shell_type, i, 1 );
           this->addReorderShell( start, target, local_shell, order  );
        }
    }
    
    if(start=="xtp"){
        order=invertOrder(order );  
    }
        neworder.insert(neworder.end(), order.begin(), order.end());
    
    } else {

        cerr << "Tried to reorder functions (neworder) from " << start << " to " << target << endl;
        throw std::runtime_error( "Reordering not implemented yet!");

    }
    return;
}


std::vector<AOShell*> AOBasis::getShellsperAtom(int AtomId){
  std::vector<AOShell*> result;
  for(const auto& aoshell:_aoshells){
    if(aoshell->_atomindex==AtomId){
      result.push_back(aoshell);
    }
  }
  return result;
}

int AOBasis::getFuncperAtom(int AtomId) const{
  int number=0;
  for(const auto& aoshell:_aoshells){
    if(aoshell->_atomindex==AtomId){
      number+=aoshell->_numFunc;
    }
  }
  
  return number;
}


void AOBasis::AOBasisFill(BasisSet* bs , vector<QMAtom* > _atoms, int _fragbreak  ) {
  Elements elementinfo;
  vector<QMAtom* > :: iterator ait;

 _AOBasisSize =0;
 _AOBasisFragA=0;
 _AOBasisFragB=0;

 

 // loop over atoms
 for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
    
    vec pos=(*ait)->getPos();
    
    // get element type of the atom
    string  name = (*ait)->getType();
    //assign its nuclear charge
    (*ait)->nuccharge=elementinfo.getNucCrg(name);
    // get the basis set entry for this element
    Element* element = bs->getElement(name);
    
    int funcperAtom=0;
              // and loop over all shells
    for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
              Shell* shell = (*its);
              int numfuncshell=NumFuncShell(shell->getType());
              AOShell* aoshell = addShell(shell->getType(), shell->getLmax(), shell->getLmin(), shell->getScale(),
                      numfuncshell, _AOBasisSize, OffsetFuncShell(shell->getType()), pos, name, (*ait)->getAtomID());
              _AOBasisSize += numfuncshell;
              funcperAtom+=numfuncshell;
              for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                  GaussianPrimitive* gaussian = *itg;
                  aoshell->addGaussian(gaussian->decay, gaussian->contraction);
              }
              aoshell->CalcMinDecay();     
          }
    if ( (*ait)->getAtomID() < _fragbreak ) _AOBasisFragA = _AOBasisSize;

   
}

 if ( _fragbreak < 0 ) {
     _AOBasisFragA = _AOBasisSize;
     _AOBasisFragB = 0;
 } else {
     _AOBasisFragB = _AOBasisSize - _AOBasisFragA;
 }
 return;
}





void AOBasis::ECPFill(BasisSet* bs , vector<QMAtom* > _atoms  ) {
        
        vector< QMAtom* > :: iterator ait;
       _AOBasisSize = 0;
      

       // loop over atoms
       for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
          
          // get coordinates of this atom and convert from Angstrom to Bohr
          vec pos=(*ait)->getPos();
          // get element type of the atom
          string  name = (*ait)->getType();
          // get the basis set entry for this element
          if(name=="H" || name=="He"){continue;}
            Element* element = bs->getElement(name);
            //update nuclear charge of QMAtom
           
            (*ait)->ecpcharge=element->getNcore();
          // and loop over all shells
           
          int lmax=0;
          for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
               Shell* shell = (*its);
               if(shell->getType().size()>1){
                   throw runtime_error("In ecps no combined shells e.g. SP are allowed");
               }
               int l=FindLmax(shell->getType() );
               if (its == element->firstShell()) lmax = l;
               // first shell is local component, identification my negative angular momentum

                //why is the shell not properly used, apparently it is only used in aoecp.cc and there it is iterated over so l and l make no difference CHECK!!
                   AOShell* aoshell = addShell( shell->getType(),l,l, shell->getScale(), lmax, l, l, pos, name, (*ait)->getAtomID() );
                   _AOBasisSize += NumFuncShell( shell->getType() );
                   for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                      GaussianPrimitive* gaussian = *itg;
                      aoshell->addGaussian(gaussian->power, gaussian->decay, gaussian->contraction);
               }
                   aoshell->CalcMinDecay();
                 
          }

         
      }
       return;
}











}}
