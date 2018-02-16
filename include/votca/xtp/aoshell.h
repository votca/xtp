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

#ifndef __XTP_AOSHELL__H
#define	__XTP_AOSHELL__H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <votca/xtp/basisset.h>
#include <votca/tools/constants.h>


using namespace votca::tools;

namespace votca { namespace xtp {
namespace ub = boost::numeric::ublas;
class AOBasis;
class AOShell;  




// Gaussian function: contraction*exp(-decay*r^2)
class AOGaussianPrimitive 
{
    friend class AOShell;
public:
    
    

    double getPowfactor()const {return powfactor;}
    int    getPower()const{return power;}
    double getDecay()const {return decay;}
    const std::vector<double>& getContraction()const {return contraction;}
    const AOShell* getShell() const{return aoshell;}
private:
     
    int power; // used in pseudopotenials only
    double decay;
    std::vector<double> contraction;
    AOShell* aoshell;
    double powfactor;//used in evalspace to speed up DFT
    // private constructor, only a shell can create a primitive
    AOGaussianPrimitive( double _decay, std::vector<double> _contraction, AOShell *_aoshell = NULL ) 
    : power(-1),decay(_decay),
            contraction(_contraction),
            aoshell(_aoshell) {powfactor=pow(2.0 * decay / boost::math::constants::pi<double>(), 0.75) ; }

    AOGaussianPrimitive( int _power, double _decay, std::vector<double> _contraction, AOShell *_aoshell = NULL ) 
    : power(_power),
    decay(_decay),
    contraction(_contraction),
    aoshell(_aoshell) {powfactor=pow(2.0 * decay / boost::math::constants::pi<double>(), 0.75) ; }
};      
    
/*
 * S, P, or D functions in a Gaussian-basis expansion
 */
class AOShell 
{
    //friend class AOElement;
    friend class AOBasis;
public:

    const std::string& getType() const{ return _type; }
    int    getNumFunc() const{ return _numFunc ;}
    int    getStartIndex() const{ return _startIndex ;}
    int    getOffset() const{ return _offset ;}
    int    getIndex() const{ return _atomindex;}
    const std::string& getName() const{ return _atomname;}
    
    int getLmax(  ) const{ return _Lmax;}
    int getLmin(  ) const{ return _Lmin;}
    
    const vec& getPos() const{ return _pos; }
    double getScale() const{ return _scale; }
    
    int getSize() const{ return _gaussians.size(); }
    
    void CalcMinDecay(){
     _mindecay=std::numeric_limits<double>::max();
     for(auto& gaussian:_gaussians){
         if(gaussian.getDecay()<_mindecay){
             _mindecay=gaussian.getDecay();
         }
     }
     return;
    }
    
    double getMinDecay() const{return _mindecay;}
    
    void EvalAOspace(ub::matrix_range<ub::matrix<double> >& AOvalues, const vec& grid_pos ) const;
    void EvalAOspace(ub::matrix_range<ub::matrix<double> >& AOvalues,ub::matrix_range<ub::matrix<double> >& AODervalues, const vec& grid_pos ) const;

    // iterator over pairs (decay constant; contraction coefficient)
    typedef std::vector< AOGaussianPrimitive >::const_iterator GaussianIterator;
    GaussianIterator firstGaussian() const{ return _gaussians.begin(); }
    GaussianIterator lastGaussian()const{ return _gaussians.end(); }
   
    // adds a Gaussian 
    void  addGaussian( double decay, std::vector<double> contraction ) 
    {
        AOGaussianPrimitive gaussian = AOGaussianPrimitive(decay, contraction, this);
        _gaussians.push_back( gaussian );
        return;
    }
    // used for ecps
    void  addGaussian( int power, double decay, std::vector<double> contraction ) 
    {                                                                                
        AOGaussianPrimitive gaussian = AOGaussianPrimitive(power, decay, contraction, this); 
        _gaussians.push_back( gaussian );                                                  
        return;                                                                 
    }                                                                                


    
private:   

    // only class aobasis can construct shells    
    AOShell( string type,int Lmax,int Lmin, double scale, int numFunc, int startIndex, 
            int offset, vec pos, string atomname, int atomindex, AOBasis* aobasis = NULL )
            : _type(type),_Lmax(Lmax),_Lmin(Lmin), _scale(scale), _numFunc(numFunc),
                    _startIndex(startIndex), _offset(offset), _pos(pos) , 
                    _atomname(atomname), _atomindex(atomindex) { ; }
    
    // only class aobasis can destruct shells
            ~AOShell(){};
    
    // shell type (S, P, D))
    string _type;
    int _Lmax;
    int _Lmin;
    // scaling factor
    double _scale;
    // number of functions in shell
    int _numFunc;
    double _mindecay;
    int _startIndex;
    int _offset;
    vec _pos;
    string _atomname;
    int _atomindex;
     

    
    // vector of pairs of decay constants and contraction coefficients
    std::vector< AOGaussianPrimitive > _gaussians;
    
};

    
}}

#endif	/* AOSHELL_H */

