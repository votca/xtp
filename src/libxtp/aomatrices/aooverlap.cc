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

#include <votca/xtp/aomatrix.h>

#include <votca/xtp/aobasis.h>

#include <vector>



using namespace votca::tools;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        void AOOverlap::FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix, const AOShell* _shell_row, const AOShell* _shell_col) {


            // shell info, only lmax tells how far to go
            int _lmax_row = _shell_row->getLmax();
            int _lmax_col = _shell_col->getLmax();

            // set size of internal block for recursion
            int _nrows = this->getBlockSize(_lmax_row);
            int _ncols = this->getBlockSize(_lmax_col);

            if (_lmax_col > 6 || _lmax_row > 6) {
                cerr << "Orbitals higher than i are not yet implemented. This should not have happened!" << flush;
                exit(1);
            }

            /* FOR CONTRACTED FUNCTIONS, ADD LOOP OVER ALL DECAYS IN CONTRACTION
             * MULTIPLY THE TRANSFORMATION MATRICES BY APPROPRIATE CONTRACTION 
             * COEFFICIENTS, AND ADD TO matrix(i,j)
             */

            // get shell positions
            const vec& _pos_row = _shell_row->getPos();
            const vec& _pos_col = _shell_col->getPos();
            const vec _diff = _pos_row - _pos_col;
            std::vector<double> _pma(3, 0.0);
            std::vector<double> _pmb(3, 0.0);

            double _distsq = (_diff * _diff);
            int n_orbitals[] = {1, 4, 10, 20, 35, 56, 84};

            int nx[] = {0,
                1, 0, 0,
                2, 1, 1, 0, 0, 0,
                3, 2, 2, 1, 1, 1, 0, 0, 0, 0,
                4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0,
                5, 4, 4, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0};

            int ny[] = {0,
                0, 1, 0,
                0, 1, 0, 2, 1, 0,
                0, 1, 0, 2, 1, 0, 3, 2, 1, 0,
                0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0,
                0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
                0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 6, 5, 4, 3, 2, 1, 0};

            int nz[] = {0,
                0, 0, 1,
                0, 0, 1, 0, 1, 2,
                0, 0, 1, 0, 1, 2, 0, 1, 2, 3,
                0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4,
                0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5,
                0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6};


            int i_less_x[] = {0,
                0, 0, 0,
                1, 2, 3, 0, 0, 0,
                4, 5, 6, 7, 8, 9, 0, 0, 0, 0,
                10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 0, 0, 0, 0, 0,
                20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 0, 0, 0, 0, 0, 0,
                35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 0, 0, 0, 0, 0, 0, 0};

            int i_less_y[] = {0,
                0, 0, 0,
                0, 1, 0, 2, 3, 0,
                0, 4, 0, 5, 6, 0, 7, 8, 9, 0,
                0, 10, 0, 11, 12, 0, 13, 14, 15, 0, 16, 17, 18, 19, 0,
                0, 20, 0, 21, 22, 0, 23, 24, 25, 0, 26, 27, 28, 29, 0, 30, 31, 32, 33, 34, 0,
                0, 35, 0, 36, 37, 0, 38, 39, 40, 0, 41, 42, 43, 44, 0, 45, 46, 47, 48, 49, 0, 50, 51, 52, 53, 54, 55, 0};

            int i_less_z[] = {0,
                0, 0, 0,
                0, 0, 1, 0, 2, 3,
                0, 0, 4, 0, 5, 6, 0, 7, 8, 9,
                0, 0, 10, 0, 11, 12, 0, 13, 14, 15, 0, 16, 17, 18, 19,
                0, 0, 20, 0, 21, 22, 0, 23, 24, 25, 0, 26, 27, 28, 29, 0, 30, 31, 32, 33, 34,
                0, 0, 35, 0, 36, 37, 0, 38, 39, 40, 0, 41, 42, 43, 44, 0, 45, 46, 47, 48, 49, 0, 50, 51, 52, 53, 54, 55};




            // iterate over Gaussians in this _shell_row
            for (AOShell::GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr) {
                // iterate over Gaussians in this _shell_col
                const double _decay_row = itr->getDecay();

                for (AOShell::GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc) {


                    // get decay constants 
                    const double _decay_col = itc->getDecay();

                    // some helpers
                    const double _fak = 0.5 / (_decay_row + _decay_col);
                    const double _fak2 = 2.0 * _fak;

                    // check if distance between postions is big, then skip step   
                    double _exparg = _fak2 * _decay_row * _decay_col *_distsq;
                    if (_exparg > 30.0) {
                        continue;
                    }
                    // initialize local matrix block for unnormalized cartesians
                    ub::matrix<double> _ol = ub::zero_matrix<double>(_nrows, _ncols);


                    // Definition of coefficients for recursive overlap formulas 
                    // A for rows (i). B for columns (j)
                    const double PmA0 = _fak2 * (_decay_row * _pos_row.getX() + _decay_col * _pos_col.getX()) - _pos_row.getX();
                    const double PmA1 = _fak2 * (_decay_row * _pos_row.getY() + _decay_col * _pos_col.getY()) - _pos_row.getY();
                    const double PmA2 = _fak2 * (_decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ()) - _pos_row.getZ();

                    const double PmB0 = _fak2 * (_decay_row * _pos_row.getX() + _decay_col * _pos_col.getX()) - _pos_col.getX();
                    const double PmB1 = _fak2 * (_decay_row * _pos_row.getY() + _decay_col * _pos_col.getY()) - _pos_col.getY();
                    const double PmB2 = _fak2 * (_decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ()) - _pos_col.getZ();



                    // calculate matrix elements
                    _ol(0, 0) = pow(4.0 * _decay_row*_decay_col, 0.75) * pow(_fak2, 1.5) * exp(-_exparg); // s-s element
                    //cout << "\t setting s-s: " << _ol(0,0) << endl;


                    //Integrals     p - s
                    if (_lmax_row > 0) {
                        _ol(Cart::x, 0) = PmA0 * _ol(0, 0);
                        _ol(Cart::y, 0) = PmA1 * _ol(0, 0);
                        _ol(Cart::z, 0) = PmA2 * _ol(0, 0);
                    }
                    //------------------------------------------------------

                    //Integrals     d - s
                    if (_lmax_row > 1) {
                        double term = _fak * _ol(0, 0);
                        _ol(Cart::xx, 0) = PmA0 * _ol(Cart::x, 0) + term;
                        _ol(Cart::xy, 0) = PmA0 * _ol(Cart::y, 0);
                        _ol(Cart::xz, 0) = PmA0 * _ol(Cart::z, 0);
                        _ol(Cart::yy, 0) = PmA1 * _ol(Cart::y, 0) + term;
                        _ol(Cart::yz, 0) = PmA1 * _ol(Cart::z, 0);
                        _ol(Cart::zz, 0) = PmA2 * _ol(Cart::z, 0) + term;
                    }
                    //------------------------------------------------------

                    //Integrals     f - s
                    if (_lmax_row > 2) {
                        _ol(Cart::xxx, 0) = PmA0 * _ol(Cart::xx, 0) + 2 * _fak * _ol(Cart::x, 0);
                        _ol(Cart::xxy, 0) = PmA1 * _ol(Cart::xx, 0);
                        _ol(Cart::xxz, 0) = PmA2 * _ol(Cart::xx, 0);
                        _ol(Cart::xyy, 0) = PmA0 * _ol(Cart::yy, 0);
                        _ol(Cart::xyz, 0) = PmA0 * _ol(Cart::yz, 0);
                        _ol(Cart::xzz, 0) = PmA0 * _ol(Cart::zz, 0);
                        _ol(Cart::yyy, 0) = PmA1 * _ol(Cart::yy, 0) + 2 * _fak * _ol(Cart::y, 0);
                        _ol(Cart::yyz, 0) = PmA2 * _ol(Cart::yy, 0);
                        _ol(Cart::yzz, 0) = PmA1 * _ol(Cart::zz, 0);
                        _ol(Cart::zzz, 0) = PmA2 * _ol(Cart::zz, 0) + 2 * _fak * _ol(Cart::z, 0);
                    }
                    //------------------------------------------------------

                    //Integrals     g - s
                    if (_lmax_row > 3) {
                        double term_xx = _fak * _ol(Cart::xx, 0);
                        double term_yy = _fak * _ol(Cart::yy, 0);
                        double term_zz = _fak * _ol(Cart::zz, 0);
                        _ol(Cart::xxxx, 0) = PmA0 * _ol(Cart::xxx, 0) + 3 * term_xx;
                        _ol(Cart::xxxy, 0) = PmA1 * _ol(Cart::xxx, 0);
                        _ol(Cart::xxxz, 0) = PmA2 * _ol(Cart::xxx, 0);
                        _ol(Cart::xxyy, 0) = PmA0 * _ol(Cart::xyy, 0) + term_yy;
                        _ol(Cart::xxyz, 0) = PmA1 * _ol(Cart::xxz, 0);
                        _ol(Cart::xxzz, 0) = PmA0 * _ol(Cart::xzz, 0) + term_zz;
                        _ol(Cart::xyyy, 0) = PmA0 * _ol(Cart::yyy, 0);
                        _ol(Cart::xyyz, 0) = PmA0 * _ol(Cart::yyz, 0);
                        _ol(Cart::xyzz, 0) = PmA0 * _ol(Cart::yzz, 0);
                        _ol(Cart::xzzz, 0) = PmA0 * _ol(Cart::zzz, 0);
                        _ol(Cart::yyyy, 0) = PmA1 * _ol(Cart::yyy, 0) + 3 * term_yy;
                        _ol(Cart::yyyz, 0) = PmA2 * _ol(Cart::yyy, 0);
                        _ol(Cart::yyzz, 0) = PmA1 * _ol(Cart::yzz, 0) + term_zz;
                        _ol(Cart::yzzz, 0) = PmA1 * _ol(Cart::zzz, 0);
                        _ol(Cart::zzzz, 0) = PmA2 * _ol(Cart::zzz, 0) + 3 * term_zz;
                    }
                    //------------------------------------------------------

                    //Integrals     h - s
                    if (_lmax_row > 4) {
                        double term_xxx = _fak * _ol(Cart::xxx, 0);
                        double term_yyy = _fak * _ol(Cart::yyy, 0);
                        double term_zzz = _fak * _ol(Cart::zzz, 0);
                        _ol(Cart::xxxxx, 0) = PmA0 * _ol(Cart::xxxx, 0) + 4 * term_xxx;
                        _ol(Cart::xxxxy, 0) = PmA1 * _ol(Cart::xxxx, 0);
                        _ol(Cart::xxxxz, 0) = PmA2 * _ol(Cart::xxxx, 0);
                        _ol(Cart::xxxyy, 0) = PmA1 * _ol(Cart::xxxy, 0) + term_xxx;
                        _ol(Cart::xxxyz, 0) = PmA1 * _ol(Cart::xxxz, 0);
                        _ol(Cart::xxxzz, 0) = PmA2 * _ol(Cart::xxxz, 0) + term_xxx;
                        _ol(Cart::xxyyy, 0) = PmA0 * _ol(Cart::xyyy, 0) + term_yyy;
                        _ol(Cart::xxyyz, 0) = PmA2 * _ol(Cart::xxyy, 0);
                        _ol(Cart::xxyzz, 0) = PmA1 * _ol(Cart::xxzz, 0);
                        _ol(Cart::xxzzz, 0) = PmA0 * _ol(Cart::xzzz, 0) + term_zzz;
                        _ol(Cart::xyyyy, 0) = PmA0 * _ol(Cart::yyyy, 0);
                        _ol(Cart::xyyyz, 0) = PmA0 * _ol(Cart::yyyz, 0);
                        _ol(Cart::xyyzz, 0) = PmA0 * _ol(Cart::yyzz, 0);
                        _ol(Cart::xyzzz, 0) = PmA0 * _ol(Cart::yzzz, 0);
                        _ol(Cart::xzzzz, 0) = PmA0 * _ol(Cart::zzzz, 0);
                        _ol(Cart::yyyyy, 0) = PmA1 * _ol(Cart::yyyy, 0) + 4 * term_yyy;
                        _ol(Cart::yyyyz, 0) = PmA2 * _ol(Cart::yyyy, 0);
                        _ol(Cart::yyyzz, 0) = PmA2 * _ol(Cart::yyyz, 0) + term_yyy;
                        _ol(Cart::yyzzz, 0) = PmA1 * _ol(Cart::yzzz, 0) + term_zzz;
                        _ol(Cart::yzzzz, 0) = PmA1 * _ol(Cart::zzzz, 0);
                        _ol(Cart::zzzzz, 0) = PmA2 * _ol(Cart::zzzz, 0) + 4 * term_zzz;
                    }
                    //------------------------------------------------------

                    //Integrals     i - s
                    if (_lmax_row > 5) {
                        double term_xxxx = _fak * _ol(Cart::xxxx, 0);
                        double term_xyyy = _fak * _ol(Cart::xyyy, 0);
                        double term_xzzz = _fak * _ol(Cart::xzzz, 0);
                        double term_yyyy = _fak * _ol(Cart::yyyy, 0);
                        double term_yyzz = _fak * _ol(Cart::yyzz, 0);
                        double term_yzzz = _fak * _ol(Cart::yzzz, 0);
                        double term_zzzz = _fak * _ol(Cart::zzzz, 0);
                        _ol(Cart::xxxxxx, 0) = PmA0 * _ol(Cart::xxxxx, 0) + 5 * term_xxxx;
                        _ol(Cart::xxxxxy, 0) = PmA1 * _ol(Cart::xxxxx, 0);
                        _ol(Cart::xxxxxz, 0) = PmA2 * _ol(Cart::xxxxx, 0);
                        _ol(Cart::xxxxyy, 0) = PmA1 * _ol(Cart::xxxxy, 0) + term_xxxx;
                        _ol(Cart::xxxxyz, 0) = PmA1 * _ol(Cart::xxxxz, 0);
                        _ol(Cart::xxxxzz, 0) = PmA2 * _ol(Cart::xxxxz, 0) + term_xxxx;
                        _ol(Cart::xxxyyy, 0) = PmA0 * _ol(Cart::xxyyy, 0) + 2 * term_xyyy;
                        _ol(Cart::xxxyyz, 0) = PmA2 * _ol(Cart::xxxyy, 0);
                        _ol(Cart::xxxyzz, 0) = PmA1 * _ol(Cart::xxxzz, 0);
                        _ol(Cart::xxxzzz, 0) = PmA0 * _ol(Cart::xxzzz, 0) + 2 * term_xzzz;
                        _ol(Cart::xxyyyy, 0) = PmA0 * _ol(Cart::xyyyy, 0) + term_yyyy;
                        _ol(Cart::xxyyyz, 0) = PmA2 * _ol(Cart::xxyyy, 0);
                        _ol(Cart::xxyyzz, 0) = PmA0 * _ol(Cart::xyyzz, 0) + term_yyzz;
                        _ol(Cart::xxyzzz, 0) = PmA1 * _ol(Cart::xxzzz, 0);
                        _ol(Cart::xxzzzz, 0) = PmA0 * _ol(Cart::xzzzz, 0) + term_zzzz;
                        _ol(Cart::xyyyyy, 0) = PmA0 * _ol(Cart::yyyyy, 0);
                        _ol(Cart::xyyyyz, 0) = PmA0 * _ol(Cart::yyyyz, 0);
                        _ol(Cart::xyyyzz, 0) = PmA0 * _ol(Cart::yyyzz, 0);
                        _ol(Cart::xyyzzz, 0) = PmA0 * _ol(Cart::yyzzz, 0);
                        _ol(Cart::xyzzzz, 0) = PmA0 * _ol(Cart::yzzzz, 0);
                        _ol(Cart::xzzzzz, 0) = PmA0 * _ol(Cart::zzzzz, 0);
                        _ol(Cart::yyyyyy, 0) = PmA1 * _ol(Cart::yyyyy, 0) + 5 * term_yyyy;
                        _ol(Cart::yyyyyz, 0) = PmA2 * _ol(Cart::yyyyy, 0);
                        _ol(Cart::yyyyzz, 0) = PmA2 * _ol(Cart::yyyyz, 0) + term_yyyy;
                        _ol(Cart::yyyzzz, 0) = PmA1 * _ol(Cart::yyzzz, 0) + 2 * term_yzzz;
                        _ol(Cart::yyzzzz, 0) = PmA1 * _ol(Cart::yzzzz, 0) + term_zzzz;
                        _ol(Cart::yzzzzz, 0) = PmA1 * _ol(Cart::zzzzz, 0);
                        _ol(Cart::zzzzzz, 0) = PmA2 * _ol(Cart::zzzzz, 0) + 5 * term_zzzz;
                    }
                    //------------------------------------------------------


                    if (_lmax_col > 0) {

                        //Integrals     s - p
                        _ol(0, Cart::x) = PmB0 * _ol(0, 0);
                        _ol(0, Cart::y) = PmB1 * _ol(0, 0);
                        _ol(0, Cart::z) = PmB2 * _ol(0, 0);
                        //------------------------------------------------------

                        //Integrals     p - p     d - p     f - p     g - p     h - p     i - p
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            _ol(_i, Cart::x) = PmB0 * _ol(_i, 0) + nx[_i] * _fak * _ol(i_less_x[_i], 0);
                            _ol(_i, Cart::y) = PmB1 * _ol(_i, 0) + ny[_i] * _fak * _ol(i_less_y[_i], 0);
                            _ol(_i, Cart::z) = PmB2 * _ol(_i, 0) + nz[_i] * _fak * _ol(i_less_z[_i], 0);
                        }
                        //------------------------------------------------------

                    } // end if (_lmax_col > 0)


                    if (_lmax_col > 1) {

                        //Integrals     s - d
                        double term = _fak * _ol(0, 0);
                        _ol(0, Cart::xx) = PmB0 * _ol(0, Cart::x) + term;
                        _ol(0, Cart::xy) = PmB0 * _ol(0, Cart::y);
                        _ol(0, Cart::xz) = PmB0 * _ol(0, Cart::z);
                        _ol(0, Cart::yy) = PmB1 * _ol(0, Cart::y) + term;
                        _ol(0, Cart::yz) = PmB1 * _ol(0, Cart::z);
                        _ol(0, Cart::zz) = PmB2 * _ol(0, Cart::z) + term;
                        //------------------------------------------------------

                        //Integrals     p - d     d - d     f - d     g - d     h - d     i - d
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            double term = _fak * _ol(_i, 0);
                            _ol(_i, Cart::xx) = PmB0 * _ol(_i, Cart::x) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::x) + term;
                            _ol(_i, Cart::xy) = PmB0 * _ol(_i, Cart::y) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::y);
                            _ol(_i, Cart::xz) = PmB0 * _ol(_i, Cart::z) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::z);
                            _ol(_i, Cart::yy) = PmB1 * _ol(_i, Cart::y) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::y) + term;
                            _ol(_i, Cart::yz) = PmB1 * _ol(_i, Cart::z) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::z);
                            _ol(_i, Cart::zz) = PmB2 * _ol(_i, Cart::z) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::z) + term;
                        }
                        //------------------------------------------------------

                    } // end if (_lmax_col > 1)


                    if (_lmax_col > 2) {

                        //Integrals     s - f
                        _ol(0, Cart::xxx) = PmB0 * _ol(0, Cart::xx) + 2 * _fak * _ol(0, Cart::x);
                        _ol(0, Cart::xxy) = PmB1 * _ol(0, Cart::xx);
                        _ol(0, Cart::xxz) = PmB2 * _ol(0, Cart::xx);
                        _ol(0, Cart::xyy) = PmB0 * _ol(0, Cart::yy);
                        _ol(0, Cart::xyz) = PmB0 * _ol(0, Cart::yz);
                        _ol(0, Cart::xzz) = PmB0 * _ol(0, Cart::zz);
                        _ol(0, Cart::yyy) = PmB1 * _ol(0, Cart::yy) + 2 * _fak * _ol(0, Cart::y);
                        _ol(0, Cart::yyz) = PmB2 * _ol(0, Cart::yy);
                        _ol(0, Cart::yzz) = PmB1 * _ol(0, Cart::zz);
                        _ol(0, Cart::zzz) = PmB2 * _ol(0, Cart::zz) + 2 * _fak * _ol(0, Cart::z);
                        //------------------------------------------------------

                        //Integrals     p - f     d - f     f - f     g - f     h - f     i - f
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            double term_x = 2 * _fak * _ol(_i, Cart::x);
                            double term_y = 2 * _fak * _ol(_i, Cart::y);
                            double term_z = 2 * _fak * _ol(_i, Cart::z);
                            _ol(_i, Cart::xxx) = PmB0 * _ol(_i, Cart::xx) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xx) + term_x;
                            _ol(_i, Cart::xxy) = PmB1 * _ol(_i, Cart::xx) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xx);
                            _ol(_i, Cart::xxz) = PmB2 * _ol(_i, Cart::xx) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::xx);
                            _ol(_i, Cart::xyy) = PmB0 * _ol(_i, Cart::yy) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yy);
                            _ol(_i, Cart::xyz) = PmB0 * _ol(_i, Cart::yz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yz);
                            _ol(_i, Cart::xzz) = PmB0 * _ol(_i, Cart::zz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::zz);
                            _ol(_i, Cart::yyy) = PmB1 * _ol(_i, Cart::yy) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::yy) + term_y;
                            _ol(_i, Cart::yyz) = PmB2 * _ol(_i, Cart::yy) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::yy);
                            _ol(_i, Cart::yzz) = PmB1 * _ol(_i, Cart::zz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::zz);
                            _ol(_i, Cart::zzz) = PmB2 * _ol(_i, Cart::zz) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::zz) + term_z;
                        }
                        //------------------------------------------------------

                    } // end if (_lmax_col > 2)


                    if (_lmax_col > 3) {

                        //Integrals     s - g
                        double term_xx = _fak * _ol(0, Cart::xx);
                        double term_yy = _fak * _ol(0, Cart::yy);
                        double term_zz = _fak * _ol(0, Cart::zz);
                        _ol(0, Cart::xxxx) = PmB0 * _ol(0, Cart::xxx) + 3 * term_xx;
                        _ol(0, Cart::xxxy) = PmB1 * _ol(0, Cart::xxx);
                        _ol(0, Cart::xxxz) = PmB2 * _ol(0, Cart::xxx);
                        _ol(0, Cart::xxyy) = PmB0 * _ol(0, Cart::xyy) + term_yy;
                        _ol(0, Cart::xxyz) = PmB1 * _ol(0, Cart::xxz);
                        _ol(0, Cart::xxzz) = PmB0 * _ol(0, Cart::xzz) + term_zz;
                        _ol(0, Cart::xyyy) = PmB0 * _ol(0, Cart::yyy);
                        _ol(0, Cart::xyyz) = PmB0 * _ol(0, Cart::yyz);
                        _ol(0, Cart::xyzz) = PmB0 * _ol(0, Cart::yzz);
                        _ol(0, Cart::xzzz) = PmB0 * _ol(0, Cart::zzz);
                        _ol(0, Cart::yyyy) = PmB1 * _ol(0, Cart::yyy) + 3 * term_yy;
                        _ol(0, Cart::yyyz) = PmB2 * _ol(0, Cart::yyy);
                        _ol(0, Cart::yyzz) = PmB1 * _ol(0, Cart::yzz) + term_zz;
                        _ol(0, Cart::yzzz) = PmB1 * _ol(0, Cart::zzz);
                        _ol(0, Cart::zzzz) = PmB2 * _ol(0, Cart::zzz) + 3 * term_zz;
                        //------------------------------------------------------

                        //Integrals     p - g     d - g     f - g     g - g     h - g     i - g
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            double term_xx = _fak * _ol(_i, Cart::xx);
                            double term_yy = _fak * _ol(_i, Cart::yy);
                            double term_zz = _fak * _ol(_i, Cart::zz);
                            _ol(_i, Cart::xxxx) = PmB0 * _ol(_i, Cart::xxx) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xxx) + 3 * term_xx;
                            _ol(_i, Cart::xxxy) = PmB1 * _ol(_i, Cart::xxx) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xxx);
                            _ol(_i, Cart::xxxz) = PmB2 * _ol(_i, Cart::xxx) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::xxx);
                            _ol(_i, Cart::xxyy) = PmB0 * _ol(_i, Cart::xyy) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xyy) + term_yy;
                            _ol(_i, Cart::xxyz) = PmB1 * _ol(_i, Cart::xxz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xxz);
                            _ol(_i, Cart::xxzz) = PmB0 * _ol(_i, Cart::xzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xzz) + term_zz;
                            _ol(_i, Cart::xyyy) = PmB0 * _ol(_i, Cart::yyy) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yyy);
                            _ol(_i, Cart::xyyz) = PmB0 * _ol(_i, Cart::yyz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yyz);
                            _ol(_i, Cart::xyzz) = PmB0 * _ol(_i, Cart::yzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yzz);
                            _ol(_i, Cart::xzzz) = PmB0 * _ol(_i, Cart::zzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::zzz);
                            _ol(_i, Cart::yyyy) = PmB1 * _ol(_i, Cart::yyy) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::yyy) + 3 * term_yy;
                            _ol(_i, Cart::yyyz) = PmB2 * _ol(_i, Cart::yyy) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::yyy);
                            _ol(_i, Cart::yyzz) = PmB1 * _ol(_i, Cart::yzz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::yzz) + term_zz;
                            _ol(_i, Cart::yzzz) = PmB1 * _ol(_i, Cart::zzz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::zzz);
                            _ol(_i, Cart::zzzz) = PmB2 * _ol(_i, Cart::zzz) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::zzz) + 3 * term_zz;
                        }
                        //------------------------------------------------------

                    } // end if (_lmax_col > 3)


                    if (_lmax_col > 4) {

                        //Integrals     s - h
                        double term_xxx = _fak * _ol(0, Cart::xxx);
                        double term_yyy = _fak * _ol(0, Cart::yyy);
                        double term_zzz = _fak * _ol(0, Cart::zzz);
                        _ol(0, Cart::xxxxx) = PmB0 * _ol(0, Cart::xxxx) + 4 * term_xxx;
                        _ol(0, Cart::xxxxy) = PmB1 * _ol(0, Cart::xxxx);
                        _ol(0, Cart::xxxxz) = PmB2 * _ol(0, Cart::xxxx);
                        _ol(0, Cart::xxxyy) = PmB1 * _ol(0, Cart::xxxy) + term_xxx;
                        _ol(0, Cart::xxxyz) = PmB1 * _ol(0, Cart::xxxz);
                        _ol(0, Cart::xxxzz) = PmB2 * _ol(0, Cart::xxxz) + term_xxx;
                        _ol(0, Cart::xxyyy) = PmB0 * _ol(0, Cart::xyyy) + term_yyy;
                        _ol(0, Cart::xxyyz) = PmB2 * _ol(0, Cart::xxyy);
                        _ol(0, Cart::xxyzz) = PmB1 * _ol(0, Cart::xxzz);
                        _ol(0, Cart::xxzzz) = PmB0 * _ol(0, Cart::xzzz) + term_zzz;
                        _ol(0, Cart::xyyyy) = PmB0 * _ol(0, Cart::yyyy);
                        _ol(0, Cart::xyyyz) = PmB0 * _ol(0, Cart::yyyz);
                        _ol(0, Cart::xyyzz) = PmB0 * _ol(0, Cart::yyzz);
                        _ol(0, Cart::xyzzz) = PmB0 * _ol(0, Cart::yzzz);
                        _ol(0, Cart::xzzzz) = PmB0 * _ol(0, Cart::zzzz);
                        _ol(0, Cart::yyyyy) = PmB1 * _ol(0, Cart::yyyy) + 4 * term_yyy;
                        _ol(0, Cart::yyyyz) = PmB2 * _ol(0, Cart::yyyy);
                        _ol(0, Cart::yyyzz) = PmB2 * _ol(0, Cart::yyyz) + term_yyy;
                        _ol(0, Cart::yyzzz) = PmB1 * _ol(0, Cart::yzzz) + term_zzz;
                        _ol(0, Cart::yzzzz) = PmB1 * _ol(0, Cart::zzzz);
                        _ol(0, Cart::zzzzz) = PmB2 * _ol(0, Cart::zzzz) + 4 * term_zzz;
                        //------------------------------------------------------

                        //Integrals     p - h     d - h     f - h     g - h     h - h     i - h
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            double term_xxx = _fak * _ol(_i, Cart::xxx);
                            double term_yyy = _fak * _ol(_i, Cart::yyy);
                            double term_zzz = _fak * _ol(_i, Cart::zzz);
                            _ol(_i, Cart::xxxxx) = PmB0 * _ol(_i, Cart::xxxx) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xxxx) + 4 * term_xxx;
                            _ol(_i, Cart::xxxxy) = PmB1 * _ol(_i, Cart::xxxx) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xxxx);
                            _ol(_i, Cart::xxxxz) = PmB2 * _ol(_i, Cart::xxxx) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::xxxx);
                            _ol(_i, Cart::xxxyy) = PmB1 * _ol(_i, Cart::xxxy) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xxxy) + term_xxx;
                            _ol(_i, Cart::xxxyz) = PmB1 * _ol(_i, Cart::xxxz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xxxz);
                            _ol(_i, Cart::xxxzz) = PmB2 * _ol(_i, Cart::xxxz) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::xxxz) + term_xxx;
                            _ol(_i, Cart::xxyyy) = PmB0 * _ol(_i, Cart::xyyy) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xyyy) + term_yyy;
                            _ol(_i, Cart::xxyyz) = PmB2 * _ol(_i, Cart::xxyy) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::xxyy);
                            _ol(_i, Cart::xxyzz) = PmB1 * _ol(_i, Cart::xxzz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xxzz);
                            _ol(_i, Cart::xxzzz) = PmB0 * _ol(_i, Cart::xzzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xzzz) + term_zzz;
                            _ol(_i, Cart::xyyyy) = PmB0 * _ol(_i, Cart::yyyy) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yyyy);
                            _ol(_i, Cart::xyyyz) = PmB0 * _ol(_i, Cart::yyyz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yyyz);
                            _ol(_i, Cart::xyyzz) = PmB0 * _ol(_i, Cart::yyzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yyzz);
                            _ol(_i, Cart::xyzzz) = PmB0 * _ol(_i, Cart::yzzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yzzz);
                            _ol(_i, Cart::xzzzz) = PmB0 * _ol(_i, Cart::zzzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::zzzz);
                            _ol(_i, Cart::yyyyy) = PmB1 * _ol(_i, Cart::yyyy) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::yyyy) + 4 * term_yyy;
                            _ol(_i, Cart::yyyyz) = PmB2 * _ol(_i, Cart::yyyy) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::yyyy);
                            _ol(_i, Cart::yyyzz) = PmB2 * _ol(_i, Cart::yyyz) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::yyyz) + term_yyy;
                            _ol(_i, Cart::yyzzz) = PmB1 * _ol(_i, Cart::yzzz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::yzzz) + term_zzz;
                            _ol(_i, Cart::yzzzz) = PmB1 * _ol(_i, Cart::zzzz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::zzzz);
                            _ol(_i, Cart::zzzzz) = PmB2 * _ol(_i, Cart::zzzz) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::zzzz) + 4 * term_zzz;
                        }
                        //------------------------------------------------------

                    } // end if (_lmax_col > 4)


                    if (_lmax_col > 5) {

                        //Integrals     s - i
                        double term_xxxx = _fak * _ol(0, Cart::xxxx);
                        double term_xyyy = _fak * _ol(0, Cart::xyyy);
                        double term_xzzz = _fak * _ol(0, Cart::xzzz);
                        double term_yyyy = _fak * _ol(0, Cart::yyyy);
                        double term_yyzz = _fak * _ol(0, Cart::yyzz);
                        double term_yzzz = _fak * _ol(0, Cart::yzzz);
                        double term_zzzz = _fak * _ol(0, Cart::zzzz);
                        _ol(0, Cart::xxxxxx) = PmB0 * _ol(0, Cart::xxxxx) + 5 * term_xxxx;
                        _ol(0, Cart::xxxxxy) = PmB1 * _ol(0, Cart::xxxxx);
                        _ol(0, Cart::xxxxxz) = PmB2 * _ol(0, Cart::xxxxx);
                        _ol(0, Cart::xxxxyy) = PmB1 * _ol(0, Cart::xxxxy) + term_xxxx;
                        _ol(0, Cart::xxxxyz) = PmB1 * _ol(0, Cart::xxxxz);
                        _ol(0, Cart::xxxxzz) = PmB2 * _ol(0, Cart::xxxxz) + term_xxxx;
                        _ol(0, Cart::xxxyyy) = PmB0 * _ol(0, Cart::xxyyy) + 2 * term_xyyy;
                        _ol(0, Cart::xxxyyz) = PmB2 * _ol(0, Cart::xxxyy);
                        _ol(0, Cart::xxxyzz) = PmB1 * _ol(0, Cart::xxxzz);
                        _ol(0, Cart::xxxzzz) = PmB0 * _ol(0, Cart::xxzzz) + 2 * term_xzzz;
                        _ol(0, Cart::xxyyyy) = PmB0 * _ol(0, Cart::xyyyy) + term_yyyy;
                        _ol(0, Cart::xxyyyz) = PmB2 * _ol(0, Cart::xxyyy);
                        _ol(0, Cart::xxyyzz) = PmB0 * _ol(0, Cart::xyyzz) + term_yyzz;
                        _ol(0, Cart::xxyzzz) = PmB1 * _ol(0, Cart::xxzzz);
                        _ol(0, Cart::xxzzzz) = PmB0 * _ol(0, Cart::xzzzz) + term_zzzz;
                        _ol(0, Cart::xyyyyy) = PmB0 * _ol(0, Cart::yyyyy);
                        _ol(0, Cart::xyyyyz) = PmB0 * _ol(0, Cart::yyyyz);
                        _ol(0, Cart::xyyyzz) = PmB0 * _ol(0, Cart::yyyzz);
                        _ol(0, Cart::xyyzzz) = PmB0 * _ol(0, Cart::yyzzz);
                        _ol(0, Cart::xyzzzz) = PmB0 * _ol(0, Cart::yzzzz);
                        _ol(0, Cart::xzzzzz) = PmB0 * _ol(0, Cart::zzzzz);
                        _ol(0, Cart::yyyyyy) = PmB1 * _ol(0, Cart::yyyyy) + 5 * term_yyyy;
                        _ol(0, Cart::yyyyyz) = PmB2 * _ol(0, Cart::yyyyy);
                        _ol(0, Cart::yyyyzz) = PmB2 * _ol(0, Cart::yyyyz) + term_yyyy;
                        _ol(0, Cart::yyyzzz) = PmB1 * _ol(0, Cart::yyzzz) + 2 * term_yzzz;
                        _ol(0, Cart::yyzzzz) = PmB1 * _ol(0, Cart::yzzzz) + term_zzzz;
                        _ol(0, Cart::yzzzzz) = PmB1 * _ol(0, Cart::zzzzz);
                        _ol(0, Cart::zzzzzz) = PmB2 * _ol(0, Cart::zzzzz) + 5 * term_zzzz;
                        //------------------------------------------------------

                        //Integrals     p - i     d - i     f - i     g - i     h - i     i - i
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            double term_xxxx = _fak * _ol(_i, Cart::xxxx);
                            double term_xyyy = _fak * _ol(_i, Cart::xyyy);
                            double term_xzzz = _fak * _ol(_i, Cart::xzzz);
                            double term_yyyy = _fak * _ol(_i, Cart::yyyy);
                            double term_yyzz = _fak * _ol(_i, Cart::yyzz);
                            double term_yzzz = _fak * _ol(_i, Cart::yzzz);
                            double term_zzzz = _fak * _ol(_i, Cart::zzzz);
                            _ol(_i, Cart::xxxxxx) = PmB0 * _ol(_i, Cart::xxxxx) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xxxxx) + 5 * term_xxxx;
                            _ol(_i, Cart::xxxxxy) = PmB1 * _ol(_i, Cart::xxxxx) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xxxxx);
                            _ol(_i, Cart::xxxxxz) = PmB2 * _ol(_i, Cart::xxxxx) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::xxxxx);
                            _ol(_i, Cart::xxxxyy) = PmB1 * _ol(_i, Cart::xxxxy) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xxxxy) + term_xxxx;
                            _ol(_i, Cart::xxxxyz) = PmB1 * _ol(_i, Cart::xxxxz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xxxxz);
                            _ol(_i, Cart::xxxxzz) = PmB2 * _ol(_i, Cart::xxxxz) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::xxxxz) + term_xxxx;
                            _ol(_i, Cart::xxxyyy) = PmB0 * _ol(_i, Cart::xxyyy) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xxyyy) + 2 * term_xyyy;
                            _ol(_i, Cart::xxxyyz) = PmB2 * _ol(_i, Cart::xxxyy) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::xxxyy);
                            _ol(_i, Cart::xxxyzz) = PmB1 * _ol(_i, Cart::xxxzz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xxxzz);
                            _ol(_i, Cart::xxxzzz) = PmB0 * _ol(_i, Cart::xxzzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xxzzz) + 2 * term_xzzz;
                            _ol(_i, Cart::xxyyyy) = PmB0 * _ol(_i, Cart::xyyyy) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xyyyy) + term_yyyy;
                            _ol(_i, Cart::xxyyyz) = PmB2 * _ol(_i, Cart::xxyyy) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::xxyyy);
                            _ol(_i, Cart::xxyyzz) = PmB0 * _ol(_i, Cart::xyyzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xyyzz) + term_yyzz;
                            _ol(_i, Cart::xxyzzz) = PmB1 * _ol(_i, Cart::xxzzz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::xxzzz);
                            _ol(_i, Cart::xxzzzz) = PmB0 * _ol(_i, Cart::xzzzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::xzzzz) + term_zzzz;
                            _ol(_i, Cart::xyyyyy) = PmB0 * _ol(_i, Cart::yyyyy) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yyyyy);
                            _ol(_i, Cart::xyyyyz) = PmB0 * _ol(_i, Cart::yyyyz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yyyyz);
                            _ol(_i, Cart::xyyyzz) = PmB0 * _ol(_i, Cart::yyyzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yyyzz);
                            _ol(_i, Cart::xyyzzz) = PmB0 * _ol(_i, Cart::yyzzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yyzzz);
                            _ol(_i, Cart::xyzzzz) = PmB0 * _ol(_i, Cart::yzzzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::yzzzz);
                            _ol(_i, Cart::xzzzzz) = PmB0 * _ol(_i, Cart::zzzzz) + nx[_i] * _fak * _ol(i_less_x[_i], Cart::zzzzz);
                            _ol(_i, Cart::yyyyyy) = PmB1 * _ol(_i, Cart::yyyyy) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::yyyyy) + 5 * term_yyyy;
                            _ol(_i, Cart::yyyyyz) = PmB2 * _ol(_i, Cart::yyyyy) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::yyyyy);
                            _ol(_i, Cart::yyyyzz) = PmB2 * _ol(_i, Cart::yyyyz) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::yyyyz) + term_yyyy;
                            _ol(_i, Cart::yyyzzz) = PmB1 * _ol(_i, Cart::yyzzz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::yyzzz) + 2 * term_yzzz;
                            _ol(_i, Cart::yyzzzz) = PmB1 * _ol(_i, Cart::yzzzz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::yzzzz) + term_zzzz;
                            _ol(_i, Cart::yzzzzz) = PmB1 * _ol(_i, Cart::zzzzz) + ny[_i] * _fak * _ol(i_less_y[_i], Cart::zzzzz);
                            _ol(_i, Cart::zzzzzz) = PmB2 * _ol(_i, Cart::zzzzz) + nz[_i] * _fak * _ol(i_less_z[_i], Cart::zzzzz) + 5 * term_zzzz;
                        }
                        //------------------------------------------------------




                    } // end if (_lmax_col > 5)



                    //cout << "Done with unnormalized matrix " << endl;

                    ub::matrix<double> _trafo_row = getTrafo(*itr);
                    ub::matrix<double> _trafo_col_tposed = ub::trans(getTrafo(*itc));

                    // cartesian -> spherical

                    ub::matrix<double> _ol_tmp = ub::prod(_trafo_row, _ol);
                    ub::matrix<double> _ol_sph = ub::prod(_ol_tmp, _trafo_col_tposed);
                    // save to _matrix
                    for (unsigned i = 0; i < _matrix.size1(); i++) {
                        for (unsigned j = 0; j < _matrix.size2(); j++) {
                            _matrix(i, j) += _ol_sph(i + _shell_row->getOffset(), j + _shell_col->getOffset());
                        }
                    }


                    _ol.clear();
                } // _shell_col Gaussians
            } // _shell_row Gaussians
        }




    }
}
