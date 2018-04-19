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
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>

#include <votca/xtp/aomatrix.h>

#include <votca/xtp/aobasis.h>
#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <votca/tools/linalg.h>
#include <votca/tools/elements.h>






namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    using namespace votca::tools;


    void AOECP::FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix, const AOShell* _shell_row, const AOShell* _shell_col) {

        /*
         *
         *  Currently supported:
         *
         *        S, P, D, F, G  functions in DFT basis and non-local ECPs with l = 0, 1, 2, 3
         *
         *    or
         *
         *        S, P, D, F  functions in DFT basis and non-local ECPs with l = 0, 1, 2, 3, 4
         *
         */



            // get shell positions       
    
            int _lmax_row = _shell_row->getLmax();
            std::vector<double> _contractions_row_full((_lmax_row + 1)*(_lmax_row + 1));


            int _lmax_col = _shell_col->getLmax();
            std::vector<double> _contractions_col_full((_lmax_col + 1)*(_lmax_col + 1));



            const vec& _pos_row = _shell_row->getPos();
            const vec& _pos_col = _shell_col->getPos();
            const vec _diff = _pos_row - _pos_col;
            // initialize some helper
            double _distsq = _diff*_diff;


            // iterate over Gaussians in this _shell_row
            for (AOShell::GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr) {
                // iterate over Gaussians in this _shell_col
                // get decay constant
                const double _decay_row = itr->getDecay();

                //if ( _decay_row > 0.08 ) continue;

                const std::vector<double>& _contractions_row = itr->getContraction();
                // shitty magic
                for (int L = 0; L <= _lmax_row; L++) {
                    for (int M = L*L; M < (L + 1)*(L + 1); M++) {
                        _contractions_row_full[M] = _contractions_row[L];
                    }
                }


                for (AOShell::GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc) {
                    //get decay constant
                    const double _decay_col = itc->getDecay();
                    const double _fak  = 0.5 / (_decay_row + _decay_col);
                    const double _fak2 = 2.0 * _fak;
                
                
                    double _exparg = _fak2 * _decay_row * _decay_col *_distsq;
        
                    // check if distance between postions is big, then skip step   
       
                    if ( _exparg > 30.0 ) {
                        continue;
                    }

                    const std::vector<double>& _contractions_col = itc->getContraction();
                    for (int L = 0; L <= _lmax_col; L++) {
                        for (int M = L * L; M < (L + 1)*(L + 1); M++) {
                            _contractions_col_full[M] = _contractions_col[L];
                        }
                    }
                    // for each atom and its pseudopotential, get a matrix
                    int _atomidx = 0;

                    ub::matrix<int> _power_matrix = ub::zero_matrix<int>(4, 5); // 4 fit components, non-local ECPs l = 0, 1, 2, 3, 4
                    ub::matrix<double> _decay_matrix = ub::zero_matrix<double>(4, 5);
                    ub::matrix<double> _coef_matrix  = ub::zero_matrix<double>(4, 5);

                    AOBasis::AOShellIterator final_iter = _ecp->lastShell();
                    --final_iter;
                    vec _ecp_eval_pos = vec(0.0);
                    int _lmax_ecp = 0;
                    for (AOBasis::AOShellIterator _ecpit = _ecp->firstShell(); _ecpit != _ecp->lastShell(); ++_ecpit) {

                        const AOShell* _shell_ecp = _ecp->getShell(_ecpit);
                        const vec& _ecp_pos = _shell_ecp->getPos();

                        int this_atom = _shell_ecp->getIndex();

                        const int _ecp_l = _shell_ecp->getOffset(); //  angular momentum l is stored in offset for ECP

                        // only do the non-local parts
                        if (_ecp_l < _shell_ecp->getNumFunc()) {
                            int _lmax_ecp_old = _lmax_ecp;
                            _lmax_ecp = _shell_ecp->getNumFunc() - 1;
                            int i_fit = -1;
                            for (AOShell::GaussianIterator itecp = _shell_ecp->firstGaussian(); itecp != _shell_ecp->lastGaussian(); ++itecp) {
                                i_fit++;

                                // get info for this angular momentum shell
                                const int _power_ecp = itecp->getPower();
                                const double _decay_ecp = itecp->getDecay();
                                const double _contraction_ecp = itecp->getContraction()[0];

                                // collect atom ECP
                                if (this_atom == _atomidx) {
                                    _ecp_eval_pos = _ecp_pos;
                                    _power_matrix(i_fit, _ecp_l ) = _power_ecp;
                                    _decay_matrix(i_fit, _ecp_l ) = _decay_ecp;
                                    _coef_matrix(i_fit, _ecp_l )  = _contraction_ecp;
                                }
                                
                                if ((this_atom != _atomidx ) || ( _ecpit == final_iter)) {
                                    if (this_atom != _atomidx) {
                                       _lmax_ecp = _lmax_ecp_old;
                                    }

                                    // evaluate collected data, returns a (10x10) matrix of already normalized matrix elements
                                    ub::matrix<double> VNL_ECP = calcVNLmatrix(_lmax_ecp, _ecp_eval_pos, *itr, *itc,  _power_matrix ,_decay_matrix, _coef_matrix);

                                    // consider contractions
                                    // cut out block that is needed. sum
                                    //                                   cout << "_matrix.size1,2()   " << _matrix.size1() << "    " << _matrix.size2() << endl;
                                    for ( unsigned i = 0; i < _matrix.size1(); i++ ) {
                                        for (unsigned j = 0; j < _matrix.size2(); j++) {
                                            _matrix(i,j) += VNL_ECP(i+_shell_row->getOffset(),j+_shell_col->getOffset()) * _contractions_row_full[i+_shell_row->getOffset()]* _contractions_col_full[j+_shell_col->getOffset()];
                                        }
                                    }


                                    // reset atom ECP containers
                                    _power_matrix = ub::zero_matrix<int>(4, 5); // 4 fit components, non-local ECPs l = 0, 1, 2, 3, 4
                                    _decay_matrix = ub::zero_matrix<double>(4, 5);
                                    _coef_matrix  = ub::zero_matrix<double>(4, 5);
                                    _atomidx++;
                                    i_fit = 0;
                                    //cout << "setting new matrix " << i_fit << " l " << _ecp_l << " alpha  " << _decay_ecp <<  " pref " << _contraction_ecp << endl;
                                    _power_matrix(i_fit, _ecp_l ) = _power_ecp;
                                    _decay_matrix(i_fit, _ecp_l ) = _decay_ecp;
                                    _coef_matrix(i_fit, _ecp_l )  = _contraction_ecp;
                                } // evaluate if new atom is found

                            } // all Gaussians in ecp_shell
                        } // only for non local parts

                    } // all ecp_shells
                  
                }// _shell_col Gaussians
            }// _shell_row Gaussians
         
            return;
        }

        ub::matrix<double> AOECP::calcVNLmatrix(int _lmax_ecp, const vec& posC, const AOGaussianPrimitive& _g_row, const AOGaussianPrimitive& _g_col,const ub::matrix<int>& _power_ecp, const ub::matrix<double>& _gamma_ecp,const ub::matrix<double>& _pref_ecp) {

            /* calculate the contribution of the nonlocal 
             *     ECP of atom at posC with 
             *       decay constants in _gamma_ecp
             *       coefficients in    _pref_ecp
             *       with angular momentum of max 4
             * 
             * to DFT basis shell pair 
             *    with decay alpha at posA
             *         decay beta  at posB
             */

            const double conv = 1.e-9;  // 1.e-8
            const int NMAX = 41;
            const double PI = boost::math::constants::pi<double>();
            double SQPI = sqrt(PI);
            double SQ2 = sqrt(2.);
            double SQ3 = sqrt(3.);
            double SQ5 = sqrt(5.);
            double SQ7 = sqrt(7.);

            double alpha = _g_row.getDecay();
            double beta = _g_col.getDecay();
            const vec& posA = _g_row.getShell()->getPos();
            const vec& posB = _g_col.getShell()->getPos();
            int _lmax_row = _g_row.getShell()->getLmax();
            int _lmax_col = _g_col.getShell()->getLmax();
            int _lmin = std::min({_lmax_row, _lmax_col, _lmax_ecp});
            int _lmax = std::max({_lmax_row, _lmax_col, _lmax_ecp});
            int _nsph_row = (_lmax_row + 1) * (_lmax_row + 1);
            int _nsph_col = (_lmax_col + 1) * (_lmax_col + 1);

            vec AVS = posA - posC;
            vec BVS = posB - posC;
            double AVS2 = AVS * AVS;
            double BVS2 = BVS * BVS;     

            int INULL = 0;
            if (AVS2 > 0.01) INULL = 2;
            if (BVS2 > 0.01) INULL++;


            ub::matrix<double> matrix = ub::zero_matrix<double>(_nsph_row,_nsph_col);
            const int nnonsep = _gamma_ecp.size1();
            int nmax;
            if (INULL == 0) {
                nmax = 2 * _lmin;
            } else if (INULL == 3) {
                nmax = 2 * NMAX + _lmax_row + _lmax_col;
            } else {
                nmax = NMAX + 2 * _lmax;
            }
            ub::matrix<double> XI = ub::zero_matrix<double>(_lmax_ecp + 1, nmax + 1);

            double f_even_r0 = .5 * SQPI;
            double f_even_r1 = .5;
            double f_even_r2 = .25 * SQPI;

            for (int N = 0; N <= nmax; N++) {

                if ( (N % 2) == 0 ) { // N even,   XI(L, odd N) never needed

                    if (N > 0) {
                       f_even_r0 = f_even_r2; // !!
                       f_even_r1 = f_even_r1 * double(N / 2);
                       f_even_r2 = .5 * f_even_r0 * double(N + 1);
                    }

                    double DFAK_r0 = .5 * double(N + 1);
                    double DFAK_r1 = .5 * double(N + 2);
                    double DFAK_r2 = .5 * double(N + 3);

                    for (int L = 0; L <= _lmax_ecp; L++) {

                        for (int I = 0; I < nnonsep; I++) {
                        int power = _power_ecp(I, L);
                        double DLI = (alpha + beta + _gamma_ecp(I, L));
                            if (power == 2) {
                                XI(L, N) += f_even_r2 * _pref_ecp(I, L) / pow(DLI, DFAK_r2); // r^2 terms
                            } else if (power == 0) {
                                XI(L, N) += f_even_r0 * _pref_ecp(I, L) / pow(DLI, DFAK_r0); // r^0 terms
                            } else if (power == 1) {
                                XI(L, N) += f_even_r1 * _pref_ecp(I, L) / pow(DLI, DFAK_r1); // r^1 terms
                            }
                        }

                    }

                } // end if ( (N % 2) == 0 )

            }






            /**** PREPARATIONS DONE, NOW START ******/

            // some limit determinations


            double G1 = 1.;
            double AVSSQ = 0.;
            int NMAX1 = 0;
            if (AVS2 > 0.01) {

                G1 = exp(-alpha * AVS2);
                AVSSQ = sqrt(AVS2);
                double AMAX = 0.0;  
                double fak = 2.0 * alpha * AVSSQ;
                double Pow = 1.;
                double factorialNN = 1;
                for (int NN = 0; NN <= NMAX; NN++ ) {

                    if (NN != 0) {
                        Pow = Pow * fak;
                        factorialNN = factorialNN * NN;
                    }
                    double AF = G1 * Pow / factorialNN;
                    if ((NN % 2) == 0) {
                        int ii = NN + 2 * _lmax;

                        switch (_lmax) {
                            case 0:
                                AMAX = AF * XI(0, ii);
                                break;
                            case 1:
                                AMAX = AF * std::max({XI(0, ii), XI(1, ii)});
                                break;
                            case 2:
                                AMAX = AF * std::max({XI(0, ii), XI(1, ii), XI(2, ii)});
                                break;
                            case 3:
                                AMAX = AF * std::max({XI(0, ii), XI(1, ii), XI(2, ii), XI(3, ii)});
                                break;
                            case 4:
                                AMAX = AF * std::max({XI(0, ii), XI(1, ii), XI(2, ii), XI(3, ii), XI(4, ii)});
                                break;
                        }

                        if (NMAX1 == 0 && AMAX <= conv) NMAX1 = NN;
                        if (NMAX1 != 0 && AMAX >  conv ) NMAX1 = 0;
                    }

                }
                if (NMAX1 == 0 && AMAX > conv ) NMAX1 = NMAX;

            }

            // same story for B
            double G2 = 1.;
            double BVSSQ = 0.;
            int NMAX2 = 0;
            if (BVS2 > 0.01) {

                G2 = exp(-beta * BVS2);
                BVSSQ = sqrt(BVS2);
                double BMAX = 0.0;  
                double fak = 2.0 * beta * BVSSQ;
                double Pow = 1.;
                double factorialNN = 1;
                for (int NN = 0; NN <= NMAX; NN++) {

                    if(NN != 0) {
                       Pow = Pow * fak;
                       factorialNN = factorialNN * NN;
                    }
                    double BF = G2 * Pow / factorialNN;
                    if ((NN % 2) == 0) {
                        int ii = NN + 2 * _lmax;

                        switch (_lmax) {
                            case 0:
                                BMAX = BF * XI(0, ii);
                                break;
                            case 1:
                                BMAX = BF * std::max({XI(0, ii), XI(1, ii)});
                                break;
                            case 2:
                                BMAX = BF * std::max({XI(0, ii), XI(1, ii), XI(2, ii)});
                                break;
                            case 3:
                                BMAX = BF * std::max({XI(0, ii), XI(1, ii), XI(2, ii), XI(3, ii)});
                                break;
                            case 4:
                                BMAX = BF * std::max({XI(0, ii), XI(1, ii), XI(2, ii), XI(3, ii), XI(4, ii)});
                                break;
                        }

                        if (NMAX2 == 0 && BMAX <= conv) NMAX2 = NN;
                        if (NMAX2 != 0 && BMAX >  conv) NMAX2 = 0;
                    }

                }
                if (NMAX2 == 0 && BMAX > conv ) NMAX2 = NMAX;

            }

            double GAUSS = G1 * G2;







            /****** ORIGINAL CKO SUBROUTINE **********/
            // get a multi dimensional array
            typedef boost::multi_array<double, 4> ma_type;
            typedef boost::multi_array_types::extent_range range;
            typedef ma_type::index index;
            ma_type::extent_gen extents;
            ma_type COEF;

            if (INULL != 0) {

                COEF.resize(extents[range(0, 5)][range(0, 5)][range(0, 9)][range(0, NMAX + 1)]);

                int _lmin_dft_ecp=0;
                int _lmax_dft_ecp=0;
                if (INULL == 2) {
                    _lmin_dft_ecp = std::min(_lmax_row, _lmax_ecp);
                    _lmax_dft_ecp = std::max(_lmax_row, _lmax_ecp);
                } else if (INULL == 1) {
                    _lmin_dft_ecp = std::min(_lmax_col, _lmax_ecp);
                    _lmax_dft_ecp = std::max(_lmax_col, _lmax_ecp);
                } else if (INULL == 3) {
                    int _lmax_dft = std::max(_lmax_row, _lmax_col);
                    _lmin_dft_ecp = std::min(_lmax_dft, _lmax_ecp);
                    _lmax_dft_ecp = std::max(_lmax_dft, _lmax_ecp);
                }

                for ( index i4 = 0; i4 <= NMAX; i4++ ) {
                /********** ORIGINAL CKOEF SUBROUTINE *************************/
                    int NU = i4 % 2;
                    int NG = (i4 + 1) % 2;
                    double FN1 = double(i4 + 1);
                    double FN2 = double(i4 + 2);
                    double FN3 = double(i4 + 3);
                    double FN4 = double(i4 + 4);
                    double FN5 = double(i4 + 5);
                    double FN6 = double(i4 + 6);
                    double FN7 = double(i4 + 7);
                    double FN8 = double(i4 + 8);

                    COEF[0][0][4][i4] = NG/FN1;  //  M0                      Mn is a modified spherical Bessel function of the first kind

                    if (_lmax_dft_ecp > 0) {

                        double COEFF = NU/FN2*SQ3;  //  SQ(3) * M1
                        COEF[0][1][4][i4] = COEFF;
                        COEF[1][0][4][i4] = COEFF;

                        if (_lmin_dft_ecp > 0) {
                            COEF[1][1][4][i4] = NG*3.0/FN3;        //  M0 + 2 * M2
                            COEFF = 3.0/2.0*NG*(1.0/FN1-1.0/FN3);  //  M0 - M2
                            COEF[1][1][3][i4] = COEFF;
                            COEF[1][1][5][i4] = COEFF;
                        }

                    }

                    if (_lmax_dft_ecp > 1) {

                        double COEFF = NG/2.0*SQ5*(3.0/FN3-1.0/FN1);  //  SQ(5) * M2
                        COEF[0][2][4][i4] = COEFF;
                        COEF[2][0][4][i4] = COEFF;

                        if (_lmin_dft_ecp > 0) {
                            COEFF = SQ3*SQ5/2.0*NU*(3.0/FN4-1.0/FN2);  //  ( SQ(15)/5 ) * ( 2 * M1 + 3 * M3 )
                            COEF[1][2][4][i4] = COEFF;
                            COEF[2][1][4][i4] = COEFF;
                            COEFF = 3.*SQ5/2.0*NU*(1.0/FN2-1.0/FN4);   //  ( (3*SQ(5))/5 ) * ( M1 - M3 )
                            COEF[1][2][3][i4] = COEFF;
                            COEF[1][2][5][i4] = COEFF;
                            COEF[2][1][3][i4] = COEFF;
                            COEF[2][1][5][i4] = COEFF;
                        }

                        if (_lmin_dft_ecp > 1) {
                            COEF[2][2][4][i4] = 5.0/4.0*NG*(9.0/FN5-6.0/FN3+1.0/FN1);  //  (1/7) * ( 7 * M0 + 10 * M2 + 18 * M4 )
                            COEFF = NG*15.0/2.0*(1.0/FN3-1.0/FN5);                     //  (1/7) * ( 7 * M0 + 5 * M2 - 12 * M4 )
                            COEF[2][2][3][i4] = COEFF;
                            COEF[2][2][5][i4] = COEFF;
                            COEFF = 15.0/8.0*NG*(1.0/FN1-2.0/FN3+1.0/FN5);             //  (1/7) * ( 7 * M0 - 10 * M2 + 3 * M4 )
                            COEF[2][2][2][i4] = COEFF;
                            COEF[2][2][6][i4] = COEFF;
                        }

                    }

                    if (_lmax_dft_ecp > 2) {

                        double COEFF = NU*.5*SQ7 * (5./FN4 - 3./FN2); //  SQ(7) * M3
                        COEF[0][3][4][i4] = COEFF;
                        COEF[3][0][4][i4] = COEFF;

                        if (_lmin_dft_ecp > 0) {
                            COEFF = NG*.5*SQ3*SQ7 * (5./FN5 - 3./FN3);              //  ( SQ(21)/7 ) * ( 3 * M2 + 4 * M4 )
                            COEF[1][3][4][i4] = COEFF;
                            COEF[3][1][4][i4] = COEFF;
                            COEFF = NG*.375*SQ2*SQ7 * (-5./FN5 + 6./FN3 - 1./FN1);  //  ( (3*SQ(14))/7 ) * ( M2 - M4 )
                            COEF[1][3][3][i4] = COEFF;
                            COEF[1][3][5][i4] = COEFF;
                            COEF[3][1][3][i4] = COEFF;
                            COEF[3][1][5][i4] = COEFF;
                        }

                        if (_lmin_dft_ecp > 1) {
                            COEFF = NU*.25*SQ5*SQ7 * (15./FN6 - 14./FN4 + 3./FN2);      //  ( SQ(35)/105 ) * ( 27 * M1 + 28 * M3 + 50 * M5 )
                            COEF[2][3][4][i4] = COEFF;
                            COEF[3][2][4][i4] = COEFF;
                            COEFF = NU*.375*SQ2*SQ5*SQ7 * (-5./FN6 + 6./FN4 - 1./FN2);  //  ( SQ(70)/105 ) * ( 18 * M1 + 7 * M3 - 25 * M5 )
                            COEF[2][3][3][i4] = COEFF;
                            COEF[2][3][5][i4] = COEFF;
                            COEF[3][2][3][i4] = COEFF;
                            COEF[3][2][5][i4] = COEFF;
                            COEFF = NU*1.875*SQ7 * (1./FN6 - 2./FN4 + 1./FN2);          //  ( SQ(7)/21 ) * ( 9 * M1 - 14 * M3 + 5 * M5 )
                            COEF[2][3][2][i4] = COEFF;
                            COEF[2][3][6][i4] = COEFF;
                            COEF[3][2][2][i4] = COEFF;
                            COEF[3][2][6][i4] = COEFF;
                        }

                        if (_lmin_dft_ecp > 2) {
                            COEF[3][3][4][i4] = NG*1.75 * (50./FN7 - 30./FN5 + 9./FN3);   //  (1/33) * ( 33 * M0 + 44 * M2 + 54 * M4 + 100 * M6 )
                            COEFF = NG*1.3125 * (-25./FN7 + 35./FN5 - 11./FN3 + 1./FN1);  //  (1/11) * ( 11 * M0 + 11 * M2 + 3 * M4 - 25 * M6 )
                            COEF[3][3][3][i4] = COEFF;
                            COEF[3][3][5][i4] = COEFF;
                            COEFF = NG*105*.125 * (1./FN7 - 2./FN5 + 1./FN3);             //  (1/11) * ( 11 * M0 - 21 * M4 + 10 * M6 )
                            COEF[3][3][2][i4] = COEFF;
                            COEF[3][3][6][i4] = COEFF;
                            COEFF = NG*35*.0625 * (-1./FN7 + 3./FN5 - 3./FN3 + 1./FN1);   //  (1/33) * ( 33 * M0 - 55 * M2 + 27 * M4 - 5 * M6 )
                            COEF[3][3][1][i4] = COEFF;
                            COEF[3][3][7][i4] = COEFF;
                        }

                    }


                    if (_lmax_dft_ecp > 3) {

                        double COEFF = NG*.375 * (35./FN5 - 30./FN3 + 3./FN1);  //  3 * M4
                        COEF[0][4][4][i4] = COEFF;
                        COEF[4][0][4][i4] = COEFF;

                        if (_lmin_dft_ecp > 0) {
                            COEFF = NU*.375*SQ3 * (35./FN6 - 30./FN4 + 3./FN2);          //  ( SQ(3)/3 ) * ( 4 * M3 + 5 * M5 )
                            COEF[1][4][4][i4] = COEFF;
                            COEF[4][1][4][i4] = COEFF;
                            COEFF = NU*.375*SQ2*SQ3*SQ5 * (-7./FN6 + 10./FN4 - 3./FN2);  //  ( SQ(30)/3 ) * ( M3 - M5 )
                            COEF[1][4][3][i4] = COEFF;
                            COEF[1][4][5][i4] = COEFF;
                            COEF[4][1][3][i4] = COEFF;
                            COEF[4][1][5][i4] = COEFF;
                        }

                        if (_lmin_dft_ecp > 1) {
                            COEFF = NG*.1875*SQ5 * (105./FN7 - 125./FN5 + 39./FN3 - 3./FN1);  //  ( (3*SQ(5))/77 ) * ( 22 * M2 + 20 * M4 + 35 * M6 )
                            COEF[2][4][4][i4] = COEFF;
                            COEF[4][2][4][i4] = COEFF;
                            COEFF = NG*1.875*SQ2*SQ3 * (-7./FN7 + 10./FN5 - 3./FN3);          //  ( (5*SQ(6))/77 ) * ( 11 * M2 + 3 * M4 - 14 * M6 )
                            COEF[2][4][3][i4] = COEFF;
                            COEF[2][4][5][i4] = COEFF;
                            COEF[4][2][3][i4] = COEFF;
                            COEF[4][2][5][i4] = COEFF;
                            COEFF = NG*.9375*SQ3 * (7./FN7 - 15./FN5 + 9./FN3 - 1./FN1);      //  ( (5*SQ(3))/77 ) * ( 11 * M2 - 18 * M4 + 7 * M6 )
                            COEF[2][4][2][i4] = COEFF;
                            COEF[2][4][6][i4] = COEFF;
                            COEF[4][2][2][i4] = COEFF;
                            COEF[4][2][6][i4] = COEFF;
                        }

                        if (_lmin_dft_ecp > 2) {
                            COEFF = NU*.1875*SQ7 * (175./FN8 - 255./FN6 + 105./FN4 - 9./FN2);        //  ( SQ(7)/1001 ) * ( 572 * M1 + 546 * M3 + 660 * M5 + 1225 * M7)
                            COEF[3][4][4][i4] = COEFF;
                            COEF[4][3][4][i4] = COEFF;
                            COEFF = NU*1.875*SQ3*SQ5*SQ7 * (-35./FN8 + 57./FN6 - 25./FN4 + 3./FN2);  //  ( SQ(105))/1001 ) * ( 143 * M1 + 91 * M3 + 11 * M5 - 245 * M7)
                            COEF[3][4][3][i4] = COEFF;
                            COEF[3][4][5][i4] = COEFF;
                            COEF[4][3][3][i4] = COEFF;
                            COEF[4][3][5][i4] = COEFF;
                            COEFF = NU*.9375*SQ3*SQ7 * (7./FN8 - 15./FN6 + 9./FN4 - 1./FN2);         //  ( SQ(21))/1001 ) * ( 286 * M1 - 91 * M3 - 440 * M5 + 245 * M7)
                            COEF[3][4][2][i4] = COEFF;
                            COEF[3][4][6][i4] = COEFF;
                            COEF[4][3][2][i4] = COEFF;
                            COEF[4][3][6][i4] = COEFF;
                            COEFF = NU*6.5625 * (-1./FN8 + 3./FN6 - 3./FN4 + 1./FN2);                //  ( 1/143 ) * ( 143 * M1 - 273 * M3 + 165 * M5 - 35 * M7)
                            COEF[3][4][1][i4] = COEFF;
                            COEF[3][4][7][i4] = COEFF;
                            COEF[4][3][1][i4] = COEFF;
                            COEF[4][3][7][i4] = COEFF;
                        }

                        if (_lmin_dft_ecp > 3) {

                             cout << "Sorry, not yet supported: Combination of G functions in DFT basis and ECPs with l = 4." << endl;
                             exit(1);

                        }

                    }

                } // i4 loop (== CKO )

            } // end if (INULL != 0)



            type_3D BLMA;
            type_3D CA;
            getBLMCOF(_lmax_ecp, _lmax_row, AVS, BLMA, CA);

            type_3D BLMB;
            type_3D CB;
            getBLMCOF(_lmax_ecp, _lmax_col, BVS, BLMB, CB);

            typedef boost::multi_array_types::extent_range range;
            typedef type_3D::index index;
            type_3D::extent_gen extents3D;


            switch (INULL) {

                case 0:  //  AVSSQ <= 0.1 && BVSSQ <= 0.1
                {

                    for (int i = 0; i < _nsph_row; i++) {
                        for (int j = 0; j < _nsph_col; j++) {
                            for (index L = 0; L <= _lmin; L++) {
                                double XI_L = XI(L, L + L);
                                for (index M = 4 - L; M <= 4 + L; M++) {
                                    matrix(i,j) += BLMA[i][L][M] * BLMB[j][L][M] * XI_L;
                                }
                            }
                        }
                    }
                    break;
                }


                case 1:  //  AVSSQ <= 0.1
                {

                    type_3D SUMCI3;
                    SUMCI3.resize(extents3D[range(0, 5)][range(0, 5)][range(0, 9)]);
                    for (index L = 0; L <= _lmax_ecp; L++) {
                        for (index L2 = 0; L2 <= _lmax_col; L2++) {
                            int range_M2 = std::min(L2, L);
                            for (index M2 = 4 - range_M2; M2 <= 4 + range_M2; M2++) {

                                double VAR2 = 0.0;
                                double fak = 2.0 * beta * BVSSQ;
                                double pow = 1;
                                double factorialNN = 1;

                                for (int NN = 0; NN <= NMAX2; NN++) {

                                    if (NN != 0) {
                                        pow = pow * fak;
                                        factorialNN = factorialNN * NN;
                                    }

                                    double XDUM = COEF[L][L2][M2][NN] * pow / factorialNN;
                                    VAR2 += XDUM * XI(L, NN + L + L2);

                                }

                                SUMCI3[L][L2][M2] = VAR2;

                            } // end M2
                        } // end L2
                    } // end L

                    for (int i = 0; i < _nsph_row; i++) {
                        for (int j = 0; j < _nsph_col; j++) {

                            for (index L = 0; L <= std::min(_lmax_row, _lmax_ecp); L++) {
                                for (index L2 = 0; L2 <= _lmax_col; L2++) {
                                    int range_M2 = std::min(L2, L);
                                    for (index M2 = 4 - range_M2; M2 <= 4 + range_M2; M2++) {

                                        for (index M1 = 4 - L; M1 <= 4 + L; M1++) {
                                            matrix(i,j) += BLMA[i][L][M1] * BLMB[j][L2][M2] * SUMCI3[L][L2][M2] * CB[L][M1][M2];
                                        }
                                    }
                                }
                            }
                        }
                    }

                    break;
                }


                case 2:  //  BVSSQ <= 0.1
                {

                    type_3D SUMCI3;
                    SUMCI3.resize(extents3D[range(0, 5)][range(0, 5)][range(0, 9)]);
                    for (index L = 0; L <= _lmax_ecp; L++) {
                        for (index L1 = 0; L1 <= _lmax_row; L1++) {
                            int range_M1 = std::min(L1, L);
                            for (index M1 = 4 - range_M1; M1 <= 4 + range_M1; M1++) {

                                double VAR1 = 0.0;
                                double fak = 2.0 * alpha * AVSSQ;
                                double pow = 1;
                                double factorialN = 1;

                                for (int N = 0; N <= NMAX1; N++) {

                                    if (N != 0) {
                                        pow = pow * fak;
                                        factorialN = factorialN * N;
                                    }

                                    double XDUM = COEF[L][L1][M1][N] * pow / factorialN;
                                    VAR1 += XDUM * XI(L, N + L1 + L);

                                }

                                SUMCI3[L][L1][M1] = VAR1;

                            } // end M1
                        } // end L1
                    } // end L


                    for (int i = 0; i < _nsph_row; i++) {
                        for (int j = 0; j < _nsph_col; j++) {

                            for (index L = 0; L <= std::min(_lmax_col, _lmax_ecp); L++) {
                                for (index L1 = 0; L1 <= _lmax_row; L1++) {
                                    int range_M1 = std::min(L1, L);
                                    for (index M1 = 4 - range_M1; M1 <= 4 + range_M1; M1++) {

                                        for (index M2 = 4 - L; M2 <= 4 + L; M2++) {
                                            matrix(i,j) += BLMA[i][L1][M1] * BLMB[j][L][M2] * SUMCI3[L][L1][M1] * CA[L][M2][M1];
                                        }
                                    }
                                }
                            }
                        }
                    }

                    break;
                }


                case 3:
                {

                    type_3D CC;
                    CC.resize(extents3D[range(0, 5)][range(0, 9)][range(0, 9)]);
                    for (index L = 0; L <= _lmax_ecp; L++) {
                        int range_M1 = std::min(_lmax_row, int(L));
                        int range_M2 = std::min(_lmax_col, int(L));
                        for ( index M1 = 4 - range_M1; M1 <= 4 + range_M1; M1++) {
                            for (index M2 = 4 - range_M2; M2 <= 4 + range_M2; M2++) {
                                CC[L][M1][M2] = 0.0;
                                for (index M = 4 - L; M <= 4 + L; M++) {
                                    CC[L][M1][M2] += CA[L][M][M1] * CB[L][M][M2];
                                }
                            }
                        }
                    }

                    typedef boost::multi_array<double, 5> type_5D;
                    type_5D::extent_gen extents5D;
                    type_5D SUMCI;
                    SUMCI.resize(extents5D[range(0, 5)][range(0, 5)][range(0, 5)][range(0, 9)][range(0, 9)]);
                    for (index L = 0; L <= _lmax_ecp; L++) {
                        for (index L1 = 0; L1 <= _lmax_row; L1++) {
                            int range_M1 = std::min(L1, L);
                            for (index L2 = 0; L2 <= _lmax_col; L2++) {
                                int range_M2 = std::min(L2, L);
                                for (index M1 = 4 - range_M1; M1 <= 4 + range_M1; M1++) {
                                    for (index M2 = 4 - range_M2; M2 <= 4 + range_M2; M2++) {

                                        SUMCI[L][L1][L2][M1][M2]  = 0.0;

                                        double fak1 = 2.0 * alpha * AVSSQ;
                                        double pow1 = 1;
                                        double factorialN = 1;

                                        for (int N = 0; N <= NMAX1; N++ ) {

                                            if (N != 0) {
                                            pow1 = pow1 * fak1;
                                            factorialN = factorialN*N;
                                            }

                                            double VAR1 = COEF[L][L1][M1][N] * pow1 / factorialN;
                                            double VAR2 = 0.0;
                                            double fak2 = 2.0 * beta * BVSSQ;
                                            double pow2 = 1;
                                            double factorialNN = 1;

                                            for (int NN = 0; NN <= NMAX2; NN++) {

                                                if (NN != 0) {
                                                    pow2 = pow2 * fak2;
                                                    factorialNN = factorialNN * NN;
                                                }
                                                double XDUM = COEF[L][L2][M2][NN] * pow2 / factorialNN;
                                                VAR2  += XDUM * XI(L,N+NN+L1+L2);

                                            }

                                            SUMCI[L][L1][L2][M1][M2]  += VAR1 * VAR2;

                                        }

                                    } // end M2
                                } // end M1
                            } // end L2
                        } // end L1 
                    } // end L

                    for (int i = 0; i < _nsph_row; i++) {
                        for (int j = 0; j < _nsph_col; j++) {

                            for (index L = 0; L <= _lmax_ecp; L++) {
                                for (index L1 = 0; L1 <= _lmax_row; L1++) {
                                    int range_M1 = std::min(L1, L);
                                    for (index L2 = 0; L2 <= _lmax_col; L2++) {
                                        int range_M2 = std::min(L2, L);
                                        for (index M1 = 4 - range_M1; M1 <= 4 + range_M1; M1++) {
                                            for (index M2 = 4 - range_M2; M2 <= 4 + range_M2; M2++) {

                                                matrix(i,j) += BLMA[i][L1][M1] * BLMB[j][L2][M2] * SUMCI[L][L1][L2][M1][M2] * CC[L][M1][M2];

                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    break;
                }


                default:
                    cout << "Wrong ECP summation mode";
                    exit(1);
            } // switch


            // GET TRAFO HERE ALREADY
            ub::vector<double> NormA = CalcNorms(alpha, _nsph_row);
            ub::vector<double> NormB = CalcNorms(beta, _nsph_col);

            for (int i = 0; i < _nsph_row; i++) {
                for (int j = 0; j < _nsph_col; j++) {

                    matrix(i,j) = matrix(i,j) * GAUSS * NormA[i] * NormB[j];

                }
            }


            return matrix;
        }


        ub::vector<double> AOECP::CalcNorms(double decay, int size) {
            ub::vector<double> Norms = ub::vector<double>(size);
            const double PI = boost::math::constants::pi<double>();
            double SQ2, SQ3, SQ5;

            double Norm_S = pow(2.0 * decay / PI, 0.75);
            double Norm_P=0.0;
            double Norm_D=0.0;
            Norms[0] = Norm_S;  //  Y 00

            if (size > 1) {
                Norm_P = 2.0 * sqrt(decay) * Norm_S;
                Norms[1] = Norm_P;  //  Y 10
                Norms[2] = Norm_P;  //  Y 1-1
                Norms[3] = Norm_P;  //  Y 11
            }

            if (size > 4) {
                SQ3 = sqrt(3.);
                double Norm_D = 4.00 * decay * Norm_S;
                Norms[4] = .5 * Norm_D / SQ3;  //  Y 20
                Norms[5] = Norm_D;             //  Y 2-1
                Norms[6] = Norm_D;             //  Y 21
                Norms[7] = Norm_D;             //  Y 2-2
                Norms[8] = .5 * Norm_D;        //  Y 22
            }

            if (size > 9) {
                SQ2 = sqrt(2.);
                SQ5 = sqrt(5.);
                double Norm_F = 4.00 * decay * Norm_P;
                double Norm_F_1 = .5 * Norm_F / (SQ2 * SQ5);
                double Norm_F_3 = .5 * Norm_F / (SQ2 * SQ3);
                Norms[9] =  .5 * Norm_F / (SQ3 * SQ5);  //  Y 30
                Norms[10] = Norm_F_1;                   //  Y 3-1
                Norms[11] = Norm_F_1;                   //  Y 31
                Norms[12] = Norm_F;                     //  Y 3-2
                Norms[13] = .5 * Norm_F;                //  Y 32
                Norms[14] = Norm_F_3;                   //  Y 3-3
                Norms[15] = Norm_F_3;                   //  Y 33
            }

            if (size > 16) {
                double SQ7 = sqrt(7.);
                double Norm_G = 4.00 * decay * Norm_D;
                double Norm_G_1 = .5 * Norm_G / (SQ2 * SQ3 * SQ7);
                double Norm_G_m2 = .5 * Norm_G / (SQ3 * SQ7);
                double Norm_G_3 = .5 * Norm_G / (SQ2 * SQ3);
                double Norm_G_m4 = .5 * Norm_G / SQ3;
                Norms[16] =  .125 * Norm_G / (SQ3 * SQ5 * SQ7);  //  Y 40
                Norms[17] = Norm_G_1;                            //  Y 4-1
                Norms[18] = Norm_G_1;                            //  Y 41
                Norms[19] = Norm_G_m2;                           //  Y 4-2
                Norms[20] = .5 * Norm_G_m2;                      //  Y 42
                Norms[21] = Norm_G_3;                            //  Y 4-3
                Norms[22] = Norm_G_3;                            //  Y 43
                Norms[23] = Norm_G_m4;                           //  Y 4-4
                Norms[24] = .25 * Norm_G_m4;                     //  Y 44
            }
            return Norms;
        }

        void AOECP::getBLMCOF(int _lmax_ecp, int _lmax_dft, const vec& pos, type_3D& BLC, type_3D& C) {

            typedef boost::multi_array_types::extent_range range;
            typedef type_3D::index index;
            type_3D::extent_gen extents;

            int _nsph = (_lmax_dft + 1) * (_lmax_dft + 1);
            int _lmax_dft_ecp = std::max(_lmax_dft, _lmax_ecp);
            int _lmin_dft_ecp = std::min(_lmax_dft, _lmax_ecp);

            BLC.resize(extents[range(0, _nsph)][range(0, 5)][range(0, 9)]);
///            C.resize(extents[range(0, _lmax_dft_ecp+1)][range(0, 9)][range(0, 9)]);
            C.resize(extents[range(0, 5)][range(0, 9)][range(0, 9)]);

            type_3D BLM;
///            BLM.resize(extents[range(0, _nsph)][range(0, _lmax_dft+1)][range(0, 9)]);
            BLM.resize(extents[range(0, _nsph)][range(0, 5)][range(0, 9)]);

            const double PI = boost::math::constants::pi<double>();
            double SQPI = sqrt(PI);
            double SQ2 = sqrt(2.);
            double SQ2PI = SQ2 * SQPI;
            double SQ3 = sqrt(3.);
            double SQ5 = sqrt(5.);
            double SQ7 = sqrt(7.);
            double XS = 2. * SQPI;

            double XP, XD, XD_0, XD_p2;
            double XF_0, XF_1, XF_m2, XF_p2, XF_3;
      
            for (index I = 0; I < _nsph; I++) {
                for (index L = 0; L <= _lmax_dft; L++) {
                    for (index M = 4 - L; M <= 4 + L; M++) {
                        BLM[I][L][M] = 0.0;
                    }
                }
            }

            double BVS_X = pos.getX();
            double BVS_Y = pos.getY(); 
            double BVS_Z = pos.getZ();
            double BVS_XX = BVS_X * BVS_X;
            double BVS_YY = BVS_Y * BVS_Y;
            double BVS_ZZ = BVS_Z * BVS_Z;
            double BVS_RR = BVS_XX + BVS_YY + BVS_ZZ;
            double BVS_XY, BVS_XZ, BVS_YZ;

            BLM[0][0][4] = XS;  //  Y 00

            if (_lmax_dft > 0) {

                XP = XS / SQ3;

                BLM[1][0][4] = -BVS_Z * XS;  //  Y 10
                BLM[1][1][4] = XP;

                BLM[2][0][4] = -BVS_Y * XS;  //  Y 1 -1
                BLM[2][1][3] = XP;

                BLM[3][0][4] = -BVS_X * XS;  //  Y 11
                BLM[3][1][5] = XP;

            }

            if (_lmax_dft > 1) {

                BVS_XY = BVS_X * BVS_Y;
                BVS_XZ = BVS_X * BVS_Z;
                BVS_YZ = BVS_Y * BVS_Z;

                XD = XP / SQ5;
                XD_0 = 4.0 * SQPI / SQ5;
                XD_p2 = 2. * XD;

                BLM[4][0][4] = (3.0 * BVS_ZZ - BVS_RR) * XS;  //  Y 20
                BLM[4][1][5] = 2.0 * BVS_X * XP;
                BLM[4][1][3] = 2.0 * BVS_Y * XP;
                BLM[4][1][4] = -4.0 * BVS_Z * XP;
                BLM[4][2][4] = XD_0;

                BLM[5][0][4] = BVS_YZ * XS;  //  Y 2 -1
                BLM[5][1][3] = -BVS_Z * XP;
                BLM[5][1][4] = -BVS_Y * XP;
                BLM[5][2][3] = XD;

                BLM[6][0][4] = BVS_XZ * XS;  //  Y 21
                BLM[6][1][5] = -BVS_Z * XP;
                BLM[6][1][4] = -BVS_X * XP;
                BLM[6][2][5] = XD;

                BLM[7][0][4] = BVS_XY * XS;  //  Y 2 -2
                BLM[7][1][5] = -BVS_Y * XP;
                BLM[7][1][3] = -BVS_X * XP;
                BLM[7][2][2] = XD;

                BLM[8][0][4] = (BVS_XX - BVS_YY) * XS;  //  Y 22
                BLM[8][1][5] = -2.0 * BVS_X * XP;
                BLM[8][1][3] = 2.0 * BVS_Y * XP;
                BLM[8][2][6] = XD_p2;
            }

            if (_lmax_dft > 2) {

                XF_0 = 4. * SQPI / SQ7;
                XF_1 = 4. * SQ2PI / (SQ3 * SQ7);
                XF_m2 = 2. * SQPI / (SQ3 * SQ5 * SQ7);
                XF_p2 = 2. * XF_m2;
                XF_3 = 4. * SQ2PI / (SQ5 * SQ7);

                BLM[9][0][4] = (3. * BVS_RR - 5. * BVS_ZZ) * BVS_Z * XS;  //  Y 30
                BLM[9][1][3] = -6.* BVS_YZ * XP;
                BLM[9][1][4] =  3. * (3. * BVS_ZZ - BVS_RR) * XP;
                BLM[9][1][5] = -6. * BVS_XZ * XP;
                BLM[9][2][3] = 6. * BVS_Y * XD;
                BLM[9][2][4] = -3. * BVS_Z * XD_0;
                BLM[9][2][5] = 6. * BVS_X * XD;
                BLM[9][3][4] = XF_0;

                BLM[10][0][4] = (BVS_RR - 5. * BVS_ZZ) * BVS_Y * XS;  //  Y 3 -1
                BLM[10][1][3] = (4. * BVS_ZZ - BVS_XX - 3. * BVS_YY)  * XP;
                BLM[10][1][4] = 8. * BVS_YZ * XP;
                BLM[10][1][5] = -2. * BVS_XY * XP;
                BLM[10][2][2] = 2. * BVS_X * XD;
                BLM[10][2][3] = -8. * BVS_Z * XD;
                BLM[10][2][4] = -2. * BVS_Y * XD_0;
                BLM[10][2][6] = -BVS_Y * XD_p2;
                BLM[10][3][3] = XF_1;

                BLM[11][0][4] = (BVS_RR - 5. * BVS_ZZ) * BVS_X * XS;  //  Y 31
                BLM[11][1][3] = -2. * BVS_XY * XP;
                BLM[11][1][4] = 8. * BVS_XZ * XP;
                BLM[11][1][5] = (4. * BVS_ZZ - 3. * BVS_XX - BVS_YY)  * XP;
                BLM[11][2][2] = 2. * BVS_Y * XD;
                BLM[11][2][4] = -2. * BVS_X * XD_0;
                BLM[11][2][5] = -8. * BVS_Z * XD;
                BLM[11][2][6] = BVS_X * XD_p2;
                BLM[11][3][5] = XF_1;

                BLM[12][0][4] = -BVS_XY * BVS_Z * XS;  //  Y 3 -2
                BLM[12][1][3] = BVS_XZ * XP;
                BLM[12][1][4] = BVS_XY * XP;
                BLM[12][1][5] = BVS_YZ * XP;
                BLM[12][2][2] = -BVS_Z * XD;
                BLM[12][2][3] = -BVS_X * XD;
                BLM[12][2][5] = -BVS_Y * XD;
                BLM[12][3][2] = XF_m2;

                BLM[13][0][4] = (BVS_YY - BVS_XX) * BVS_Z * XS;  //  Y 32
                BLM[13][1][3] = -2. * BVS_YZ * XP;
                BLM[13][1][4] = (BVS_XX - BVS_YY)* XP;
                BLM[13][1][5] = 2. * BVS_XZ * XP;
                BLM[13][2][3] = 2. * BVS_Y * XD;
                BLM[13][2][5] = -2. * BVS_X * XD;
                BLM[13][2][6] = -BVS_Z * XD_p2;
                BLM[13][3][6] = XF_p2;

                BLM[14][0][4] = (BVS_YY - 3. * BVS_XX) * BVS_Y * XS;  //  Y 3 -3
                BLM[14][1][3] = 3. * (BVS_XX - BVS_YY) * XP;
                BLM[14][1][5] = 6. * BVS_XY * XP;
                BLM[14][2][2] = -6. * BVS_X * XD;
                BLM[14][2][6] = -3. * BVS_Y * XD_p2;
                BLM[14][3][1] = XF_3;

                BLM[15][0][4] = (3. * BVS_YY - BVS_XX) * BVS_X * XS;  //  Y 33
                BLM[15][1][3] = -6. * BVS_XY * XP;
                BLM[15][1][5] = 3. * (BVS_XX - BVS_YY) * XP;
                BLM[15][2][2] = 6. * BVS_Y * XD;
                BLM[15][2][6] = -3. * BVS_X * XD_p2;
                BLM[15][3][7] = XF_3;
            }

            if (_lmax_dft > 3) {

                double XG_0 = 16. * SQPI / 3.;
                double XG_1 = 4. * SQ2PI / (3. * SQ5);
                double XG_m2 = 4. * SQPI / (3. * SQ5);
                double XG_p2 = 2. * XG_m2;
                double XG_3 = 4. * SQ2PI / (3. * SQ5 * SQ7);
                double XG_m4 = 4. * SQPI / (3. * SQ5 * SQ7);
                double XG_p4 = 4. * XG_m4;

                BLM[16][0][4] = (35. * BVS_ZZ - 30. * BVS_ZZ * BVS_RR + 3. * BVS_RR * BVS_RR) * XS;  //  Y 40
                BLM[16][1][3] = 12. * (5. * BVS_ZZ - BVS_RR) * BVS_Y * XP;
                BLM[16][1][4] = 16. * (3. * BVS_RR - 5. * BVS_ZZ) * BVS_Z * XP;
                BLM[16][1][5] = 12. * (5. * BVS_ZZ - BVS_RR) * BVS_X * XP;
                BLM[16][2][2] = -24. * BVS_XY * XD;
                BLM[16][2][3] = -96. * BVS_YZ * XD;
                BLM[16][2][4] = 12. * (3. * BVS_ZZ - BVS_RR) * XD_0;
                BLM[16][2][5] = -96. * BVS_XZ * XD;
                BLM[16][2][6] = 6. * (BVS_XX - BVS_YY) * XD_p2;
                BLM[16][3][3] = 12. * BVS_Y * XF_1;
                BLM[16][3][4] = -16. * BVS_Z * XF_0;
                BLM[16][3][5] = 12. * BVS_X * XF_1;
                BLM[16][4][4] = XG_0;

                BLM[17][0][4] = (7. * BVS_ZZ - 3. * BVS_RR) * BVS_YZ * XS;  //  Y  4 -1
                BLM[17][1][3] = (3. * BVS_XX + 9. * BVS_YY - 4. * BVS_ZZ) * BVS_Z * XP;
                BLM[17][1][4] = 3. * (BVS_RR - 5. * BVS_ZZ) * BVS_Y * XP;
                BLM[17][1][5] = 6. * BVS_XY * BVS_Z * XP;
                BLM[17][2][2] = -6. * BVS_XZ * XD;
                BLM[17][2][3] = 3. * (4. * BVS_ZZ - BVS_XX - 3. * BVS_YY) * XD;
                BLM[17][2][4] = 6. * BVS_YZ * XD_0;
                BLM[17][2][5] = -6. * BVS_XY * XD;
                BLM[17][2][6] = 3. * BVS_YZ * XD_p2;
                BLM[17][3][2] = 6. * BVS_X * XF_m2;
                BLM[17][3][3] = -3. * BVS_Z * XF_1;
                BLM[17][3][4] = -2. * BVS_Y * XF_0;
                BLM[17][3][6] = -3. * BVS_Y * XF_p2;
                BLM[17][4][3] = XG_1;

                BLM[18][0][4] = (7. * BVS_ZZ - 3. * BVS_RR) * BVS_XZ * XS;  //  Y 41
                BLM[18][1][3] = 6. * BVS_XY * BVS_Z * XP;
                BLM[18][1][4] = 3. * (BVS_RR - 5. * BVS_ZZ) * BVS_Y * XP;
                BLM[18][1][5] = (9. * BVS_XX + 3. * BVS_YY - 4. * BVS_ZZ) * BVS_Z * XP;
                BLM[18][2][2] = -6. * BVS_YZ * XD;
                BLM[18][2][3] = -6. * BVS_XY * XD;
                BLM[18][2][4] = 6. * BVS_XZ * XD_0;
                BLM[18][2][5] = 3. * (4. * BVS_ZZ - 3. * BVS_XX - BVS_YY) * XD;
                BLM[18][2][6] = -3. * BVS_XZ * XD_p2;
                BLM[18][3][2] = 6. * BVS_Y * XF_m2;
                BLM[18][3][4] = -2. * BVS_X * XF_0;
                BLM[18][3][5] = -3. * BVS_Z * XF_1;
                BLM[18][3][6] = 3. * BVS_X * XF_p2;
                BLM[18][4][5] = XG_1;

                BLM[19][0][4] = (7. * BVS_ZZ - BVS_RR) * BVS_XY * XS;  //  Y 4 -2
                BLM[19][1][3] = (BVS_XX + 3. * BVS_YY - 6. * BVS_ZZ) * BVS_X * XP;
                BLM[19][1][4] = -12. * BVS_XY *BVS_Z * XP;
                BLM[19][1][5] = (3. * BVS_XX + BVS_YY - 6. * BVS_ZZ) * BVS_Y * XP;
                BLM[19][2][2] = 3.* (3. * BVS_ZZ - BVS_RR) * XD;
                BLM[19][2][3] = 12. * BVS_XZ * XD;
                BLM[19][2][4] = 3. * BVS_XY * XD_0;
                BLM[19][2][5] = 12. * BVS_YZ * XD;
                BLM[19][3][1] = .5 * BVS_X * XF_3;
                BLM[19][3][2] = -12. * BVS_Z * XF_m2;
                BLM[19][3][3] = -1.5 * BVS_X * XF_1;
                BLM[19][3][5] = -1.5 * BVS_Y * XF_1;
                BLM[19][3][7] = -.5 * BVS_Y * XF_3;
                BLM[19][4][2] = XG_m2;

                BLM[20][0][4] = (7. * BVS_ZZ - BVS_RR) * (BVS_XX - BVS_YY) * XS;  //  Y 42
                BLM[20][1][3] = 4. * (3. * BVS_ZZ - BVS_YY) * BVS_Y * XP;
                BLM[20][1][4] = 12. * (BVS_YY - BVS_XX) * BVS_Z * XP;
                BLM[20][1][5] = 4. *(BVS_XX - 3. * BVS_ZZ) * BVS_X * XP;
                BLM[20][2][3] = -24. * BVS_YZ * XD;
                BLM[20][2][4] = 3. * (BVS_XX - BVS_YY) * XD_0;
                BLM[20][2][5] = 24. * BVS_XZ * XD;
                BLM[20][2][6] = 3. * (3. * BVS_ZZ - BVS_RR) * XD_p2;
                BLM[20][3][1] = BVS_Y * XF_3;
                BLM[20][3][3] = 3. * BVS_Y * XF_1;
                BLM[20][3][5] = -3. * BVS_X * XF_1;
                BLM[20][3][6] = -12. * BVS_Z * XF_m2;
                BLM[20][3][7] = BVS_X * XF_3;
                BLM[20][4][6] = XG_p2;

                BLM[21][0][4] = (3. * BVS_XX - BVS_YY) * BVS_YZ * XS;  //  Y 4 -3
                BLM[21][1][3] = 3. * (BVS_YY - BVS_XX) *BVS_Z * XP;
                BLM[21][1][4] = (BVS_YY - 3. * BVS_XX) *BVS_Y * XP;
                BLM[21][1][5] = -6. * BVS_XY * BVS_Z * XP;
                BLM[21][2][2] = 6. * BVS_XZ * XD;
                BLM[21][2][3] = 3. * (BVS_XX - BVS_YY) * XD;
                BLM[21][2][5] = 6. * BVS_XY * XD;
                BLM[21][2][6] = 3. * BVS_YZ * XD_p2;
                BLM[21][3][1] = -BVS_Z * XF_3;
                BLM[21][3][2] = -6. * BVS_X * XF_m2;
                BLM[21][3][6] = -3. * BVS_Y * XF_p2;
                BLM[21][4][1] = XG_3;

                BLM[22][0][4] = (BVS_XX - 3. * BVS_YY) * BVS_XZ * XS;  //  Y 43
                BLM[22][1][3] = 6. * BVS_XY * BVS_Z * XP;
                BLM[22][1][4] = (3. * BVS_YY - BVS_XX) * BVS_X * XP;
                BLM[22][1][5] = 3. * (BVS_YY - BVS_XX)* BVS_Z * XP;
                BLM[22][2][2] = -6. * BVS_YZ * XD;
                BLM[22][2][3] = -6. * BVS_XY * XD;
                BLM[22][2][5] = 3. * (BVS_XX - BVS_YY) * XD;
                BLM[22][2][6] = 3. * BVS_XZ * XD_p2;
                BLM[22][3][2] = 6. * BVS_Y * XF_m2;
                BLM[22][3][6] = -3. * BVS_X * XF_p2;
                BLM[22][3][7] = -BVS_Z * XF_3;
                BLM[22][4][7] = XG_3;

                BLM[23][0][4] = (BVS_XX - BVS_YY) * BVS_XY * XS;  //  Y 4 -4
                BLM[23][1][3] = (3. * BVS_YY - BVS_XX) * BVS_X * XP;
                BLM[23][1][5] = (BVS_YY - 3.* BVS_XX) * BVS_Y * XP;
                BLM[23][2][2] = 3. * (BVS_XX - BVS_YY) * XD;
                BLM[23][2][6] = 3. * BVS_XY * XD_p2;
                BLM[23][3][1] = -BVS_X * XF_3;
                BLM[23][3][7] = -BVS_Y * XF_3;
                BLM[23][4][0] = XG_m4;

                BLM[24][0][4] = (BVS_XX * BVS_XX + BVS_YY * BVS_YY - 6. * BVS_XX * BVS_YY) * XS;  //  Y 4 4
                BLM[24][1][3] = 4. * (3. * BVS_XX - BVS_YY) * BVS_Y * XP;
                BLM[24][1][5] = 4. * (3. * BVS_YY - BVS_XX) * BVS_X * XP;
                BLM[24][2][2] = -24. * BVS_XY * XD;
                BLM[24][2][6] = 6. * (BVS_XX - BVS_YY) * XD_p2;
                BLM[24][3][1] = 4. * BVS_Y * XF_3;
                BLM[24][3][7] = -4. * BVS_X * XF_3;
                BLM[24][4][8] = XG_p4;

            }



            for (index L = 0; L <= _lmax_dft_ecp; L++) {
                for (index M = 4 - L; M <= 4 + L; M++) {
                    for (index MM = 4 - L; MM <= 4 + L; MM++) {
                        C[L][M][MM] = 0.0;
                    }
                }
            }
            double SXY = sqrt(BVS_XX + BVS_YY);  // SXY = r * sin(theta)
            double SXYZ = sqrt(BVS_RR);  // SXYZ = r

            double CP = 1.0;
            double SP = 0.0;
            if (SXY > 1.e-4) {
                CP = BVS_X / SXY;  // CP = cos(phi)
                SP = BVS_Y  /SXY;  // SP = sin(phi)
            }


            if (SXYZ > 1.e-4) {

                double CT = BVS_Z / SXYZ;  // CT = cos(theta)
                double ST = SXY / SXYZ;    // ST = sin(theta)

                C[0][4][4] = 1.0;          // 2*SQ(pi) * (Y 00)

                if (_lmax_dft_ecp > 0) {
   
                    C[1][3][4] = ST * SP;       // 2*SQ(pi/3) * (Y 1-1)
                    C[1][4][4] = CT;            // 2*SQ(pi/3) * (Y 10)
                    C[1][5][4] = ST * CP;       // 2*SQ(pi/3) * (Y 11)                  Definition of (Z lm) :
                    if (_lmin_dft_ecp > 0) {
                        C[1][3][3] = CP;        // 2*SQ(pi/3) * (Z 1-1)                 (Z lm)  :=  ( 1 / sin(theta) ) * ( d (Y lm) / d phi )
                        C[1][3][5] = CT * SP;   // 2*SQ(pi/3) * (Y 1-1)'
                        C[1][4][3] = 0.0;
                        C[1][4][5] = -ST;       // 2*SQ(pi/3) * (Y 10)'                 Differentiation with respect to theta is denoted by an apostrophe:
                        C[1][5][3] = -SP;       // 2*SQ(pi/3) * (Z 11)
                        C[1][5][5] = CT * CP;   // 2*SQ(pi/3) * (Y 11)'                 f'  :=  d f / d theta
                     }

                }

                if (_lmax_dft_ecp > 1) {

                    C[2][2][4] = SQ3 * ST * ST * CP * SP;              // 2*SQ(pi/5) * (Y 2-2)
                    C[2][3][4] = SQ3 * CT * ST * SP;                   // 2*SQ(pi/5) * (Y 2-1)
                    C[2][4][4] = 1.5 * CT * CT - 0.5;                  // 2*SQ(pi/5) * (Y 20)
                    C[2][5][4] = SQ3 * CT * ST * CP;                   // 2*SQ(pi/5) * (Y 21)
                    C[2][6][4] = SQ3 * ST * ST * (CP * CP - .5);       // 2*SQ(pi/5) * (Y 22)
                    if (_lmin_dft_ecp > 0) {
                        C[2][2][3] = ST * (2.0 * CP * CP - 1.0);       // 2*SQ(pi/15) * (Z 2-2)
                        C[2][2][5] = 2.0 * CT * ST * CP * SP;          // 2*SQ(pi/15) * (Y 2-2)'
                        C[2][3][3] = CT * CP;                          // 2*SQ(pi/15) * (Z 2-1)
                        C[2][3][5] = (2.0 * CT * CT - 1.0) * SP;       // 2*SQ(pi/15) * (Y 2-1)'
                        C[2][4][3] = 0.0;
                        C[2][4][5] = -SQ3 * CT * ST;                   // 2*SQ(pi/15) * (Y 20)'
                        C[2][5][3] = -CT * SP;                         // 2*SQ(pi/15) * (Z 21)
                        C[2][5][5] = (2.0 * CT * CT - 1.0) * CP;       // 2*SQ(pi/15) * (Y 21)'
                        C[2][6][3] = -2.0 * ST * CP * SP;              // 2*SQ(pi/15) * (Z 22)
                        C[2][6][5] = CT * ST * (2.0 * CP * CP - 1.);   // 2*SQ(pi/15) * (Y 22)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[2][2][2] = CT * (2.0 * CP * CP - 1.0);       // 2*SQ(pi/15) * (Z 2-2)'
                        C[2][2][6] = (1.0 + CT * CT) * CP * SP;        // 2*SQ(pi/15) * ( (Y 2-2)'' + 3*(Y 2-2) )
                        C[2][3][2] = -ST * CP;                         // 2*SQ(pi/15) * (Z 2-1)'
                        C[2][3][6] = -CT * ST * SP;                    // 2*SQ(pi/15) * ( (Y 2-1)'' + 3*(Y 2-1) )
                        C[2][4][2] = 0.0;
                        C[2][4][6] = .5 * SQ3 * ST * ST;               // 2*SQ(pi/15) * ( (Y 20)'' + 3*(Y 20) )
                        C[2][5][2] = ST * SP;                          // 2*SQ(pi/15) * (Z 21)'
                        C[2][5][6] = -CT * ST * CP;                    // 2*SQ(pi/15) * ( (Y 21)'' + 3*(Y 21) )
                        C[2][6][2] = -2.0 * CT * CP * SP;              // 2*SQ(pi/15) * (Z 22)'
                        C[2][6][6] = (1. + CT * CT) * (CP * CP - .5);  // 2*SQ(pi/15) * ( (Y 22)'' + 3*(Y 22) )
                    }

                }

                if (_lmax_dft_ecp > 2) {

                    double f_phi, df_dphi;

                    f_phi = (4. * CP * CP - 1.) * SP;  // sin(3*phi)
                    C[3][1][4] = (.5*SQ5/SQ2) *         ST * ST * ST                   * f_phi;   // 2*SQ(pi/7) * (Y 3-3)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = 3. * (1. - 4. * SP * SP) * CP;
                        C[3][1][3] = (.25*SQ5/SQ3) *    ST * ST                        * df_dphi; // SQ((2*pi)/21) * (Z 3-3)
                        C[3][1][5] = (.25*SQ5/SQ3) *    3. * ST * ST * CT              * f_phi;   // SQ((2*pi)/21) * (Y 3-3)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[3][1][2] = (.5/(SQ2*SQ3)) *   2. * ST * CT                   * df_dphi; // 2*SQ(pi/105) * (Z 3-3)'
                        C[3][1][6] = (.5/(SQ2*SQ3)) *   3. * (1. + CT * CT) * ST       * f_phi;   // 2*SQ(pi/105) * ( (Y 3-3)'' + 6*(Y 3-3) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[3][1][1] = (1./6.) *          (.5 + 1.5 * CT * CT)           * df_dphi; // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Z 3-3) + (Z 3-3)'' )
                        C[3][1][7] = (1./6.) *          1.5 * (3. + CT * CT) * CT      * f_phi;   // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Y 3-3)' + ( (Y 3-3)'' + 6*(Y 3-3) )' )
                    }

                    f_phi = CP * SP;  // (1/2)*sin(2*phi)
                    C[3][2][4] = (SQ3*SQ5) *            ST * ST * CT                   * f_phi;   // 2*SQ(pi/7) * (Y 3-2)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = 1. - 2. * SP * SP;
                        C[3][2][3] = (SQ5/SQ2) *        ST * CT                        * df_dphi; // SQ((2*pi)/21) * (Z 3-2)
                        C[3][2][5] = (SQ5/SQ2) *        (3. * CT * CT - 1.) * ST       * f_phi;   // SQ((2*pi)/21) * (Y 3-2)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[3][2][2] =                    (1. - 2. * ST * ST)            * df_dphi; // 2*SQ(pi/105) * (Z 3-2)'
                        C[3][2][6] =                    (3. * CT * CT - 1.) * CT       * f_phi;   // 2*SQ(pi/105) * ( (Y 3-2)'' + 6*(Y 3-2) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[3][2][1] = (SQ2/SQ3) *        (-1.5) * ST * CT               * df_dphi; // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Z 3-2) + (Z 3-2)'' )
                        C[3][2][7] = (SQ2/SQ3) *        (-1.5) * (1. + CT * CT) * ST   * f_phi;   // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Y 3-2)' + ( (Y 3-2)'' + 6*(Y 3-2) )' )
                    }

                    f_phi = SP;
                    C[3][3][4] = (.5*SQ3/SQ2) *         (5. * CT * CT - 1.) * ST       * f_phi;   // 2*SQ(pi/7) * (Y 3-1)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = CP;
                        C[3][3][3] = .25 *              (5. * CT * CT - 1.)            * df_dphi; // SQ((2*pi)/21) * (Z 3-1)
                        C[3][3][5] = .25 *              (4. - 15. * ST * ST) * CT      * f_phi;   // SQ((2*pi)/21) * (Y 3-1)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[3][3][2] = (.5/(SQ2*SQ5)) *   (-10.) * ST * CT               * df_dphi; // 2*SQ(pi/105) * (Z 3-1)'
                        C[3][3][6] = (.5/(SQ2*SQ5)) *   (5. - 15. * CT * CT) * ST      * f_phi;   // 2*SQ(pi/105) * ( (Y 3-1)'' + 6*(Y 3-1) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[3][3][1] = (.5/(SQ3*SQ5)) *   7.5 * ST * ST                  * df_dphi; // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Z 3-1) + (Z 3-1)'' )
                        C[3][3][7] = (.5/(SQ3*SQ5)) *   7.5 * ST * ST * CT             * f_phi;   // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Y 3-1)' + ( (Y 3-1)'' + 6*(Y 3-1) )' )
                    }

                    C[3][4][4] = .5 *                   (5. * CT * CT - 3.) * CT;                 // 2*SQ(pi/7) * (Y 30)
                    if (_lmin_dft_ecp > 0) {
                        C[3][4][3] = 0.;
                        C[3][4][5] = (.5/(SQ2*SQ3)) *   (3. - 15. * CT * CT) * ST;                // SQ((2*pi)/21) * (Y 30)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[3][4][2] = 0.;
                        C[3][4][6] = (.5/(SQ3*SQ5)) *   15. * ST * ST * CT;                       // 2*SQ(pi/105) * ( (Y 30)'' + 6*(Y 30) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[3][4][1] = 0.;
                        C[3][4][7] = (SQ2/(6.*SQ5)) *   (-7.5) * ST * ST * ST;                    // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Y 30)' + ( (Y 30)'' + 6*(Y 30) )' )
                    }

                    f_phi = CP;
                    C[3][5][4] = (.5*SQ3/SQ2) *         (5. * CT * CT - 1.) * ST       * f_phi;   // 2*SQ(pi/7) * (Y 31)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = -SP; 
                        C[3][5][3] = .25 *              (5. * CT * CT - 1.)            * df_dphi; // SQ((2*pi)/21) * (Z 31)
                        C[3][5][5] = .25 *              (4. - 15. * ST * ST) * CT      * f_phi;   // SQ((2*pi)/21) * (Y 31)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[3][5][2] = (.5/(SQ2*SQ5)) *   (-10.) * ST * CT               * df_dphi; // 2*SQ(pi/105) * (Z 31)'
                        C[3][5][6] = (.5/(SQ2*SQ5)) *   (5. - 15. * CT * CT) * ST      * f_phi;   // 2*SQ(pi/105) * ( (Y 31)'' + 6*(Y 31) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[3][5][1] = (.5/(SQ3*SQ5)) *   7.5 * ST * ST                  * df_dphi; // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Z 31) + (Z 31)'' )
                        C[3][5][7] = (.5/(SQ3*SQ5)) *   7.5 * ST * ST * CT             * f_phi;   // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Y 31)' + ( (Y 31)'' + 6*(Y 31) )' )
                    }

                    f_phi = 1. - 2. * SP * SP;  // cos(2*phi)
                    C[3][6][4] = (.5*SQ3*SQ5) *         ST * ST * CT                   * f_phi;   // 2*SQ(pi/7) * (Y 32)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = - 4. * CP * SP;
                        C[3][6][3] = (.5*SQ5/SQ2) *     ST * CT                        * df_dphi; // SQ((2*pi)/21) * (Z 32)
                        C[3][6][5] = (.5*SQ5/SQ2) *     (3. * CT * CT - 1.) * ST       * f_phi;   // SQ((2*pi)/21) * (Y 32)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[3][6][2] = .5 *               (1. - 2. * ST * ST)            * df_dphi; // 2*SQ(pi/105) * (Z 32)'
                        C[3][6][6] = .5 *               (3. * CT * CT - 1.) * CT       * f_phi;   // 2*SQ(pi/105) * ( (Y 32)'' + 6*(Y 32) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[3][6][1] = (1./(SQ2*SQ3)) *   (-1.5) * ST * CT               * df_dphi; // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Z 32) + (Z 32)'' )
                        C[3][6][7] = (1./(SQ2*SQ3)) *   (-1.5) * (1. + CT * CT) * ST   * f_phi;   // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Y 32)' + ( (Y 32)'' + 6*(Y 32) )' )
                    }

                    f_phi = (1. - 4. * SP * SP) * CP;  // cos(3*phi)
                    C[3][7][4] = (.5*SQ5/SQ2) *         ST * ST * ST                   * f_phi;   // 2*SQ(pi/7) * (Y 33)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = 3. * (1. - 4. * CP * CP) * SP;
                        C[3][7][3] = (.25*SQ5/SQ3) *    ST * ST                        * df_dphi; // SQ((2*pi)/21) * (Z 33)
                        C[3][7][5] = (.25*SQ5/SQ3) *    3. * ST * ST * CT              * f_phi;   // SQ((2*pi)/21) * (Y 33)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[3][7][2] = (.5/(SQ2*SQ3)) *   2. * ST * CT                   * df_dphi; // 2*SQ(pi/105) * (Z 33)'
                        C[3][7][6] = (.5/(SQ2*SQ3)) *   3. * (1. + CT * CT) * ST       * f_phi;   // 2*SQ(pi/105) * ( (Y 33)'' + 6*(Y 33) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[3][7][1] = (1./6.) *          (.5 + 1.5 * CT * CT)           * df_dphi; // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Z 33) + (Z 33)'' )
                        C[3][7][7] = (1./6.) *          1.5 * (3. + CT * CT) * CT      * f_phi;   // (2/3)*SQ((2*pi)/35) * ( (5/2)*(Y 33)' + ( (Y 33)'' + 6*(Y 33) )' )
                    }

                }


                if (_lmax_dft_ecp > 3) {

                    double S2T = ST * ST;
                    double C2T = CT * CT;
                    double C2P = CP * CP;
                    double S2P = SP * SP;
                    double f_phi, df_dphi;

                    f_phi = CP * SP * (C2P - S2P);  // (1/4)*sin(4*phi)
                    C[4][0][4] =                          (.5*SQ5*SQ7) *  S2T * S2T                  * f_phi;   // (2/3)*SQ(pi) * (Y 4-4)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = (C2P - S2P) * (C2P - S2P) - 4. * C2P * S2P;
                        C[4][0][3] = (1./(SQ2*SQ5))     * (.5*SQ5*SQ7) *  S2T * ST                   * df_dphi; // (1/3)*SQ((2*pi)/5) * (Z 4-4)
                        C[4][0][5] = (1./(SQ2*SQ5))     * (.5*SQ5*SQ7) *  4. * S2T * ST * CT         * f_phi;   // (1/3)*SQ((2*pi)/5) * (Y 4-4)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[4][0][2] = (1./(3.*SQ5))      * (.5*SQ5*SQ7) *  3. * S2T * CT              * df_dphi; // (2/9)*SQ(pi/5) * (Z 4-4)'
                        C[4][0][6] = (1./(3.*SQ5))      * (.5*SQ5*SQ7) *  6. * (C2T + 1.) * S2T      * f_phi;   // (2/9)*SQ(pi/5) * ( (Y 4-4)'' + 10*(Y 4-4) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[4][0][1] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5*SQ7) *  (4.5 * C2T + 1.5) * ST     * df_dphi; // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Z 4-4) + (Z 4-4)'' )
                        C[4][0][7] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5*SQ7) *  6. * (C2T + 3.) * ST * CT  * f_phi;   // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Y 4-4)' + ( (Y 4-4)'' + 10*(Y 4-4) )' )
                    }

                    f_phi = (4. * C2P - 1.) * SP;  // sin(3*phi)
                    C[4][1][4] =                          (.5*SQ5*SQ7/SQ2) *  S2T * ST * CT                 * f_phi;   // (2/3)*SQ(pi) * (Y 4-3)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = 3. * (1. - 4. * SP * SP) * CP;
                        C[4][1][3] = (1./(SQ2*SQ5))     * (.5*SQ5*SQ7/SQ2) *  S2T * CT                      * df_dphi; // (1/3)*SQ((2*pi)/5) * (Z 4-3)
                        C[4][1][5] = (1./(SQ2*SQ5))     * (.5*SQ5*SQ7/SQ2) *  (4. * C2T - 1.) * S2T         * f_phi;   // (1/3)*SQ((2*pi)/5) * (Y 4-3)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[4][1][2] = (1./(3.*SQ5))      * (.5*SQ5*SQ7/SQ2) *  (3. * C2T - 1.) * ST          * df_dphi; // (2/9)*SQ(pi/5) * (Z 4-3)'
                        C[4][1][6] = (1./(3.*SQ5))      * (.5*SQ5*SQ7/SQ2) *  6. * C2T * CT * ST            * f_phi;   // (2/9)*SQ(pi/5) * ( (Y 4-3)'' + 10*(Y 4-3) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[4][1][1] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5*SQ7/SQ2) *  (2. - 4.5 * S2T) * CT         * df_dphi; // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Z 4-3) + (Z 4-3)'' )
                        C[4][1][7] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5*SQ7/SQ2) *  (6. * C2T * C2T - 4.5 * S2T)  * f_phi;   // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Y 4-3)' + ( (Y 4-3)'' + 10*(Y 4-3) )' )
                    }

                    f_phi = CP * SP;  // (1/2)*sin(2*phi)
                    C[4][2][4] =                          (.5*SQ5) *  (7. * C2T - 1.) * S2T               * f_phi;   // (2/3)*SQ(pi) * (Y 4-2)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = 1. - 2. * SP * SP;
                        C[4][2][3] = (1./(SQ2*SQ5))     * (.5*SQ5) *  (7. * C2T - 1.) * ST                * df_dphi; // (1/3)*SQ((2*pi)/5) * (Z 4-2)
                        C[4][2][5] = (1./(SQ2*SQ5))     * (.5*SQ5) *  (14. * (C2T - S2T) - 2.) * ST * CT  * f_phi;   // (1/3)*SQ((2*pi)/5) * (Y 4-2)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[4][2][2] = (1./(3.*SQ5))      * (.5*SQ5) *  (6. - 21. * S2T) * CT               * df_dphi; // (2/9)*SQ(pi/5) * (Z 4-2)'
                        C[4][2][6] = (1./(3.*SQ5))      * (.5*SQ5) *  6. * (1. + C2T - 7. * S2T * C2T)    * f_phi;   // (2/9)*SQ(pi/5) * ( (Y 4-2)'' + 10*(Y 4-2) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[4][2][1] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5) *  10.5 * (1. - 3. * C2T) * ST         * df_dphi; // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Z 4-2) + (Z 4-2)'' )
                        C[4][2][7] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5) *  (-42.) * C2T * CT * ST              * f_phi;   // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Y 4-2)' + ( (Y 4-2)'' + 10*(Y 4-2) )' )
                    }

                    f_phi = SP;
                    C[4][3][4] =                          (.5*SQ5/SQ2) *  (7. * C2T - 3.) * ST * CT       * f_phi;   // (2/3)*SQ(pi) * (Y 4-1)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = CP;
                        C[4][3][3] = (1./(SQ2*SQ5))     * (.5*SQ5/SQ2) *  (7. * C2T - 3.) * CT            * df_dphi; // (1/3)*SQ((2*pi)/5) * (Z 4-1)
                        C[4][3][5] = (1./(SQ2*SQ5))     * (.5*SQ5/SQ2) *  (3. + (1. - 28. * S2T) * C2T)   * f_phi;   // (1/3)*SQ((2*pi)/5) * (Y 4-1)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[4][3][2] = (1./(3.*SQ5))      * (.5*SQ5/SQ2) *  3. * (1. - 7. * C2T) * ST       * df_dphi; // (2/9)*SQ(pi/5) * (Z 4-1)'
                        C[4][3][6] = (1./(3.*SQ5))      * (.5*SQ5/SQ2) *  6. * (7. * S2T - 3.) * ST * CT  * f_phi;   // (2/9)*SQ(pi/5) * ( (Y 4-1)'' + 10*(Y 4-1) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[4][3][1] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5/SQ2) *  .5 * 63. * S2T * CT             * df_dphi; // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Z 4-1) + (Z 4-1)'' )
                        C[4][3][7] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5/SQ2) *  21. * (1.5 - 2. * S2T) * S2T    * f_phi;   // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Y 4-1)' + ( (Y 4-1)'' + 10*(Y 4-1) )' )
                    }


                    C[4][4][4] =                          .125 *       ((35. * C2T - 30.) * C2T + 3.);         // (2/3)*SQ(pi) * (Y 40)
                    if (_lmin_dft_ecp > 0) {
                        C[4][4][3] = 0.;
                        C[4][4][5] = (1./(SQ2*SQ5))     * .125 *       20. * (3.  - 7. * C2T) * CT * ST;     // (1/3)*SQ((2*pi)/5) * (Y 40)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[4][4][2] = 0.;
                        C[4][4][6] = (1./(3.*SQ5))      * .125 *       30. * (7. * C2T - 1.) * S2T;          // (2/9)*SQ(pi/5) * ( (Y 40)'' + 10*(Y 40) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[4][4][1] = 0.;
                        C[4][4][7] = (SQ2/(3.*SQ5*SQ7)) * .125 *       (-210.) * S2T * ST * CT;              // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Y 40)' + ( (Y 40)'' + 10*(Y 40) )' )
                    }

                    f_phi = CP;
                    C[4][5][4] =                          (.5*SQ5/SQ2) *  (7. * C2T - 3.) * ST * CT       * f_phi;   // (2/3)*SQ(pi) * (Y 41)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = -SP;
                        C[4][5][3] = (1./(SQ2*SQ5))     * (.5*SQ5/SQ2) *  (7. * C2T - 3.) * CT            * df_dphi; // (1/3)*SQ((2*pi)/5) * (Z 41)
                        C[4][5][5] = (1./(SQ2*SQ5))     * (.5*SQ5/SQ2) *  (3. + (1. - 28. * S2T) * C2T)   * f_phi;   // (1/3)*SQ((2*pi)/5) * (Y 41)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[4][5][2] = (1./(3.*SQ5))      * (.5*SQ5/SQ2) *  3. * (1. - 7. * C2T) * ST       * df_dphi; // (2/9)*SQ(pi/5) * (Z 41)'
                        C[4][5][6] = (1./(3.*SQ5))      * (.5*SQ5/SQ2) *  6. * (7. * S2T - 3.) * ST * CT  * f_phi;   // (2/9)*SQ(pi/5) * ( (Y 41)'' + 10*(Y 41) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[4][5][1] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5/SQ2) *  .5 * 63. * S2T * CT             * df_dphi; // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Z 41) + (Z 41)'' )
                        C[4][5][7] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5/SQ2) *  21. * (1.5 - 2. * S2T) * S2T    * f_phi;   // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Y 41)' + ( (Y 41)'' + 10*(Y 41) )' )
                    }

                    f_phi = C2P - S2P;  // cos(2*phi)
                    C[4][6][4] =                          (.25*SQ5) *  (7. * C2T - 1.) * S2T               * f_phi;   // (2/3)*SQ(pi) * (Y 42)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = - 4. * CP * SP;
                        C[4][6][3] = (1./(SQ2*SQ5))     * (.25*SQ5) *  (7. * C2T - 1.) * ST                * df_dphi; // (1/3)*SQ((2*pi)/5) * (Z 42)
                        C[4][6][5] = (1./(SQ2*SQ5))     * (.25*SQ5) *  (14. * (C2T - S2T) - 2.) * ST * CT  * f_phi;   // (1/3)*SQ((2*pi)/5) * (Y 42)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[4][6][2] = (1./(3.*SQ5))      * (.25*SQ5) *  (6. - 21. * S2T) * CT               * df_dphi; // (2/9)*SQ(pi/5) * (Z 42)'
                        C[4][6][6] = (1./(3.*SQ5))      * (.25*SQ5) *  6. * (1. + C2T - 7. * S2T * C2T)    * f_phi;   // (2/9)*SQ(pi/5) * ( (Y 42)'' + 10*(Y 42) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[4][6][1] = (SQ2/(3.*SQ5*SQ7)) * (.25*SQ5) *  10.5 * (1. - 3. * C2T) * ST         * df_dphi; // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Z 42) + (Z 42)'' )
                        C[4][6][7] = (SQ2/(3.*SQ5*SQ7)) * (.25*SQ5) *  (-42.) * C2T * CT * ST              * f_phi;   // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Y 42)' + ( (Y 42)'' + 10*(Y 42) )' )
                    }

                    f_phi = (1. - 4. * S2P) * CP;  // cos(3*phi)
                    C[4][7][4] =                          (.5*SQ5*SQ7/SQ2) *  S2T * ST * CT                 * f_phi;   // (2/3)*SQ(pi) * (Y 43)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = 3. * (1. - 4. * CP * CP) * SP;
                        C[4][7][3] = (1./(SQ2*SQ5))     * (.5*SQ5*SQ7/SQ2) *  S2T * CT                      * df_dphi; // (1/3)*SQ((2*pi)/5) * (Z 43)
                        C[4][7][5] = (1./(SQ2*SQ5))     * (.5*SQ5*SQ7/SQ2) *  (4. * C2T - 1.) * S2T         * f_phi;   // (1/3)*SQ((2*pi)/5) * (Y 43)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[4][7][2] = (1./(3.*SQ5))      * (.5*SQ5*SQ7/SQ2) *  (3. * C2T - 1.) * ST          * df_dphi; // (2/9)*SQ(pi/5) * (Z 43)'
                        C[4][7][6] = (1./(3.*SQ5))      * (.5*SQ5*SQ7/SQ2) *  6. * C2T * CT * ST            * f_phi;   // (2/9)*SQ(pi/5) * ( (Y 43)'' + 10*(Y 43) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[4][7][1] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5*SQ7/SQ2) *  (2. - 4.5 * S2T) * CT         * df_dphi; // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Z 43) + (Z 43)'' )
                        C[4][7][7] = (SQ2/(3.*SQ5*SQ7)) * (.5*SQ5*SQ7/SQ2) *  (6. * C2T * C2T - 4.5 * S2T)  * f_phi;   // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Y 43)' + ( (Y 43)'' + 10*(Y 43) )' )
                    }

                    f_phi = (C2P - S2P) * (C2P - S2P) - 4. * C2P * S2P;  // cos(4*phi)
                    C[4][8][4] =                          (.125*SQ5*SQ7) *  S2T * S2T                  * f_phi;   // (2/3)*SQ(pi) * (Y 44)
                    if (_lmin_dft_ecp > 0) {
                        df_dphi = 16. * CP * SP * (S2P - C2P);
                        C[4][8][3] = (1./(SQ2*SQ5))     * (.125*SQ5*SQ7) *  S2T * ST                   * df_dphi; // (1/3)*SQ((2*pi)/5) * (Z 44)
                        C[4][8][5] = (1./(SQ2*SQ5))     * (.125*SQ5*SQ7) *  4. * S2T * ST * CT         * f_phi;   // (1/3)*SQ((2*pi)/5) * (Y 44)'
                    }
                    if (_lmin_dft_ecp > 1) {
                        C[4][8][2] = (1./(3.*SQ5))      * (.125*SQ5*SQ7) *  3. * S2T * CT              * df_dphi; // (2/9)*SQ(pi/5) * (Z 44)'
                        C[4][8][6] = (1./(3.*SQ5))      * (.125*SQ5*SQ7) *  6. * (C2T + 1.) * S2T      * f_phi;   // (2/9)*SQ(pi/5) * ( (Y 44)'' + 10*(Y 44) )
                    }
                    if (_lmin_dft_ecp > 2) {
                        C[4][8][1] = (SQ2/(3.*SQ5*SQ7)) * (.125*SQ5*SQ7) *  (4.5 * C2T + 1.5) * ST     * df_dphi; // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Z 44) + (Z 44)'' )
                        C[4][8][7] = (SQ2/(3.*SQ5*SQ7)) * (.125*SQ5*SQ7) *  6. * (C2T + 3.) * ST * CT  * f_phi;   // (2/9)*SQ((2*pi)/35) * ( (9/2)*(Y 44)' + ( (Y 44)'' + 10*(Y 44) )' )
                    }

                }


                for (index I = 0; I < _nsph; I++) {
                    for (index L = 0; L <= _lmax_dft; L++) {
                        for ( index M = 4 - L; M <= 4 + L; M++) {
                            BLC[I][L][M] = 0.0;
                            for (index M1 = 4 - L; M1 <= 4 + L; M1++) {
                                BLC[I][L][M] += BLM[I][L][M1] * C[L][M1][M];
                            }
                        }
                    }

                }


            } else {

                for (index L = 0; L <= _lmax_dft_ecp; L++) {
                    for (index M = 4 - L; M <= 4 + L; M++) {
                        C[L][M][M] = 1.;
                    }
                }


                for (index I = 0; I < _nsph; I++) {
                    for (index L = 0; L <= _lmax_dft; L++) {
                        for ( index M = 4 - L; M <= 4 + L; M++) {
                            BLC[I][L][M] = BLM[I][L][M];
                        }
                    }
                }

            }

            return;
        } // getBLMCOF




}}

