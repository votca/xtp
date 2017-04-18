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
// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>
//for libxc
#include <votca/xtp/votca_config.h>

#include <votca/xtp/numerical_integrations.h>
#include <boost/math/constants/constants.hpp>
#include <votca/xtp/radial_euler_maclaurin_rule.h>
#include <votca/xtp/sphere_lebedev_rule.h>
#include <votca/xtp/aoshell.h>
#include <votca/tools/constants.h>

#ifdef LIBXC
#include <xc.h>
#endif
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/exchange_correlation.h>
#include <fstream>
#include <boost/timer/timer.hpp>
#include <boost/algorithm/string.hpp>
#include <votca/xtp/vxc_functionals.h>
#include <iterator>
#include <string>




namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        double NumericalIntegration::getExactExchange(const string _functional) {
#ifdef LIBXC            

            double exactexchange = 0.0;
            Vxc_Functionals map;
            std::vector<string> strs;

            boost::split(strs, _functional, boost::is_any_of(" "));
            if (strs.size() > 2) {
                throw std::runtime_error("Too many functional names");
            } else if (strs.size() < 1) {
                throw std::runtime_error("Specify at least one funcitonal");
            }

            for (unsigned i = 0; i < strs.size(); i++) {

                int func_id = map.getID(strs[i]);
                if (func_id < 0) {
                    exactexchange = 0.0;
                    break;
                }
                xc_func_type func;
                if (xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0) {
                    fprintf(stderr, "Functional '%d' not found\n", func_id);
                    exit(1);
                }
                if (exactexchange > 0 && func.cam_alpha > 0) {
                    throw std::runtime_error("You have specified two functionals with exact exchange");
                }
                exactexchange += func.cam_alpha;



            }
            return exactexchange;

#else
            return 0.0;
#endif

        }

        ub::matrix<double> NumericalIntegration::IntegrateExternalPotential_Atomblock(AOBasis* basis, std::vector<double> Potentialvalues) {
            if (_significant_atoms.size() < 1) {
                throw runtime_error("NumericalIntegration::IntegrateExternalPotential_Atomblock:significant atoms not found yet.");
            }
            ub::matrix<double> ExternalMat = ub::zero_matrix<double>(basis->_AOBasisSize, basis->_AOBasisSize);

            // parallelization: distribute over threads inside one atom
            int nthreads = 1;
#ifdef _OPENMP
            nthreads = omp_get_max_threads();
#endif

            // separate storage for each thread
            std::vector< ub::matrix<double> > expot_thread;

            for (int i_thread = 0; i_thread < nthreads; i_thread++) {
                expot_thread.push_back(ub::zero_matrix<double>(basis->_AOBasisSize, basis->_AOBasisSize));
            }

            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
                // for each point in atom grid

                // number of points in this atomgrid
                int atom_points = _grid[i].size();
                // divide among threads
                int atom_points_per_thread = atom_points / nthreads;
                std::vector<int> _thread_start;
                std::vector<int> _thread_stop;
                for (int i_thread = 0; i_thread < nthreads; i_thread++) {
                    _thread_start.push_back(i_thread * atom_points_per_thread);
                    _thread_stop.push_back((i_thread + 1) * atom_points_per_thread);
                }
                // final stop must be size
                _thread_stop[nthreads - 1] = atom_points;


#pragma omp parallel for
                for (int i_thread = 0; i_thread < nthreads; i_thread++) {
                    for (int j = _thread_start[i_thread]; j < _thread_stop[i_thread]; j++) {

                        // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )

                        ub::matrix<double> AOgrid = ub::zero_matrix<double>(1, basis->_AOBasisSize); // TRY MORE USEFUL DATA

                        // for each significant atom for this grid point
                        for (unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size(); sigrow++) {

                            // this atom
                            int rowatom = _significant_atoms[i][j][sigrow];

                            // for each shell in this atom
                            for (unsigned ishell = 0; ishell < _atomshells[rowatom].size(); ishell++) {

                                //   boost::timer::cpu_times tstartshells = cpu_t.elapsed();
                                AOShellIterator _row = _atomshells[rowatom][ishell];
                                // for density, fill sub-part of AOatgrid
                                //ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 1);
                                ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, 0, 1, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());
                                // (*_row)->EvalAOspace(_AOgridsub, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);


                                (*_row)->EvalAOspace(_AOgridsub, _grid[i][j].grid_pos);


                            } // shell in atom


                        } // row shells 


                        // yeah storin potential values in a vector is weird but I did not want to cram it into gridpoint, because that will blow the structure more than necessary
                        ub::matrix<double> _addExt = 0.5 * _grid[i][j].grid_weight * AOgrid * Potentialvalues[i * _grid[i].size() + j];

                        // combine/sum atom-block wise, only trigonal part, symmetrize later
                        // for each significant atom for this grid point
                        // parallelization only accesses atomblock information (_addXC, AOgrid -> XCmatblock), so no trouble with shared memory access )
                        // #pragma omp parallel for
                        for (unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size(); sigrow++) {

                            // this atom
                            int rowatom = _significant_atoms[i][j][sigrow];

                            ub::matrix_range< ub::matrix<double> > _rowExt = ub::subrange(_addExt, 0, 1, _startIdx[rowatom], _startIdx[rowatom] + _blocksize[rowatom]);

                            for (unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size(); sigcol++) {
                                int colatom = _significant_atoms[i][j][sigcol];
                                // if (colatom > rowatom) break;

                                ub::matrix_range< ub::matrix<double> > _AOcol = ub::subrange(AOgrid, 0, 1, _startIdx[colatom], _startIdx[colatom] + _blocksize[colatom]);

                                // update block reference of XCMAT
                                ub::matrix_range<ub::matrix<double> > _expotmatblock = ub::subrange(expot_thread[i_thread], _startIdx[rowatom], _startIdx[rowatom] + _blocksize[rowatom], _startIdx[colatom], _startIdx[colatom] + _blocksize[colatom]);
                                //_XCmatblock += ub::prod( _rowXC, ub::trans(_AOcol)  );
                                _expotmatblock += ub::prod(ub::trans(_rowExt), _AOcol);

                                // update the other block

                            } // significant col
                        } // significant row 

                    } // j: for each point in atom grid
                }// each thread
            } // i: for each atom grid


            // sum thread matrices
            for (int i_thread = 0; i_thread < nthreads; i_thread++) {
#pragma omp parallel for
                for (unsigned _i = 0; _i < ExternalMat.size1(); _i++) {
                    //for (int _j = 0; _j <= _i; _j++) {
                    for (unsigned _j = 0; _j < ExternalMat.size2(); _j++) {
                        ExternalMat(_i, _j) += expot_thread[i_thread](_i, _j);
                    }
                }
            }

            ExternalMat += ub::trans(ExternalMat);

            return ExternalMat;

        }

        ub::matrix<double> NumericalIntegration::IntegrateVXC_Atomblock(const ub::matrix<double>& _density_matrix, AOBasis* basis, const string _functional) {
            EXC = 0;
            if (_significant_atoms.size() < 1) {
                throw runtime_error("NumericalIntegration::IntegrateVXC_Atomblock:significant atoms not found yet.");
            }
            // TODO: switch XC functionals implementation from LIBXC to base own calculation
            ExchangeCorrelation _xc;
            Vxc_Functionals map;
            std::vector<string> strs;
            boost::split(strs, _functional, boost::is_any_of(" "));
            int xfunc_id = 0;

#ifdef LIBXC
            bool _use_votca = false;
            bool _use_separate = false;
            int cfunc_id = 0;

            if (strs.size() == 1) {
                xfunc_id = map.getID(strs[0]);
                if (xfunc_id < 0) _use_votca = true;
            }
            else if (strs.size() == 2) {
                cfunc_id = map.getID(strs[0]);
                xfunc_id = map.getID(strs[1]);
                _use_separate = true;
            } else {
                cout << "LIBXC " << strs.size() << endl;
                throw std::runtime_error("With LIBXC. Please specify one combined or an exchange and a correlation functionals");

            }
            xc_func_type xfunc; // handle for exchange functional
            xc_func_type cfunc; // handle for correlation functional
            if (!_use_votca) {
                if (xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED) != 0) {
                    fprintf(stderr, "Functional '%d' not found\n", xfunc_id);
                    exit(1);
                }

                xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED);
                if (xfunc.info->kind != 2 && !_use_separate) {
                    throw std::runtime_error("Your functional misses either correlation or exchange, please specify another functional, separated by whitespace");
                }

                if (_use_separate) {
                    if (xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED) != 0) {
                        fprintf(stderr, "Functional '%d' not found\n", cfunc_id);
                        exit(1);
                    }
                    xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED);
                    xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED);
                    if ((xfunc.info->kind + cfunc.info->kind) != 1) {
                        throw std::runtime_error("Your functionals are not one exchange and one correlation");
                    }
                }
            }
#else
            if (strs.size() == 1) {
                xfunc_id = map.getID(strs[0]);
            }
            else {
                throw std::runtime_error("Running without LIBXC, Please specify one combined or an exchange and a correlation functionals");
            }
#endif

            //split dmat into atomsize protions so that access is faster later on
#pragma omp parallel for
            for (unsigned rowatom = 0; rowatom < _grid.size(); rowatom++) {
                for (unsigned colatom = 0; colatom <= rowatom; colatom++) {

                    dmat_vector[rowatom][colatom] = ub::subrange(_density_matrix, _startIdx[colatom], _startIdx[colatom] + _blocksize[colatom], _startIdx[rowatom], _startIdx[rowatom] + _blocksize[rowatom]);
                }
            }



            // parallelization: distribute over threads inside one atom
            int nthreads = 1;
#ifdef _OPENMP
            nthreads = omp_get_max_threads();
#endif

            // separate storage for each thread
            //same as for dmat for vxc mat
            std::vector< ub::matrix<double> > XCMAT_thread;
            std::vector<double> EXC_thread;
            for (int i_thread = 0; i_thread < nthreads; i_thread++) {
                EXC_thread.push_back(0.0);

            }
#pragma omp parallel for
            for (int i_thread = 0; i_thread < nthreads; i_thread++) {

                for (unsigned rowatom = 0; rowatom < _grid.size(); rowatom++) {
                    for (unsigned colatom = 0; colatom < _grid.size(); colatom++) {
                        xcmat_vector_thread[i_thread][rowatom][colatom] = ub::zero_matrix<double>(_blocksize[rowatom], _blocksize[colatom]);
                    }
                }
            }
#pragma omp parallel for
            for (unsigned rowatom = 0; rowatom < _grid.size(); rowatom++) {
                for (unsigned colatom = 0; colatom < _grid.size(); colatom++) {
                    xcmat_vector[rowatom][colatom] = ub::zero_matrix<double>(_blocksize[rowatom], _blocksize[colatom]);
                }
            }

            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
                // for each point in atom grid

                // number of points in this atomgrid
                int atom_points = _grid[i].size();
                // divide among threads
                int atom_points_per_thread = atom_points / nthreads;
                std::vector<int> _thread_start;
                std::vector<int> _thread_stop;
                for (int i_thread = 0; i_thread < nthreads; i_thread++) {
                    _thread_start.push_back(i_thread * atom_points_per_thread);
                    _thread_stop.push_back((i_thread + 1) * atom_points_per_thread);
                }
                // final stop must be size
                _thread_stop[nthreads - 1] = atom_points;


#pragma omp parallel for
                for (int i_thread = 0; i_thread < nthreads; i_thread++) {
                    for (int j = _thread_start[i_thread]; j < _thread_stop[i_thread]; j++) {

                        // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )
                        ub::matrix<double> AOgrid = ub::zero_matrix<double>(1, basis->_AOBasisSize); // TRY MORE USEFUL DATA
                        // get value of density gradient at each gridpoint

                        ub::matrix<double> gradAOgrid = ub::zero_matrix<double>(3, basis->_AOBasisSize); // for Gradients of AOs

                        ub::matrix<double> rho_mat = ub::zero_matrix<double>(1, 1);
                        //ub::matrix<double> grad_rho = ub::zero_matrix<double>(3,1);

                        ub::matrix<double> grad_rho = ub::zero_matrix<double>(1, 3);
                        // evaluate AO Functions for all shells, NOW BLOCKWISE

                        // for each significant atom for this grid point
                        for (unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size(); sigrow++) {

                            // this atom
                            int rowatom = _significant_atoms[i][j][sigrow];

                            // for each shell in this atom
                            for (unsigned ishell = 0; ishell < _atomshells[rowatom].size(); ishell++) {


                                AOShellIterator _row = _atomshells[rowatom][ishell];
                                // for density, fill sub-part of AOatgrid

                                ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, 0, 1, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());


                                // gradient of density

                                ub::matrix_range< ub::matrix<double> > _gradAO = ub::subrange(gradAOgrid, 0, 3, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());

                                (*_row)->EvalAOspace(_AOgridsub, _gradAO, _grid[i][j].grid_pos);


                            } // shell in atom

                            ub::matrix<double> _temp = ub::zero_matrix<double>(1, _blocksize[rowatom]);
                            ub::matrix<double> _tempgrad = ub::zero_matrix<double>(3, _blocksize[rowatom]);

                            ub::matrix_range< ub::matrix<double> > _AOgridrow = ub::subrange(AOgrid, 0, 1, _startIdx[rowatom], _startIdx[rowatom] + _blocksize[rowatom]);

                            // for each atom
                            // for all significant atoms of triangular matrix
                            for (unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size(); sigcol++) {

                                int colatom = _significant_atoms[i][j][sigcol];
                                if (colatom > rowatom) break;

                                // get the already calculated AO values

                                ub::matrix_range< ub::matrix<double> > _AOgridcol = ub::subrange(AOgrid, 0, 1, _startIdx[colatom], _startIdx[colatom] + _blocksize[colatom]);
                                ub::matrix_range< ub::matrix<double> > _gradAOgridcol = ub::subrange(gradAOgrid, 0, 3, _startIdx[colatom], _startIdx[colatom] + _blocksize[colatom]);


                                const ub::matrix<double> & DMAT_here = dmat_vector[rowatom][colatom];

                                if (colatom == rowatom) {
                                    _temp += 0.5 * ub::prod(_AOgridcol, DMAT_here);
                                    _tempgrad += 0.5 * ub::prod(_gradAOgridcol, DMAT_here);
                                } else {

                                    _temp += ub::prod(_AOgridcol, DMAT_here);
                                    _tempgrad += ub::prod(_gradAOgridcol, DMAT_here);
                                }

                            } //col shells


                            ub::matrix_range< ub::matrix<double> > _gradAOgridrow = ub::subrange(gradAOgrid, 0, 3, _startIdx[rowatom], _startIdx[rowatom] + _blocksize[rowatom]);


                            rho_mat += ub::prod(_temp, ub::trans(_AOgridrow));
                            grad_rho += ub::prod(_temp, ub::trans(_gradAOgridrow)) + ub::prod(_AOgridrow, ub::trans(_tempgrad));

                        } // row shells 

                        double rho = 2.0 * rho_mat(0, 0);


                        if (rho < 1.e-15) continue; // skip the rest, if density is very small
                        grad_rho = 2.0 * grad_rho;

                        // get XC for this density_at_grid
                        double f_xc; // E_xc[n] = int{n(r)*eps_xc[n(r)] d3r} = int{ f_xc(r) d3r }
                        double df_drho; // v_xc_rho(r) = df/drho
                        double df_dsigma; // df/dsigma ( df/dgrad(rho) = df/dsigma * dsigma/dgrad(rho) = df/dsigma * 2*grad(rho))

#ifdef LIBXC                   
                        if (_use_votca) {
#endif                  
                            cout << "Warning: VOTCA_PBE does give correct Vxc but incorrect E_xc" << endl;
                            _xc.getXC(xfunc_id, rho, grad_rho(0, 0), grad_rho(0, 1), grad_rho(0, 2), f_xc, df_drho, df_dsigma);
#ifdef LIBXC
                        }// evaluate via LIBXC, if compiled, otherwise, go via own implementation

                        else {


                            double sigma = ub::prod(grad_rho, ub::trans(grad_rho))(0, 0);

                            double exc[1];
                            double vsigma[1]; // libxc 
                            double vrho[1]; // libxc df/drho
                            switch (xfunc.info->family) {
                                case XC_FAMILY_LDA:
                                    xc_lda_exc_vxc(&xfunc, 1, &rho, exc, vrho);
                                    break;
                                case XC_FAMILY_GGA:
                                case XC_FAMILY_HYB_GGA:
                                    xc_gga_exc_vxc(&xfunc, 1, &rho, &sigma, exc, vrho, vsigma);
                                    break;
                            }
                            f_xc = exc[0];
                            df_drho = vrho[0];
                            df_dsigma = vsigma[0];
                            if (_use_separate) {
                                // via libxc correlation part only
                                switch (cfunc.info->family) {
                                    case XC_FAMILY_LDA:
                                        xc_lda_exc_vxc(&cfunc, 1, &rho, exc, vrho);
                                        break;
                                    case XC_FAMILY_GGA:
                                    case XC_FAMILY_HYB_GGA:
                                        xc_gga_exc_vxc(&cfunc, 1, &rho, &sigma, exc, vrho, vsigma);
                                        break;
                                }

                                f_xc += exc[0];
                                df_drho += vrho[0];
                                df_dsigma += vsigma[0];
                            }
                        }
#endif

                        ub::matrix<double> _addXC = _grid[i][j].grid_weight * df_drho * AOgrid * 0.5;

                        _addXC += 2.0 * df_dsigma * _grid[i][j].grid_weight * ub::prod(grad_rho, gradAOgrid);

                        // Exchange correlation energy
                        EXC_thread[i_thread] += _grid[i][j].grid_weight * rho * f_xc;

                        // combine/sum atom-block wise, only trigonal part, symmetrize later
                        // for each significant atom for this grid point
                        // parallelization only accesses atomblock information (_addXC, AOgrid -> XCmatblock), so no trouble with shared memory access )
                        // #pragma omp parallel for
                        for (unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size(); sigrow++) {

                            // this atom
                            int rowatom = _significant_atoms[i][j][sigrow];

                            const ub::matrix_range< ub::matrix<double> > _rowXC = ub::subrange(_addXC, 0, 1, _startIdx[rowatom], _startIdx[rowatom] + _blocksize[rowatom]);


                            std::vector< ub::matrix<double> >& _XCmatblock = xcmat_vector_thread[i_thread][rowatom];
                            for (unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size(); sigcol++) {
                                int colatom = _significant_atoms[i][j][sigcol];

                                const ub::matrix_range< ub::matrix<double> > _AOcol = ub::subrange(AOgrid, 0, 1, _startIdx[colatom], _startIdx[colatom] + _blocksize[colatom]);
                                _XCmatblock[colatom] += ub::prod(ub::trans(_rowXC), _AOcol);

                            } // significant col
                        } // significant row 

                    } // j: for each point in atom grid
                }// each thread
            } // i: for each atom grid


            // sum thread matrices
            for (int i_thread = 0; i_thread < nthreads; i_thread++) {
                EXC += EXC_thread[i_thread];

            }

            for (int i_thread = 0; i_thread < nthreads; i_thread++) {
#pragma omp parallel for
                for (unsigned _i = 0; _i < xcmat_vector.size(); _i++) {
                    for (unsigned _j = 0; _j < xcmat_vector[_i].size(); _j++) {


                        xcmat_vector[_i][_j] += xcmat_vector_thread[i_thread][_i][_j];
                    }
                }
            }


            ub::matrix<double> XCMAT = ub::zero_matrix<double>(basis->_AOBasisSize, basis->_AOBasisSize);

#pragma omp parallel for
            for (unsigned rowatom = 0; rowatom < xcmat_vector.size(); rowatom++) {
                for (unsigned colatom = 0; colatom < xcmat_vector[rowatom].size(); colatom++) {



                    ub::subrange(XCMAT, _startIdx[rowatom], _startIdx[rowatom] + _blocksize[rowatom], _startIdx[colatom], _startIdx[colatom] + _blocksize[colatom])
                            = xcmat_vector[rowatom][colatom];

                }
            }

            XCMAT += ub::trans(XCMAT);

            return XCMAT;
        }

        double NumericalIntegration::IntegratePotential(const vec& rvector) {

            double result = 0.0;

            if (density_set) {
                for (unsigned i = 0; i < _grid.size(); i++) {
                    for (unsigned j = 0; j < _grid[i].size(); j++) {
                        double dist = abs((_grid[i][j].grid_pos - rvector));
                        result -= _grid[i][j].grid_weight * _grid[i][j].grid_density / dist;
                    }
                }
            }
            else {
                throw std::runtime_error("Density not calculated");
            }

            return result;
        }

        void NumericalIntegration::FindsignificantAtoms(AOBasis* basis) {

            int _atomindex = 0;
            int _Idx = 0;
            int _size = 0;

            for (vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++) {

                if ((*_row)->getIndex() == _atomindex) {

                    _singleatom.push_back(_row);
                    _size += (*_row)->getNumFunc();

                } else {

                    // append _singleatom to _atomshells
                    _atomshells.push_back(_singleatom);
                    _startIdx.push_back(_Idx);
                    _blocksize.push_back(_size);
                    // reset _singleatom
                    _singleatom.clear();
                    _size = (*_row)->getNumFunc();
                    _Idx = (*_row)->getStartIndex();
                    _singleatom.push_back(_row);
                    _atomindex = (*_row)->getIndex();

                }
            }

            _atomshells.push_back(_singleatom);
            _startIdx.push_back(_Idx);
            _blocksize.push_back(_size);


            // setup a list of min decay constants per atom
            // for every shell
            _atomindex = 0;
            double _decaymin = 1e7;
            vector< double > _minimal_decay;
            vector < vec > _positions;
            vec _localpos = (*basis->firstShell())->getPos();
            for (vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++) {

                if ((*_row)->getIndex() == _atomindex) {

                    // check all decay constants in this shell
                    for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                        AOGaussianPrimitive* gaussian = *itg;
                        double _decay = gaussian->decay;
                        if (_decay < _decaymin) {
                            _decaymin = _decay;
                        } // decay min check

                    } // Gaussian Primitives 

                } else { // if shell belongs to the actual atom
                    // add to mininal_decay vector
                    _minimal_decay.push_back(_decaymin);
                    _positions.push_back(_localpos);
                    // reset counters
                    _decaymin = 1e7;
                    _localpos = (*_row)->getPos();

                    _atomindex++;

                    // check all decay constants in this shell
                    for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                        AOGaussianPrimitive* gaussian = *itg;
                        double _decay = gaussian->decay;
                        if (_decay < _decaymin) {
                            _decaymin = _decay;
                        } // decay min check

                    } // Gaussian Primitives                                       
                }
            } // all shells

            // push final atom
            _minimal_decay.push_back(_decaymin);
            _positions.push_back(_localpos);


            // for each gridpoint, check the value of exp(-a*(r-R)^2) < 1e-10
            //                             = alpha*(r-R)^2 >~ 20.7

            // each atomic grid
            for (unsigned i = 0; i < _grid.size(); i++) {

                vector< vector<int> > _significant_atoms_atomgrid;

                // each point of the atomic grid
                for (unsigned j = 0; j < _grid[i].size(); j++) {

                    vector<int> _significant_atoms_gridpoint;
                    const vec& grid = _grid[i][j].grid_pos;


                    // check all atoms
                    for (unsigned iatom = 0; iatom < _minimal_decay.size(); iatom++) {

                        vec dist = grid - _positions[iatom];
                        double distsq = dist*dist;

                        // if contribution is smaller than -ln(1e-10), add atom to list
                        if ((_minimal_decay[iatom] * distsq) < 20.7) {
                            _significant_atoms_gridpoint.push_back(iatom);
                        }

                    } // check all atoms

                    _significant_atoms_atomgrid.push_back(_significant_atoms_gridpoint);

                } // all points of this atom grid

                _significant_atoms.push_back(_significant_atoms_atomgrid);

            } // atomic grids



            int total_grid = 0;
            int significant_grid = 0;
            for (unsigned i = 0; i < _significant_atoms.size(); i++) {

                total_grid += _grid[i].size();

                for (unsigned j = 0; j < _significant_atoms[i].size(); j++) {

                    int gridpointsize = _significant_atoms[i][j].size();
                    significant_grid += gridpointsize * (gridpointsize + 1);

                }
            }
            int natoms = _grid.size();

            total_grid = total_grid * (natoms * (natoms + 1)) / 2;

            for (unsigned rowatom = 0; rowatom < _grid.size(); rowatom++) {
                std::vector< ub::matrix<double> > rowmatrix;
                for (unsigned colatom = 0; colatom <= rowatom; colatom++) {
                    rowmatrix.push_back(ub::zero_matrix<double>(_blocksize[colatom], _blocksize[rowatom]));
                }
                dmat_vector.push_back(rowmatrix);
            }



            // parallelization: distribute over threads inside one atom
            int nthreads = 1;
#ifdef _OPENMP
            nthreads = omp_get_max_threads();
#endif


            for (int i = 0; i < nthreads; i++) {

                std::vector< std::vector< ub::matrix<double> > > matrix;
                for (unsigned rowatom = 0; rowatom < _grid.size(); rowatom++) {
                    std::vector< ub::matrix<double> > rowmatrix;
                    for (unsigned colatom = 0; colatom < _grid.size(); colatom++) {
                        rowmatrix.push_back(ub::zero_matrix<double>(_blocksize[rowatom], _blocksize[colatom]));
                    }
                    matrix.push_back(rowmatrix);
                }
                xcmat_vector_thread.push_back(matrix);
            }

            for (unsigned rowatom = 0; rowatom < _grid.size(); rowatom++) {
                std::vector< ub::matrix<double> > rowmatrix;
                for (unsigned colatom = 0; colatom < _grid.size(); colatom++) {
                    rowmatrix.push_back(ub::zero_matrix<double>(_blocksize[colatom], _blocksize[rowatom]));
                }
                xcmat_vector.push_back(rowmatrix);
            }
            return;
        }

        double NumericalIntegration::IntegratePotential_w_PBC(vec rvector, vec boxLen) {

            double result = 0.0;
            double cutoff = min(min(boxLen(0), boxLen(1)), boxLen(2)) / 2.0;
            double vol = boxLen(0) * boxLen(1) * boxLen(2);

            if (density_set) {
                //                #pragma omp parallel for reduction(+:result)
                for (unsigned i = 0; i < _grid.size(); i++) {
                    for (unsigned j = 0; j < _grid[i].size(); j++) {

                        //charge at this point
                        double q = -_grid[i][j].grid_weight * _grid[i][j].grid_density; //density is neg of charge

                        //r-space sum
                        vec dif = _grid[i][j].grid_pos - rvector;

                        for (int k = 0; k < 3; k++) {
                            if (std::abs(dif(k) > boxLen(k)*0.5)) //correct for for bond crossing PBC, if it exists
                                if (dif(k) > 0) //i.x>j.x
                                    dif(k) -= boxLen(k);
                                else //i.x<j.x
                                    dif(k) += boxLen(k);
                        }
                        double dist = abs(dif); //in bohr
                        double potR = q * ((erfc(alpha * dist) / dist) - tools::conv::Pi / (vol * vol * alpha * alpha)); //erfc and shift average pot to 0

                        if (dist < 1.0e-12) { //point is in the same spot as we are evaluating potential
                            result -= 2.0 * (alpha / sqrt(tools::conv::Pi)) * q; //self correction
                        } else if (dist <= cutoff) {
                            //double potR=q/dist;
                            result += potR;
                        }

                    }//j
                }//i
                E_rspace = result;

                //cout<< "r-space sum: " <<result<<endl;

                //k-space sum
                vec r(rvector(0), rvector(1), rvector(2));
                int nKpoints = numK[0] * numK[1] * numK[2];
                double ksum = 0.0;
                //exclude k={0,0,0} (index=0)
                //                #pragma omp parallel for reduction(+:ksum)
                for (int index = 1; index < nKpoints; index++) {
                    double* kp = &(Kcoord[index * 3]);
                    vec K(kp[0], kp[1], kp[2]);
                    std::complex<double> Kr(0, K * r); // ik dot r

                    double potK = prefactor[index]*(Rho_k[index] * std::exp(Kr)).real();
                    ksum += potK;

                }//index
                result += ksum;
                E_kspace = ksum;


                //cout<< "k-space sum: " <<ksum<<endl;
                //cout<< "final Ewald sum: "<<result<<endl;

                //dipole moment correction
                //This is only usefull if we have a countably infinite
                //sphere of periodic cells embeded in a dielectric.
                //If material is trully periodic ("metallic boundary
                //conditions"), then this correction is 0



            } else {
                throw std::runtime_error("Density not calculated");
            }

            return result;
        }

        double NumericalIntegration::CalcDipole_w_PBC(vec rvector, vec boxLen) {
            //as a check, compute the dipole moment of the density distribution.
            //center it on the first atom, for example.

            double z[3] = {0};
            vec dip(z);
            if (density_set) {
                //                #pragma omp parallel for reduction(+:result)
                for (unsigned i = 0; i < _grid.size(); i++) {
                    for (unsigned j = 0; j < _grid[i].size(); j++) {
                        double q = -_grid[i][j].grid_weight * _grid[i][j].grid_density; //density is neg of charge
                        vec dif = _grid[i][j].grid_pos - rvector;
                        for (int k = 0; k < 3; k++) {
                            if (std::abs(dif[k]) > boxLen[k]*0.5) //correct for for bond crossing PBC, if it exists
                                if (dif[k] > 0) //i.x>j.x
                                    dif[k] -= boxLen[k];
                                else //i.x<j.x
                                    dif[k] += boxLen[k];
                        }
                        dip += dif * q;
                    }//i
                }//j
            }//density
            else {
                throw std::runtime_error("Density not calculated");
            }
            return (abs(dip)); //in Bohr*elementary charge
        }

        void NumericalIntegration::FreeKspace() {
            delete[] Rho_k;
            delete[] Kcoord;
            delete[] prefactor;
        }

        /**
         *	Calculate and return the Ewald coefficient (alpha) from cutoff distance
         *	and requested tolerance.
         */
        void NumericalIntegration::findAlpha(double Rc, double dtol) {
            double x = 1.0;
            int i = 20;
            while (erfc(x * Rc) / Rc > dtol) {
                x *= 2.0;
                i++;
            }
            double lastx = x * 2.0;
            double newx;
            while (true) {
                double reldif = (erfc(x * Rc) / Rc - dtol) / dtol;
                if (std::abs(reldif) < 1.0e-3 || i == 0)
                    break;

                if (reldif > 0)
                    newx = (x + lastx)*0.5;
                else
                    newx = x * 0.5;
                lastx = x;
                x = newx;
                i--;
            }
            alpha = x;
        }

        /**
         *
         * @param _atoms
         * @param boxLen
         * @param Kspacing in Angstroms
         */
        void NumericalIntegration::PrepKspaceDensity(vec boxLen, double ext_alpha, std::vector< ctp::QMAtom* > & _local_atomlist, bool ECP, int nK = 0) {

            cout << "box is " << boxLen[0] << " " << boxLen[1] << " " << boxLen[2] << endl;

            //add nuclear charges to _grid to avoid having to do periodic calculations on them separately;
            //_Madelung_grid.clear();
            std::vector< GridContainers::integration_grid > _Nuc;
            Elements _elements;
            for (std::vector< ctp::QMAtom* >::iterator it = _local_atomlist.begin(); it != _local_atomlist.end(); ++it) {
                GridContainers::integration_grid el;
                ctp::QMAtom* atom = *it;
                if (ECP) {
                    el.grid_density = -_elements.getNucCrgECP(atom->type);
                } else {
                    el.grid_density = -_elements.getNucCrg(atom->type);
                }
                el.grid_weight = 1.0;
                el.grid_pos = vec(atom->x, atom->y, atom->z) * tools::conv::ang2bohr;
                //cout<< atom->charge << "@ "<<el.grid_x<<'\t'<<el.grid_y<<'\t'<<el.grid_z<<endl;
                _Nuc.push_back(el);
            }
            if (_Nuc.size() > 0) {
                _grid.push_back(_Nuc);
            }

            //
            //            //fill Madelung grid with density grid, all in to _Madelung_grid[0], so that energy calculation foesn't have a quadrupple nested loop
            //            _Madelung_grid=_grid;



            //cout<<"points in _Mad: " <<_Mad.size() << endl;
            //cout<<"_Mads in _Madelung_grid: " <<_Madelung_grid.size() << endl;



            //compute alpha
            double cutoff = min(min(boxLen[0], boxLen[1]), boxLen[2]) / 2.0;
            //            findAlpha(cutoff, 1.0e-10);
            alpha = ext_alpha;
            double fourasq = 4.0 * alpha*alpha;
            cout << "found alpha = " << alpha << "\t rel err of r-sum ~ " << erfc(alpha * cutoff) / cutoff << endl;
            //            cout<<"found alpha = "<< alpha <<"\t";


            //this is going to be slow, as points we have density for are not on a periodic grid,
            //so will use simple Ewald summation, not any of the FFT methods.

            //find number of k-vectors in each lattice direction
            double minSq = 1.0e15;
            for (int i = 0; i < 3; i++) {
                if (nK == 0) {
                    //numK[i]=1+(boxLen[i]/Kspacing);
                    //if(numK[i]%2 == 0) numK[i]++; //keep it odd

                    //find maxK
                    int maxK;
                    double maxKsq;
                    double err;
                    double twoPiL = 2.0 * tools::conv::Pi / boxLen[i];
                    for (maxK = 2; true; maxK++) {
                        maxKsq = maxK * maxK * twoPiL*twoPiL;
                        err = std::exp(-maxKsq / fourasq) / maxKsq;
                        if (err < 1.0e-7) {
                            break;
                        }
                    }

                    maxK = max(16, maxK);

                    numK[i] = 2 * maxK + 1;
                    minSq = min(maxKsq, minSq);
                } else {
                    numK[i] = nK;
                }
            }
            cout << "numK={" << numK[0] << ", " << numK[1] << ", " << numK[2] << "}" << endl;
            cout << "rel err of k-sum ~ " << std::exp(-minSq / fourasq) / minSq << endl;


            //allocate
            int nKpoints = numK[0] * numK[1] * numK[2];
            Rho_k = new std::complex<double>[nKpoints];
            Kcoord = new double[nKpoints * 3];
            prefactor = new double[nKpoints];

            //fill Kcoord and prefactor
            // [0, 1, 2, 3, ..., N/2 -1, N/2, -N/2, -N/2 +1, -N/2 +1, ..., -3, -2, -1] (note no 0)
            // k={0,0,0} is the very first entry (index=0)


            double invvolume = 1.0 / (boxLen[0] * boxLen[1] * boxLen[2]);
            double pre = invvolume * 4.0 * tools::conv::Pi;

            //cout<<"pre=4*Pi/Volume="<<pre<<endl;
            //cout<<"fourasq=4*alpha^2="<<fourasq<<endl;

#pragma omp parallel for
            for (int l = 0; l < numK[0]; l++) {
                int L, M, N;
                double ksq;
                if (l > numK[0] / 2) {
                    L = l - numK[0];
                } else {
                    L = l;
                }
                //cout<<L<<endl;
                for (int m = 0; m < numK[1]; m++) {
                    if (m > numK[1] / 2) {
                        M = m - numK[1];
                    } else {
                        M = m;
                    }
                    for (int n = 0; n < numK[2]; n++) {
                        if (n > numK[2] / 2) {
                            N = n - numK[2];
                        } else {
                            N = n;
                        }

                        int index = l * (numK[1] * numK[2]) + m * numK[2] + n;

                        double* kp = &(Kcoord[index * 3]);
                        kp[0] = 2.0 * tools::conv::Pi * L / boxLen[0];
                        kp[1] = 2.0 * tools::conv::Pi * M / boxLen[1];
                        kp[2] = 2.0 * tools::conv::Pi * N / boxLen[2];

                        ksq = (kp[0] * kp[0])+(kp[1] * kp[1])+(kp[2] * kp[2]);
                        //cout<<"{L,M,N}={"<<L<<", "<<M<<", "<<N<<"}\n";
                        //                        cout<<"{l,m,n}={"<<l<<", "<<m<<", "<<n<<"}\t";
                        //                        cout<<"k={"<<kp[0]<<", "<<kp[1]<<", "<<kp[2]<<"}\tksq="<<ksq<<"\t";
                        //                        cout<<"prefactor="<<(pre/ksq) * std::exp(-ksq/fourasq)<<endl;
                        prefactor[index] = (pre / ksq) * std::exp(-ksq / fourasq);

                        //index++;
                    }
                }
            }


            //fill Rho_k
#pragma omp parallel for
            for (int index = 0; index < nKpoints; index++) {
                double* kp = &(Kcoord[index * 3]);
                vec K(kp[0], kp[1], kp[2]);
                Rho_k[index] = 0.0;
                for (unsigned i = 0; i < _grid.size(); i++) {
                    for (unsigned j = 0; j < _grid[i].size(); j++) {
                        vec r = _grid[i][j].grid_pos;
                        std::complex<double> nKr(0, -K * r); // -ik dot r
                        Rho_k[index] -= _grid[i][j].grid_weight * _grid[i][j].grid_density * std::exp(nKr); //density is neg of charge
                    }//j
                }//i
            }//index

        }

        void NumericalIntegration::PrepKspaceDensity_gromacs_like(vec boxLen, double ext_alpha, std::vector< ctp::QMAtom* > & _local_atomlist, bool ECP, Grid &eval_grid, int nK) {
            cout << "box is " << boxLen[0] << " " << boxLen[1] << " " << boxLen[2] << endl; //already in Bohr

            alpha = ext_alpha;
            double fourasq = 4.0 * alpha*alpha;

            //add nuclear charges to _grid to avoid having to do periodic calculations on them separately;
            std::vector< GridContainers::integration_grid > _Nuc;
            Elements _elements;
            for (std::vector< ctp::QMAtom* >::iterator it = _local_atomlist.begin(); it != _local_atomlist.end(); ++it) {
                GridContainers::integration_grid el;
                ctp::QMAtom* atom = *it;
                if (ECP) {
                    el.grid_density = -_elements.getNucCrgECP(atom->type);
                } else {
                    el.grid_density = -_elements.getNucCrg(atom->type);
                }
                el.grid_weight = 1.0;
                el.grid_pos = vec(atom->x, atom->y, atom->z) * tools::conv::ang2bohr;
                //cout<< atom->charge << "@ "<<el.grid_x<<'\t'<<el.grid_y<<'\t'<<el.grid_z<<endl;
                _Nuc.push_back(el);
            }
            if (_Nuc.size() > 0) {
                _grid.push_back(_Nuc);
            }


            //this is going to be slow, as points we have density for are not on a periodic grid,
            //so will use simple Ewald summation, not any of the FFT methods.

            //find number of k-vectors in each lattice direction
            double minSq = 1.0e15;
            for (int i = 0; i < 3; i++) {
                /*
                //numK[i]=1+(boxLen[i]/Kspacing);
                //if(numK[i]%2 == 0) numK[i]++; //keep it odd

                //find maxK
                int maxK;
                double maxKsq;
                double err;
                double twoPiL=2.0*tools::conv::Pi/boxLen[i];
                for(maxK=2; true; maxK++){
                    maxKsq=maxK*maxK*twoPiL*twoPiL;
                    err=std::exp(-maxKsq/fourasq)/maxKsq;
                    if(err<1.0e-7){
                            break;
                    }
                }

                maxK=max(16, maxK);
                
                
                numK[i]=maxK;
                minSq=min(maxKsq, minSq);
                 */
                numK[i] = nK;
            }
            cout << "numK={" << numK[0] << ", " << numK[1] << ", " << numK[2] << "}" << endl;
            cout << "rel err of k-sum ~ " << std::exp(-minSq / fourasq) / minSq << endl;


            //allocate space for eikr
            eikr.resize(3); //3 dim
            for (unsigned m = 0; m < 3; m++) {
                eikr[m].resize(numK[m]);
                for (unsigned k = 0; k < numK[m]; k++) {
                    eikr[m][k].resize(_grid.size());
                    for (unsigned i = 0; i < _grid.size(); i++) {
                        eikr[m][k][i].resize(_grid[i].size());
                    }
                }
                lll[m] = 2.0 * boost::math::constants::pi<double>() / boxLen[m]; //units are 1/Bohr
            }
            //eikr indexing is [DIM][k][i][j]

            //based on gromacs 4.6 tabulate_eir()
            for (unsigned i = 0; i < _grid.size(); i++) {
#pragma omp parallel for
                for (unsigned j = 0; j < _grid[i].size(); j++) {
                    vec r = _grid[i][j].grid_pos; //grid is in Bohr
                    for (unsigned m = 0; (m < 3); m++) {
                        eikr[m][0][i][j] = std::complex<double>(1, 0);
                        eikr[m][1][i][j] = std::complex<double>(cos(r[m] * lll[m]), sin(r[m] * lll[m]));
                        for (unsigned k = 2; k < numK[m]; k++) {
                            eikr[m][k][i][j] = eikr[m][k - 1][i][j] * eikr[m][1][i][j];
                        }//k
                    }//m
                }//j
            }//i


            //now do the same for the eval_grid (eikR), the points were potential is to be evaluated
            //allocate space for eikR
            eikR.resize(3); //3 dim
            double lll[3];
            for (unsigned m = 0; m < 3; m++) {
                eikR[m].resize(numK[m]);
                for (unsigned k = 0; k < numK[m]; k++) {
                    eikR[m][k].resize(_grid.size());
                }
            }
            //eikR indexing is [DIM][k][i]

            //based on gromacs 4.6 tabulate_eir()
#pragma omp parallel for
            for (unsigned i = 0; i < _grid.size(); i++) {
                vec R = eval_grid.getGrid()[i] * tools::conv::nm2bohr; //this is in Bohr, eval_grid is in nm
                for (unsigned m = 0; (m < 3); m++) {
                    eikR[m][0][i] = std::complex<double>(1, 0);
                    eikR[m][1][i] = std::complex<double>(cos(R(m) * lll[m]), sin(R(m) * lll[m]));
                    for (unsigned k = 2; k < numK[m]; k++) {
                        eikR[m][k][i] = eikR[m][k - 1][i] * eikR[m][1][i];
                    }//k
                }//m
            }//i

        }

        void NumericalIntegration::FillMadelungGrid(vec boxLen, int natomsonside) {
            //fill _Madelung_grid;
            _Madelung_grid.clear();
            std::vector< GridContainers::integration_grid > _Mad;
            double a = boxLen[0]; //in bohr
            for (int l = 0; l < natomsonside; l++) {
                for (int m = 0; m < natomsonside; m++) {
                    for (int n = 0; n < natomsonside; n++) {
                        GridContainers::integration_grid el;
                        el.grid_density = std::pow(-1.0, l + m + n);
                        el.grid_weight = 1.0;
                        el.grid_pos = vec(l, m, n) * (a / natomsonside);
                        _Mad.push_back(el);
                        //cout<< "q= "<< el.grid_density << "\t @ " << el.grid_x << " " << el.grid_y << " " << el.grid_z << "\t box size:"<< boxLen[0] << endl;
                    }
                }
            }
            _Madelung_grid.push_back(_Mad);
            _grid = _Madelung_grid;
            //exit(0);
        }

        void NumericalIntegration::IntegratePotential_w_PBC_gromacs_like(Grid &eval_grid, vec boxLen, ub::vector<double>& _ESPatGrid) {

            double cutoff = min(min(boxLen[0], boxLen[1]), boxLen[2]) / 2.0; //Bohr
            double vol = boxLen[0] * boxLen[1] * boxLen[2]; //Bohr^3
            //eval_grid is in nm

            if (density_set) {
#pragma omp parallel for
                for (unsigned p = 0; p < _ESPatGrid.size(); p++) {
                    vec rvector = eval_grid.getGrid()[p] * tools::conv::nm2bohr; //Bohr
                    for (unsigned i = 0; i < _grid.size(); i++) {
                        for (unsigned j = 0; j < _grid[i].size(); j++) {

                            //charge at this point
                            double q = -_grid[i][j].grid_weight * _grid[i][j].grid_density; //density is neg of charge
                            //double q = -_Madelung_grid[i][j].grid_weight * _Madelung_grid[i][j].grid_density; //density is neg of charge

                            //r-space sum
                            vec dif = _grid[i][j].grid_pos - rvector; //Bohr

                            for (unsigned k = 0; k < 3; k++) {
                                if (std::abs(dif[k]) > boxLen[k]*0.5) //correct for for bond crossing PBC, if it exists
                                    if (dif[k] > 0) //i.x>j.x
                                        dif[k] -= boxLen[k];
                                    else //i.x<j.x
                                        dif[k] += boxLen[k];
                            }
                            double dist = abs(dif); //in bohr
                            double potR = q * ((erfc(alpha * dist) / dist) - tools::conv::Pi / (vol * vol * alpha * alpha)); //erfc and shift average pot to 0

                            if (dist < 1.0e-12) { //point is in the same spot as we are evaluating potential
                                _ESPatGrid[p] -= 2.0 * (alpha / sqrt(tools::conv::Pi)) * q; //self correction
                            } else if (dist <= cutoff) {
                                //double potR=q/dist;
                                _ESPatGrid[p] += potR;
                            }
                        }//j
                    }//i
                }//p
                E_rspace = _ESPatGrid[0];
                E_erfc = erfc(alpha);


                /*
                                //k-space sum
                                //adapted from gromacs 4.6 do_ewald_pot()
                                int      lowiy, lowiz;
                                double   tmp, cs, ss;

                                lowiy        = 0;
                                lowiz        = 1;
                                //for (ix = 0; ix < numK[0]; ix++)
                                unsigned ix = 0; //do ix=0 first, then do everything else in parallel
                                {
                                    double mx = ix*lll[0];
                                    for (unsigned iy = lowiy; iy < numK[1]; iy++)
                                    {
                                        double my = iy*lll[1];
                                        std::vector <std::vector <std::complex<double>> > tab_xy(_grid.size());
                                        std::vector <std::vector <std::complex<double>> > tab_qxyz(_grid.size());
                                        if (iy >= 0)
                                        {
                                            for (unsigned i = 0; i < _grid.size(); i++)
                                            {
                                                tab_xy[i].resize(_grid[i].size());
                                                tab_qxyz[i].resize(_grid[i].size());
                                                for(unsigned j = 0; j< _grid[i].size(); j++){
                                                    tab_xy[i][j] = eikr[0][ix][i][j]*eikr[1][iy][i][j];
                                                }//j
                                            }//i
                                        }
                                        else
                                        {
                                            for (unsigned i = 0; i < _grid.size(); i++)
                                            {
                                                tab_xy[i].resize(_grid[i].size());
                                                tab_qxyz[i].resize(_grid[i].size());
                                                for(unsigned j = 0; j< _grid[i].size(); j++){
                                                    tab_xy[i][j] = eikr[0][ix][i][j]*std::conj(eikr[1][iy][i][j]);
                                                }//j
                                            }//i
                                        }

                                        for (unsigned iz = lowiz; iz < numK[2]; iz++)
                                        {
                                            double mz  = iz*lll[2];
                                            double m2  = mx*mx+my*my+mz*mz;
                                            double ak  = exp(m2*factor)/m2;
                                            if (iz >= 0)
                                            {
                                                for (unsigned i = 0; i < _grid.size(); i++)
                                                    for(unsigned j = 0; j< _grid[i].size(); j++)
                                                        tab_qxyz[i][j] = charge[i][j]*tab_xy[i][j]*eikr[2][iz][i][j];
                                            }
                                        }//iz

                                    }//iy
                                }//ix=0
                 */



                //adapted from gromacs 4.6 do_ewald_pot()

                int ix, iy, iz;
                double tmp, cs, ss, ak, mx, my, mz, m2;
                double factor = -1.0 / (4 * alpha * alpha);
                double scaleRecip = 1.0 / (vol); //final units of potential are Hartree/e (atomic units)
                std::vector<std::vector<std::complex<double> > > tab_xy;
                std::vector<std::vector<std::complex<double> > > tab_qxyz; //charge distrib.
                std::vector<std::complex<double> > tab_R_xyz; //where to evaluate
                std::vector<std::complex<double> > tab_R_xy;
                //allocate space
                tab_xy.resize(_grid.size());
                tab_qxyz.resize(_grid.size());
                tab_R_xy.resize(eval_grid.getsize());
                tab_R_xyz.resize(eval_grid.getsize());
                for (unsigned i = 0; i < _grid.size(); i++) {
                    tab_xy[i].resize(_grid[i].size());
                    tab_qxyz[i].resize(_grid[i].size());
                }


                unsigned int lowiy = 0;
                unsigned int lowiz = 1;
                for (ix = 0; ix < numK[0]; ix++) {
                    mx = ix * lll[0];
                    for (iy = lowiy; iy < numK[1]; iy++) {
                        my = iy * lll[1];
                        if (iy >= 0) {
                            for (unsigned i = 0; i < _grid.size(); i++) //charge distrib.
                            {
                                for (unsigned j = 0; j < _grid[i].size(); j++) {
                                    tab_xy[i][j] = eikr[0][ix][i][j] * eikr[1][iy][i][j]; //eikr indexing is [DIM][k][i][j]
                                }//j
                            }//i
                            for (unsigned n = 0; n < eval_grid.getsize(); n++) //where to evaluate
                            {
                                tab_R_xy[n] = eikR[0][ix][n] * eikR[1][iy][n];
                            }//n
                        } else //negative k component
                        {
                            for (unsigned i = 0; i < _grid.size(); i++) //charge distrib.
                            {
                                for (unsigned j = 0; j < _grid[i].size(); j++) {
                                    tab_xy[i][j] = eikr[0][ix][i][j] * std::conj(eikr[1][-iy][i][j]); //eikr indexing is [DIM][k][i][j]
                                }//j
                            }//i
                            for (unsigned n = 0; n < eval_grid.getsize(); n++) //where to evaluate
                            {
                                tab_R_xyz[n] = eikR[0][ix][n] * std::conj(eikR[1][-iy][n]);
                            }//n
                        }
                        for (iz = lowiz; iz < numK[2]; iz++) {
                            mz = iz * lll[2];
                            m2 = mx * mx + my * my + mz*mz;
                            ak = exp(m2 * factor) / m2;
                            if (iz >= 0) {
                                for (unsigned i = 0; i < _grid.size(); i++) //charge distrib.
                                {
                                    for (unsigned j = 0; j < _grid[i].size(); j++) {
                                        tab_qxyz[i][j] = -_grid[i][j].grid_weight * _grid[i][j].grid_density * tab_xy[i][j] * eikr[2][iz][i][j]; //eikr indexing is [DIM][k][i][j]
                                    }//j
                                }//i
                                for (unsigned n = 0; n < eval_grid.getsize(); n++) //where to evaluate
                                {
                                    tab_R_xyz[n] = tab_R_xy[n] * eikR[2][iz][n];
                                }//n
                            } else {
                                for (unsigned i = 0; i < _grid.size(); i++) //charge distrib.
                                {
                                    for (unsigned j = 0; j < _grid[i].size(); j++) {
                                        tab_qxyz[i][j] = -_grid[i][j].grid_weight * _grid[i][j].grid_density * tab_xy[i][j] * std::conj(eikr[2][-iz][i][j]); //eikr indexing is [DIM][k][i][j]
                                    }//j
                                }//i
                                for (unsigned n = 0; n < eval_grid.getsize(); n++) //where to evaluate
                                {
                                    tab_R_xyz[n] = tab_R_xy[n] * std::conj(eikR[2][-iz][n]);
                                }//n
                            }

                            cs = ss = 0;
                            for (unsigned i = 0; i < _grid.size(); i++) {
                                for (unsigned j = 0; j < _grid[i].size(); j++) {
                                    cs += tab_qxyz[i][j].real();
                                    ss += tab_qxyz[i][j].imag();
                                }
                            }

                            for (int n = 0; n < eval_grid.getsize(); n++) {
                                tmp = ak * (cs * tab_R_xyz[n].real() + ss * tab_R_xyz[n].imag());
                                _ESPatGrid(n) += tmp * 2 * scaleRecip;
                            }

                            lowiz = 1 - numK[2];
                        }
                        lowiy = 1 - numK[1];
                    }
                }



            } else {
                throw std::runtime_error("Density not calculated");
            }
            E_kspace = _ESPatGrid[0] - E_rspace;
        }

        double NumericalIntegration::IntegrateDensity_Atomblock(const ub::matrix<double>& _density_matrix, AOBasis* basis) {
            if (_significant_atoms.size() < 1) {
                throw runtime_error("NumericalIntegration::IntegrateDensity_Atomblock:significant atoms not found yet.");
            }
            double result = 0.0;

            // parallelization: distribute over threads inside one atom
            int nthreads = 1;
#ifdef _OPENMP
            nthreads = omp_get_max_threads();
#endif

            std::vector<double> Density_thread;
            for (int i_thread = 0; i_thread < nthreads; i_thread++) {
                Density_thread.push_back(0.0);
            }

            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
                // for each point in atom grid

                // number of points in this atomgrid
                int atom_points = _grid[i].size();
                // divide among threads
                int atom_points_per_thread = atom_points / nthreads;
                std::vector<int> _thread_start;
                std::vector<int> _thread_stop;
                for (int i_thread = 0; i_thread < nthreads; i_thread++) {
                    _thread_start.push_back(i_thread * atom_points_per_thread);
                    _thread_stop.push_back((i_thread + 1) * atom_points_per_thread);
                }
                // final stop must be size
                _thread_stop[nthreads - 1] = atom_points;



#pragma omp parallel for
                for (int i_thread = 0; i_thread < nthreads; i_thread++) {
                    for (int j = _thread_start[i_thread]; j < _thread_stop[i_thread]; j++) {
                        //boost::timer::cpu_times t0 = cpu_t.elapsed();

                        // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )
                        //ub::matrix<double> AOgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 1);

                        ub::matrix<double> AOgrid = ub::zero_matrix<double>(1, basis->_AOBasisSize); // TRY MORE USEFUL DATA

                        ub::matrix<double> rho_mat = ub::zero_matrix<double>(1, 1);


                        // evaluate AO Functions for all shells, NOW BLOCKWISE

                        // for each significant atom for this grid point
                        for (unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size(); sigrow++) {

                            // this atom
                            int rowatom = _significant_atoms[i][j][sigrow];


                            // for each shell in this atom
                            for (unsigned ishell = 0; ishell < _atomshells[rowatom].size(); ishell++) {
                                //boost::timer::cpu_times tstartshells = cpu_t.elapsed();
                                AOShellIterator _row = _atomshells[rowatom][ishell];
                                // for density, fill sub-part of AOatgrid

                                ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, 0, 1, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());

                                (*_row)->EvalAOspace(_AOgridsub, _grid[i][j].grid_pos);


                            } // shell in atom
                        }

                        for (unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size(); sigrow++) {

                            // this atom
                            int rowatom = _significant_atoms[i][j][sigrow];
                            ub::matrix<double> _temp = ub::zero_matrix<double>(1, _blocksize[rowatom]);

                            ub::matrix_range< ub::matrix<double> > _AOgridrow = ub::subrange(AOgrid, 0, 1, _startIdx[rowatom], _startIdx[rowatom] + _blocksize[rowatom]);

                            // for each atom

                            for (unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size(); sigcol++) {
                                int colatom = _significant_atoms[i][j][sigcol];


                                // get the already calculated AO values

                                ub::matrix_range< ub::matrix<double> > _AOgridcol = ub::subrange(AOgrid, 0, 1, _startIdx[colatom], _startIdx[colatom] + _blocksize[colatom]);

                                ub::matrix_range<const ub::matrix<double> > DMAT_here = ub::subrange(_density_matrix, _startIdx[colatom], _startIdx[colatom] + _blocksize[colatom], _startIdx[rowatom], _startIdx[rowatom] + _blocksize[rowatom]);

                                _temp += ub::prod(_AOgridcol, DMAT_here);


                            } //col shells

                            rho_mat += ub::prod(_temp, ub::trans(_AOgridrow));

                        } // row shells 


                        _grid[i][j].grid_density = rho_mat(0, 0);
                        Density_thread[i_thread] += _grid[i][j].grid_weight * _grid[i][j].grid_density;


                    } // j: for each point in atom grid
                }// each thread
            } // i: for each atom grid

            for (int i_thread = 0; i_thread < nthreads; i_thread++) {

                result += Density_thread[i_thread];
            }
            density_set = true;
            return result;
        }

        
        
        
         /**
         * Computes electron density belonging to atoms listed in AtomIndeces.
         * Based on IntegrateDensity_Atomblock()
         *
         * @param _density_matrix
         * @param basis
         * @param AtomIndeces of the atoms in the molecule
         * @return electron density
         */
        double NumericalIntegration::IntegrateDensity_Molecule(ub::matrix<double>& _density_matrix, AOBasis* basis, std::vector<int> AtomIndeces){

            double result=0.0;

            // generate a list of shells for each atom
            typedef vector< AOShell* >::iterator AOShellIterator;
            vector< vector< AOShellIterator > > _atomshells;
            vector< AOShellIterator > _singleatom;

            vector < int > _startIdx;
            vector < int > _blocksize;

            int _atomindex = 0;
            int _Idx       = 0;
            int _size      = 0;

            for (vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++) {


                if ( (*_row)->getIndex() == _atomindex ){

                    _singleatom.push_back(_row);
                    _size += (*_row)->getNumFunc();


                } else {

                    // append _singleatom to _atomshells
                    _atomshells.push_back(_singleatom);
                    _startIdx.push_back( _Idx );
                    _blocksize.push_back(_size);
                    // reset _singleatom
                    _singleatom.clear();
                    _size = (*_row)->getNumFunc();
                    _Idx       = (*_row)->getStartIndex();
                    _singleatom.push_back(_row);
                    _atomindex = (*_row)->getIndex();

                }

            }

            _atomshells.push_back(_singleatom);
            _startIdx.push_back( _Idx );
                    _blocksize.push_back(_size);



            // setup a list of min decay constants per atom
            // for every shell
            _atomindex = 0;
            double _decaymin = 1e7;
            vector< double > _minimal_decay;
            vector < vec > _positions;
            vec _localpos = (*basis->firstShell())->getPos();
            for ( vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++   ) {

                 if ( (*_row)->getIndex() == _atomindex ){

                     // check all decay constants in this shell
                     for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                         AOGaussianPrimitive* gaussian = *itg;
                         double _decay = gaussian->decay;
                         if (_decay < _decaymin) {
                             _decaymin = _decay;
                         } // decay min check

                     } // Gaussian Primitives

                 } else {  // if shell belongs to the actual atom
                     // add to mininal_decay vector
                     _minimal_decay.push_back(_decaymin);
                     _positions.push_back( _localpos );
                     // reset counters
                     _decaymin = 1e7;
                     _localpos = (*_row)->getPos();

                     _atomindex++;

                     // check all decay constants in this shell
                     for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                         AOGaussianPrimitive* gaussian = *itg;
                         double _decay = gaussian->decay;
                         if (_decay < _decaymin) {
                             _decaymin = _decay;
                         } // decay min check

                     } // Gaussian Primitives
                 }
            } // all shells

            // push final atom
            _minimal_decay.push_back(_decaymin);
            _positions.push_back( _localpos );



            // for each gridpoint, check the value of exp(-a*(r-R)^2) < 1e-10
            //                             = alpha*(r-R)^2 >~ 20.7

            vector< vector< vector<int> > > _significant_atoms;

            // each atomic grid
            for (unsigned i = 0; i < _grid.size(); i++) {

                vector< vector<int> > _significant_atoms_atomgrid;

                // each point of the atomic grid
                for (unsigned j = 0; j < _grid[i].size(); j++) {

                    vector<int> _significant_atoms_gridpoint;
                    vec grid=_grid[i][j].grid_pos;

                    // check all atoms
                    for ( unsigned iatom = 0 ; iatom < _minimal_decay.size(); iatom++){

                        vec dist = grid - _positions[iatom];
                        double distsq = dist.getX()*dist.getX() + dist.getY()*dist.getY()  + dist.getZ()*dist.getZ() ;

                        // if contribution is smaller than -ln(1e-10), add atom to list
                        if ( (_minimal_decay[iatom] * distsq) < 20.7 ){
                            _significant_atoms_gridpoint.push_back(iatom);
                        }

                    } // check all atoms

                    _significant_atoms_atomgrid.push_back(  _significant_atoms_gridpoint );

                } // all points of this atom grid

                _significant_atoms.push_back(_significant_atoms_atomgrid);

            } // atomic grids





            // parallelization: distribute over threads inside one atom
            int nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif

            std::vector<double> Density_thread;
            for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                Density_thread.push_back(0.0);
            }

            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
	      // for each point in atom grid

                // number of points in this atomgrid
                int atom_points = _grid[i].size();
                // divide among threads
                int atom_points_per_thread = atom_points/nthreads;
                std::vector<int> _thread_start;
                std::vector<int> _thread_stop;
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                    _thread_start.push_back( i_thread * atom_points_per_thread );
                    _thread_stop.push_back( (i_thread + 1) * atom_points_per_thread );
                }
                // final stop must be size
                _thread_stop[nthreads-1] = atom_points;


                #pragma omp parallel for
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                for (int j = _thread_start[i_thread]; j < _thread_stop[i_thread]; j++) {

                   // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )
                   ub::matrix<double> AOgrid = ub::zero_matrix<double>(1, basis->_AOBasisSize); //basis function values at grid for all significant atoms
                   ub::matrix<double> rho_mat = ub::zero_matrix<double>(1,1);

		    // evaluate AO Functions for all shells, NOW BLOCKWISE

                    // for each significant atom for this grid point
                    for ( unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size() ; sigrow++){

                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];

                        // for each shell in this atom
                        for ( unsigned ishell = 0 ; ishell < _atomshells[rowatom].size() ; ishell++ ){

                            AOShellIterator _row = _atomshells[rowatom][ishell];
                            AOShell* _shell = *_row;
                            // for density, fill sub-part of AOatgrid
                            ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, 0, 1, _shell->getStartIndex(), _shell->getStartIndex()+_shell->getNumFunc());
                            _shell->EvalAOspace(_AOgridsub, _grid[i][j].grid_pos);

                        }  // shell in atom
                    }

                   for ( unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size() ; sigrow++){

                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                        //only do things if this atom belongs to the molecule of interest and is significant
                        if(std::find(AtomIndeces.begin(), AtomIndeces.end(), rowatom) != AtomIndeces.end()){

                            ub::matrix<double> _temp     = ub::zero_matrix<double>(1,_blocksize[rowatom]);

                            ub::matrix_range< ub::matrix<double> > _AOgridrow     = ub::subrange(    AOgrid, 0,1, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);

                            // for each atom

                            for ( unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size() ; sigcol++){
                                int colatom = _significant_atoms[i][j][sigcol];
                                ub::matrix_range< ub::matrix<double> >     _AOgridcol = ub::subrange(    AOgrid, 0, 1, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);
                                ub::matrix_range< ub::matrix<double> > DMAT_here = ub::subrange( _density_matrix, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom], _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);

                                _temp     += ub::prod( _AOgridcol, DMAT_here);


                            } //col shells

                            rho_mat  += ub::prod(_temp, ub::trans( _AOgridrow) );
                        } // atom belongs to molecule

                    } // row shells


                    _grid[i][j].grid_density  =rho_mat(0,0);
                    Density_thread[i_thread] += _grid[i][j].grid_weight * _grid[i][j].grid_density;
                } // j: for each point in atom grid
                }// each thread
            } // i: for each atom grid

             for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                 //cout << result << endl;
                result += Density_thread[i_thread];
                }
            density_set=true;

            return(result);
        }
        
        
        
        
        
        double NumericalIntegration::StupidIntegrate(std::vector<double>& _data) {


            double integral = 0.0;
            int _i_point = 0;
            for (unsigned i = 0; i < _grid.size(); i++) {
                for (unsigned j = 0; j < _grid[i].size(); j++) {

                    integral += _data[_i_point] * _grid[i][j].grid_weight;

                    _i_point++;

                }
            }

            return integral;

        }

        std::vector<vec const *> NumericalIntegration::getGridpoints() {

            std::vector<vec const *> gridpoints;


            for (unsigned i = 0; i < _grid.size(); i++) {
                for (unsigned j = 0; j < _grid[i].size(); j++) {
                    gridpoints.push_back(&_grid[i][j].grid_pos);

                }
            }
            return gridpoints;
        }

        void NumericalIntegration::GridSetup(string type, BasisSet* bs, vector<ctp::QMAtom*> _atoms, AOBasis* basis) {

            const double pi = boost::math::constants::pi<double>();
            // get GridContainer
            GridContainers _grids;

            // get radial grid per element
            EulerMaclaurinGrid _radialgrid;
            _radialgrid.getRadialGrid(bs, _atoms, type, _grids); // this checks out 1:1 with NWChem results! AWESOME


            map<string, GridContainers::radial_grid>::iterator it;

            LebedevGrid _sphericalgrid;

            for (it = _grids._radial_grids.begin(); it != _grids._radial_grids.end(); ++it) {
                _sphericalgrid.getSphericalGrid(_atoms, type, _grids);

            }


            // for the partitioning, we need all inter-center distances later, stored in one-directional list
            int ij = 0;
            Rij.push_back(0.0); // 1st center "self-distance"

            vector< ctp::QMAtom* > ::iterator ait;
            vector< ctp::QMAtom* > ::iterator bit;
            int i = 1;
            for (ait = _atoms.begin() + 1; ait != _atoms.end(); ++ait) {
                // get center coordinates in Bohr
                double x_a = (*ait)->x * tools::conv::ang2bohr;
                double y_a = (*ait)->y * tools::conv::ang2bohr;
                double z_a = (*ait)->z * tools::conv::ang2bohr;
                int j = 0;
                for (bit = _atoms.begin(); bit != ait; ++bit) {
                    ij++;
                    // get center coordinates in Bohr
                    double x_b = (*bit)->x * tools::conv::ang2bohr;
                    double y_b = (*bit)->y * tools::conv::ang2bohr;
                    double z_b = (*bit)->z * tools::conv::ang2bohr;

                    Rij.push_back(1.0 / sqrt((x_a - x_b)*(x_a - x_b) + (y_a - y_b)*(y_a - y_b) + (z_a - z_b)*(z_a - z_b)));




                    j++;
                } // atoms
                Rij.push_back(0.0); // self-distance again
                i++;
            } // atoms



            int i_atom = 0;
            _totalgridsize = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                // get center coordinates in Bohr
                std::vector< GridContainers::integration_grid > _atomgrid;
                const vec atomA_pos = vec((*ait)->x * tools::conv::ang2bohr, (*ait)->y * tools::conv::ang2bohr, (*ait)->z * tools::conv::ang2bohr);

                string name = (*ait)->type;

                // get radial grid information for this atom type
                GridContainers::radial_grid _radial_grid = _grids._radial_grids.at(name);


                // get spherical grid information for this atom type
                GridContainers::spherical_grid _spherical_grid = _grids._spherical_grids.at(name);

                // maximum order (= number of points) in spherical integration grid
                int maxorder = _sphericalgrid.Type2MaxOrder(name, type);
                int maxindex = _sphericalgrid.getIndexFromOrder(maxorder);

                // for pruning of integration grid, get interval boundaries for this element
                std::vector<double> PruningIntervals = _radialgrid.getPruningIntervals(name);
                //  cout << " Pruning Intervals: " << PruningIntervals[0] << " " << PruningIntervals[1] << " " << PruningIntervals[2] << " " << PruningIntervals[3] << endl;

                int current_order = 0;
                // get spherical grid
                std::vector<double> _theta;
                std::vector<double> _phi;
                std::vector<double> _weight;

                // for each radial value
                for (unsigned _i_rad = 0; _i_rad < _radial_grid.radius.size(); _i_rad++) {
                    double r = _radial_grid.radius[_i_rad];
                    int order;
                    // which Lebedev order for this point?
                    if (maxindex == 1) {
                        // smallest possible grid anyway, nothing to do
                        order = maxorder;
                    } else if (maxindex == 2) {
                        // only three intervals
                        if (r < PruningIntervals[0]) {
                            order = _sphericalgrid.getOrderFromIndex(1); //1;
                        } else if ((r >= PruningIntervals[0]) && (r < PruningIntervals[3])) {
                            order = _sphericalgrid.getOrderFromIndex(2);
                        } else {
                            order = _sphericalgrid.getOrderFromIndex(1);
                        } // maxorder == 2
                    } else {
                        // five intervals
                        if (r < PruningIntervals[0]) {
                            order = _sphericalgrid.getOrderFromIndex(int(2));
                        } else if ((r >= PruningIntervals[0]) && (r < PruningIntervals[1])) {
                            order = _sphericalgrid.getOrderFromIndex(4);
                        } else if ((r >= PruningIntervals[1]) && (r < PruningIntervals[2])) {
                            order = _sphericalgrid.getOrderFromIndex(max(maxindex - 1, 4));
                        } else if ((r >= PruningIntervals[2]) && (r < PruningIntervals[3])) {
                            order = maxorder;
                        } else {
                            order = _sphericalgrid.getOrderFromIndex(max(maxindex - 1, 1));
                        }
                    }



                    // get new spherical grid, if order changed
                    if (order != current_order) {
                        _theta.clear();
                        _phi.clear();
                        _weight.clear();

                        _sphericalgrid.getUnitSphereGrid(order, _theta, _phi, _weight);
                        current_order = order;
                    }

                    // for each (theta,phi)
                    // for (int _i_sph = 0; _i_sph < _spherical_grid.phi.size(); _i_sph++) {

                    for (unsigned _i_sph = 0; _i_sph < _phi.size(); _i_sph++) {

                        double p = _phi[_i_sph] * pi / 180.0; // back to rad
                        double t = _theta[_i_sph] * pi / 180.0; // back to rad
                        double ws = _weight[_i_sph];

                        const vec s = vec(sin(p) * cos(t), sin(p) * sin(t), cos(p));

                        GridContainers::integration_grid _gridpoint;
                        _gridpoint.grid_pos = atomA_pos + r*s;

                        _gridpoint.grid_weight = _radial_grid.weight[_i_rad] * ws;

                        _atomgrid.push_back(_gridpoint);


                    } // spherical gridpoints
                } // radial gridpoint


                // get all distances from grid points to centers
                std::vector< std::vector<double> > rq;
                // for each center
                for (bit = _atoms.begin(); bit < _atoms.end(); ++bit) {
                    // get center coordinates
                    const vec atom_pos = vec((*bit)->x * tools::conv::ang2bohr, (*bit)->y * tools::conv::ang2bohr, (*bit)->z * tools::conv::ang2bohr);


                    std::vector<double> temp;
                    // for each gridpoint
                    for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end(); ++git) {

                        temp.push_back(abs(git->grid_pos - atom_pos));

                    } // gridpoint of _atomgrid
                    rq.push_back(temp); // rq[center][gridpoint]

                } // centers
                // cout << " Calculated all gridpoint distances to centers for " << i_atom << endl;

                // find nearest-neighbor of this atom
                double distNN = 1e10;

                vector< ctp::QMAtom* > ::iterator NNit;
                //int i_NN;

                // now check all other centers
                int i_b = 0;
                for (bit = _atoms.begin(); bit != _atoms.end(); ++bit) {

                    if (bit != ait) {
                        // get center coordinates

                        const vec atomB_pos = vec((*bit)->x * tools::conv::ang2bohr, (*bit)->y * tools::conv::ang2bohr, (*bit)->z * tools::conv::ang2bohr);
                        double distSQ = (atomA_pos - atomB_pos)*(atomA_pos - atomB_pos);

                        // update NN distance and iterator
                        if (distSQ < distNN) {
                            distNN = distSQ;
                            NNit = bit;
                            //i_NN = i_b;
                        }

                    } // if ( ait != bit) 
                    i_b++;
                }// bit centers

                for (unsigned i_grid = 0; i_grid < _atomgrid.size(); i_grid++) {
                    //cout << " modifying point " << i_grid << endl;
                    // call some shit called grid_ssw0 in NWChem
                    std::vector<double> _p = SSWpartition(i_grid, _atoms.size(), rq);
                    //cout << " partition for gridpoint " << i_grid << endl;
                    // check weight sum
                    double wsum = 0.0;
                    for (unsigned i = 0; i < _p.size(); i++) {
                        wsum += _p[i];
                    }
                    //cout << " sum of partition weights " << wsum << endl;
                    if (wsum != 0.0) {

                        // update the weight of this grid point
                        _atomgrid[i_grid].grid_weight = _atomgrid[i_grid].grid_weight * _p[i_atom] / wsum;
                        //cout << " adjusting gridpoint weight "  << endl;
                    } else {

                        cerr << "\nSum of partition weights of grid point " << i_grid << " of atom " << i_atom << " is zero! ";
                        throw std::runtime_error("\nThis should never happen!");

                    }


                } // partition weight for each gridpoint

                // now remove points from the grid with negligible weights

                for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end();) {
                    if (git->grid_weight < 1e-13) {
                        git = _atomgrid.erase(git);
                    } else {
                        ++git;
                    }
                }



                // cout << " Total size of integration grid for atom: " << i_atom << " : " << _atomgrid.size() << " from " << fullsize << endl;

                _totalgridsize += _atomgrid.size();
                _grid.push_back(_atomgrid);
                i_atom++;
            } // atoms




            FindsignificantAtoms(basis);
            return;
        }

        std::vector<double> NumericalIntegration::SSWpartition(int igrid, int ncenters, std::vector< std::vector<double> >& rq) {
            const double ass = 0.725;
            // initialize partition vector to 1.0
            std::vector<double> p(ncenters, 1.0);

            const double tol_scr = 1e-10;
            const double leps = 1e-6;
            // go through centers
            for (int i = 1; i < ncenters; i++) {

                int ij = i * (i + 1) / 2 - 1; // indexing magic
                double rag = rq[i][igrid];

                // through all other centers (one-directional)
                for (int j = 0; j < i; j++) {

                    ij++;
                    if ((std::abs(p[i]) > tol_scr) || (std::abs(p[j]) > tol_scr)) {



                        double mu = (rag - rq[j][igrid]) * Rij[ij];
                        if (mu > ass) {
                            p[i] = 0.0;
                        } else if (mu < -ass) {
                            p[j] = 0.0;
                        } else {

                            double sk;
                            if (std::abs(mu) < leps) {
                                sk = -1.88603178008 * mu + 0.5;
                            } else {
                                sk = erf1c(mu);
                            }
                            if (mu > 0.0) sk = 1.0 - sk;
                            p[j] = p[j] * sk;
                            p[i] = p[i] * (1.0 - sk);

                        }
                    }
                }

            }

            return p;
        }

        double NumericalIntegration::erf1c(double x) {

            const static double alpha_erf1 = 1.0 / 0.30;
            return 0.5 * erfcc((x / (1.0 - x * x)) * alpha_erf1);
        }

        double NumericalIntegration::erfcc(double x) {

            double tau = 1.0 / (1.0 + 0.5 * std::abs(x));

            return tau * exp(-x * x - 1.26551223 + 1.00002368 * tau + 0.37409196 * tau * tau
                    + 0.09678418 * pow(tau, 3) - 0.18628806 * pow(tau, 4) + 0.27886807 * pow(tau, 5)
                    - 1.13520398 * pow(tau, 6) + 1.48851587 * pow(tau, 7) - 0.82215223 * pow(tau, 8)
                    + 0.17087277 * pow(tau, 9));
        }

    }
}
