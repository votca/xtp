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
//for libxc
#include <votca/xtp/votca_config.h>

#include <votca/xtp/numerical_integrations_periodic.h>
#include <boost/math/constants/constants.hpp>
#include <votca/xtp/aoshell.h>
#include <votca/tools/constants.h>

#ifdef LIBXC
#include <xc.h>
#endif
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/elements.h>
#include <fstream>
#include <iterator>
#include <string>




namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        double NumericalIntegrationPeriodic::IntegratePotential_w_PBC(vec rvector) {
            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }

            double result = 0.0;
            double cutoff = min(min(boxLen[0], boxLen[1]), boxLen[2]) / 2.0;
            double vol = boxLen[0] * boxLen[1] * boxLen[2];

            if (density_set) {
                //                #pragma omp parallel for reduction(+:result)
                for (unsigned i = 0; i < _grid.size(); i++) {
                    for (unsigned j = 0; j < _grid[i].size(); j++) {

                        //charge at this point
                        double q = -_grid[i][j].grid_weight * _grid[i][j].grid_density; //density is neg of charge

                        //r-space sum
                        vec dif = _grid[i][j].grid_pos - rvector;

                        for (int k = 0; k < 3; k++) {
                            if (std::abs(dif[k] > boxLen[k]*0.5)) //correct for for bond crossing PBC, if it exists
                                if (dif[k] > 0) //i.x>j.x
                                    dif[k] -= boxLen[k];
                                else //i.x<j.x
                                    dif[k] += boxLen[k];
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
                int nKpoints = numK[0] * numK[1] * numK[2];
                double ksum = 0.0;
                //exclude k={0,0,0} (index=0)
                //                #pragma omp parallel for reduction(+:ksum)
                for (int index = 1; index < nKpoints; index++) {
                    double* kp = &(Kcoord[index * 3]);
                    vec K(kp[0], kp[1], kp[2]);
                    std::complex<double> Kr(0, K * rvector); // ik dot r

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

        double NumericalIntegrationPeriodic::CalcDipole_w_PBC(vec rvector) {
            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }
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

        void NumericalIntegrationPeriodic::FreeKspace() {
            delete[] Rho_k;
            delete[] Kcoord;
            delete[] prefactor;
        }

        /**
         *	Calculate and return the Ewald coefficient (alpha) from cutoff distance
         *	and requested tolerance.
         */
        void NumericalIntegrationPeriodic::findAlpha(double Rc, double dtol) {
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
        void NumericalIntegrationPeriodic::PrepKspaceDensity(double ext_alpha, std::vector< ctp::QMAtom* > & _local_atomlist, bool ECP, int nK = 0) {

            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }
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
            //cout << "numK={" << numK[0] << ", " << numK[1] << ", " << numK[2] << "}" << endl;
            //cout << "rel err of k-sum ~ " << std::exp(-minSq / fourasq) / minSq << endl;


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

        void NumericalIntegrationPeriodic::PrepKspaceDensity_gromacs_like(double ext_alpha, std::vector< ctp::QMAtom* > & _local_atomlist, bool ECP, Grid &eval_grid, int nK) {
            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }

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
            //cout << "numK={" << numK[0] << ", " << numK[1] << ", " << numK[2] << "}" << endl;
            //cout << "rel err of k-sum ~ " << std::exp(-minSq / fourasq) / minSq << endl;


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

                //also compute the inverse box dimensions
                lll[m] = 2.0 * boost::math::constants::pi<double>() / boxLen[m]; //units are 1/Bohr
            }
            //eikr indexing is [DIM][k][i][j]

            //based on gromacs 4.6 tabulate_eir()
            for (unsigned i = 0; i < _grid.size(); i++) {
#pragma omp parallel for
                for (unsigned j = 0; j < _grid[i].size(); j++) {
                    vec r = _grid[i][j].grid_pos; //grid is in Bohr, need it in reduced units (*2pi/box))
                    for (unsigned m = 0; (m < 3); m++) {
                        eikr[m][0][i][j] = std::complex<double>(1, 0);
                        eikr[m][1][i][j] = std::complex<double>(cos(r[m] * lll[m]), sin(r[m] * lll[m])); //this is where reduced units are applied
                        for (unsigned k = 2; k < numK[m]; k++) {
                            eikr[m][k][i][j] = eikr[m][k - 1][i][j] * eikr[m][1][i][j];
                        }//k
                    }//m
                }//j
            }//i


            //now do the same for the eval_grid (eikR), the points were potential is to be evaluated
            //allocate space for eikR
            eikR.resize(3); //3 dim
            for (unsigned m = 0; m < 3; m++) {
                eikR[m].resize(numK[m]);
                for (unsigned k = 0; k < numK[m]; k++) {
                    eikR[m][k].resize(eval_grid.getGrid().size());
                }
            }
            //eikR indexing is [DIM][k][i]

            //based on gromacs 4.6 tabulate_eir()
#pragma omp parallel for
            //printf("eval grid size = %d\n", eval_grid.getGrid().size());
            //printf("lll_x = %f\n", lll[0]);
            for (unsigned i = 0; i < eval_grid.getGrid().size(); i++) {
                vec R = eval_grid.getGrid()[i] * tools::conv::nm2bohr; //this is in Bohr, eval_grid is in nm
                for (unsigned m = 0; (m < 3); m++) {
                    eikR[m][0][i] = std::complex<double>(1, 0);
                    eikR[m][1][i] = std::complex<double>(cos(R[m] * lll[m]), sin(R[m] * lll[m]));
                    //unsigned k=0;
                    //printf("%d\t%f+i%f\t R_x=%f\n", k, eikR[m][k][i].real(),eikR[m][k][i].imag(), R[m]);
                    //k=1;
                    //printf("%d\t%f+i%f\t R_x=%f\n", k, eikR[m][k][i].real(),eikR[m][k][i].imag(), R[m]);
                    for (unsigned k = 2; k < numK[m]; k++) {
                        eikR[m][k][i] = eikR[m][k - 1][i] * eikR[m][1][i];
                        //printf("%d\t%f+i%f\t R_x=%f\n", k, eikR[m][k][i].real(),eikR[m][k][i].imag(), R[m]);
                    }//k
                    //exit(0);
                }//m
            }//i

        }

        void NumericalIntegrationPeriodic::FillMadelungGrid(vec box, int natomsonside) {
            //fill _Madelung_grid;
            _Madelung_grid.clear();
            std::vector< GridContainers::integration_grid > _Mad;
            boxLen = box;
            double a = boxLen[0]; //in bohr
            for (int l = 0; l < 2; l++) {
                for (int m = 0; m < 2; m++) {
                    for (int n = 0; n < 2; n++) {
                        GridContainers::integration_grid el;
                        el.grid_density = std::pow(-1.0, l + m + n);
                        el.grid_weight = 1.0;
                        el.grid_pos = vec(l, m, n) * (a / 2);
                        _Mad.push_back(el);
                        //cout<< "q= "<< el.grid_density << "\t @ " << el.grid_x << " " << el.grid_y << " " << el.grid_z << "\t box size:"<< boxLen[0] << endl;
                    }
                }
            }
            _Madelung_grid.push_back(_Mad);
            _grid = _Madelung_grid;
            //exit(0);
        }

        void NumericalIntegrationPeriodic::IntegratePotential_w_PBC_gromacs_like(Grid &eval_grid, ub::vector<double>& _ESPatGrid) {
            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }
            
            double cutoff = min(min(boxLen[0], boxLen[1]), boxLen[2]) / 2.0; //Bohr
            double vol = boxLen[0] * boxLen[1] * boxLen[2]; //Bohr^3
            //eval_grid is in nm


            ub::vector<double> k_pot;
            k_pot.resize(_ESPatGrid.size());
            for (unsigned p = 0; p < k_pot.size(); p++) {
                k_pot[p] = 0;
            }


            if (density_set) {
                //#pragma omp parallel for
                for (unsigned p = 0; p < _ESPatGrid.size(); p++) {
                    _ESPatGrid[p] = 0;
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
                E_rspace = _ESPatGrid(0);
                E_erfc = erfc(alpha);


                //adapted from gromacs 4.6 do_ewald_pot()

                //zero pot
                //                for (unsigned p = 0; p < _ESPatGrid.size(); p++) {
                //                    _ESPatGrid[p] = 0;
                //                }



                int ix, iy, iz;
                double tmp, cs, ss, ak, mx, my, mz, m2;
                double factor = -1.0 / (4 * alpha * alpha);
                double scaleRecip = 4.0 * boost::math::constants::pi<double>() / (vol); //final units of potential are Hartree/e (atomic units)
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

                //printf("factor = %f\n", factor);
                //printf("scaleRecip = %f\n", scaleRecip);
                //printf("vol = %f\n", vol);
                //printf("actual box for Madelung: %f\t%f\t%f\n", boxLen[0], boxLen[1], boxLen[2]);
                //exit(0);


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
                            m2 = mx * mx + my * my + mz * mz;
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
                                    //                                    if(ix==4 && iy==4 && iz==4){
                                    //                                        printf("cs+= %f = %f * (%f+i%f) * (%f+i%f)\n", tab_qxyz[i][j].real(), -_grid[i][j].grid_weight * _grid[i][j].grid_density, tab_xy[i][j].real(), tab_xy[i][j].imag(), eikr[2][iz][i][j].real(), eikr[2][iz][i][j].imag());
                                    //                                        printf("\teikr=(%f+i%f) * (%f+i%f) * (%f+i%f)\n",eikr[0][ix][i][j].real(), eikr[0][ix][i][j].imag(),eikr[1][iy][i][j].real(), eikr[1][iy][i][j].imag(),eikr[2][iz][i][j].real(), eikr[2][iz][i][j].imag());
                                    //                                        printf("\tr=%f\t%f\t%f\n", _grid[i][j].grid_pos[0],_grid[i][j].grid_pos[1], _grid[i][j].grid_pos[2]);
                                    //                                        std::complex<double> krx = std::complex<double>(0, -_grid[i][j].grid_pos[0]*mz);
                                    //                                        krx=exp(krx);
                                    //                                        printf("\tkr_x=exp(-i%f)=%f+i%f\n",_grid[i][j].grid_pos[0]*mz, krx.real(), krx.imag());
                                    //                                    }
                                }
                            }
                            //exit(0);
                            //                            if(abs(cs)>1e-6 || abs(ss)>1e-6)
                            //                            if(ix==4 && iy==4 && iz==4)
                            //                                printf("k=[%d,%d,%d]\tcs=%f\t ss=%f \tak=%f=exp(%f * %f) / %f\n", ix, iy, iz, cs, ss, ak, m2, factor, m2);
                            for (int n = 0; n < eval_grid.getsize(); n++) {
                                tmp = ak * (cs * tab_R_xyz[n].real() + ss * tab_R_xyz[n].imag());
                                //                                if((abs(cs)>1e-6 || abs(ss)>1e-6) && n==0){
                                //                                    printf("ak=%f=exp(%f * %f) / %f\n", ak, m2, factor, m2);
                                //                                    printf("tab_R_xyz=%f+i%f= (%f+i%f) * (%f+i%f) * (%f+i%f)\n", tab_R_xyz[n].real(), tab_R_xyz[n].imag(), eikR[0][ix][n].real(), eikR[0][ix][n].imag(), eikR[1][abs(iy)][n].real(), eikR[1][abs(iy)][n].imag(), eikR[2][abs(iz)][n].real(), eikR[2][abs(iz)][n].imag());
                                //                                    printf("\ttmp=%f=%f*(%f*%f + %f*%f)\t m2=%f\n", tmp, ak, cs, tab_R_xyz[n].real(), ss, tab_R_xyz[n].imag(), m2);
                                //                                    printf("\tk=[%d,%d,%d] of [%d,%d,%d]\ttmp=%f \t sum=%f\n", ix, iy, iz, numK[0], numK[1], numK[2], tmp, k_pot(n));
                                //                                }
                                _ESPatGrid(n) += tmp * 2 * scaleRecip;
                                k_pot(n) += tmp * 2 * scaleRecip;

                            }

                        }//iz
                        lowiz = 1 - numK[2];
                        //                        if(ix==1)
                        //                            exit(0);                        
                    }//iy
                    lowiy = 1 - numK[1];
                }//ix
                //printf("%f=%f+%f\n", _ESPatGrid(0), E_rspace, k_pot(0));
                //                exit(0);


            } else {
                throw std::runtime_error("Density not calculated");
            }
            E_kspace = k_pot(0);
            //E_kspace = _ESPatGrid(1);
            //_ESPatGrid(1)+=E_rspace;
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
        double NumericalIntegrationPeriodic::IntegrateDensity_Molecule(ub::matrix<double>& _density_matrix, AOBasis* basis, std::vector<int> AtomIndeces) {




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

                        // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )
                        ub::matrix<double> AOgrid = ub::zero_matrix<double>(1, basis->AOBasisSize()); //basis function values at grid for all significant atoms
                        ub::matrix<double> rho_mat = ub::zero_matrix<double>(1, 1);

                        // evaluate AO Functions for all shells, NOW BLOCKWISE

                        // for each significant atom for this grid point
                        for (unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size(); sigrow++) {

                            // this atom
                            int rowatom = _significant_atoms[i][j][sigrow];

                            // for each shell in this atom
                            for (unsigned ishell = 0; ishell < _atomshells[rowatom].size(); ishell++) {

                                AOBasis::AOShellIterator _row = _atomshells[rowatom][ishell];
                                AOShell* _shell = *_row;
                                // for density, fill sub-part of AOatgrid
                                ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, 0, 1, _shell->getStartIndex(), _shell->getStartIndex() + _shell->getNumFunc());
                                _shell->EvalAOspace(_AOgridsub, _grid[i][j].grid_pos);

                            } // shell in atom
                        }

                        for (unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size(); sigrow++) {

                            // this atom
                            int rowatom = _significant_atoms[i][j][sigrow];
                            //only do things if this atom belongs to the molecule of interest and is significant
                            if (std::find(AtomIndeces.begin(), AtomIndeces.end(), rowatom) != AtomIndeces.end()) {

                                ub::matrix<double> _temp = ub::zero_matrix<double>(1, _blocksize[rowatom]);

                                ub::matrix_range< ub::matrix<double> > _AOgridrow = ub::subrange(AOgrid, 0, 1, _startIdx[rowatom], _startIdx[rowatom] + _blocksize[rowatom]);

                                // for each atom

                                for (unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size(); sigcol++) {
                                    int colatom = _significant_atoms[i][j][sigcol];
                                    ub::matrix_range< ub::matrix<double> > _AOgridcol = ub::subrange(AOgrid, 0, 1, _startIdx[colatom], _startIdx[colatom] + _blocksize[colatom]);
                                    ub::matrix_range< ub::matrix<double> > DMAT_here = ub::subrange(_density_matrix, _startIdx[colatom], _startIdx[colatom] + _blocksize[colatom], _startIdx[rowatom], _startIdx[rowatom] + _blocksize[rowatom]);

                                    _temp += ub::prod(_AOgridcol, DMAT_here);


                                } //col shells

                                rho_mat += ub::prod(_temp, ub::trans(_AOgridrow));
                            } // atom belongs to molecule

                        } // row shells


                        _grid[i][j].grid_density = rho_mat(0, 0);
                        Density_thread[i_thread] += _grid[i][j].grid_weight * _grid[i][j].grid_density;
                    } // j: for each point in atom grid
                }// each thread
            } // i: for each atom grid

            double result = 0.0;
            for (int i_thread = 0; i_thread < nthreads; i_thread++) {
                //cout << result << endl;
                result += Density_thread[i_thread];
            }
            density_set = true;

            return (result);
        }

        /*
         * Finds distances between integration centers.
         * Overloads the same function from the parent
         */
        void NumericalIntegrationPeriodic::FindCenterCenterDist(vector<ctp::QMAtom*> _atoms){
        
            if(boxLen*boxLen==0){
                throw std::runtime_error("NumericalIntegrationPeriodic: periodic box not set.");
            }
            
            int ij = 0;
            Rij.push_back(0.0); // 1st center "self-distance"
            
            vector< ctp::QMAtom* > ::iterator ait;
            vector< ctp::QMAtom* > ::iterator bit;
            int i = 1;
            for (ait = _atoms.begin() + 1; ait != _atoms.end(); ++ait) {
                // get center coordinates in Bohr
                vec pos_a = (*ait)->getPos() * tools::conv::ang2bohr;
                
                int j = 0;
                for (bit = _atoms.begin(); bit != ait; ++bit) {
                    ij++;
                    // get center coordinates in Bohr
                    vec pos_b = (*bit)->getPos() * tools::conv::ang2bohr;
                    vec dif = WrapDisplacement(pos_a, pos_b, boxLen);
                    Rij.push_back(1.0 / abs(dif));
                                        
                    j++;
                } // atoms
                Rij.push_back(0.0); // self-distance again
                i++;
            } // atoms
            return;
        }
        
        
        /**
         * \brief Wraps a vector to periodic boundary conditions. This eventually needs to be moved to tools::vec.
         * 
         * @param vector r to wrap
         * @param vector box specifying the periodic box
         * @return the wrapped vector
         */
        tools::vec WrapPoint(const tools::vec r, const tools::vec box) {
            tools::vec ret;
            for(int k=0; k<3; k++)
            {
                ret[k] = fmod(r[k], box[k]);
            }
            return(ret);
        }
        
        /**
         * \brief Finds the shortest displacement between images of two vectors due to periodic boundary conditions.
         * This eventually needs to be moved to tools::vec.
         * 
         * @param vector a
         * @param vector b
         * @param vector box specifying the periodic box
         * @return the shortest displacement
         */
        tools::vec WrapDisplacement(const tools::vec a, const tools::vec b , const tools::vec box) {
            tools::vec ret = a-b;
            for(int k=0; k<3; k++)
            {
                ret[k] = fmod(ret[k], box[k]);
                if(ret[k] > box[k]*0.5){
                    ret[k] = box[k]-ret[k];
                }
            }
            return(ret);
        }
        
    }//xtp
}//votca
