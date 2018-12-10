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

#include <votca/ctp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/threecenter.h>
#include <votca/xtp/gwbse_exact.h>
#include <votca/xtp/rpa_spectral.h>
#include <votca/xtp/sigma_spectral.h>

#include <votca/xtp/numerical_integrations.h>
#include <votca/xtp/ppm.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/sigma.h>

namespace votca {
    namespace xtp {
        
        Eigen::MatrixXd GWBSE_Exact::CalculateVXC(const AOBasis& dftbasis) {
            bool _doVxc = true;
            std::string _functional = "XC_HYB_GGA_XC_PBEH";
            std::string _grid = "medium";
            Eigen::MatrixXd vxc_ao;
            if (_orbitals.hasAOVxc()) {
                if (_doVxc) {
                    if (_orbitals.getQMpackage() == "xtp") {
                        CTP_LOG(ctp::logDEBUG, *_pLog)
                                << ctp::TimeStamp()
                                << " Taking VXC from xtp DFT run."
                                << std::flush;
                    } else {
                        CTP_LOG(ctp::logDEBUG, *_pLog)
                                << ctp::TimeStamp()
                                << " There is already a Vxc matrix loaded from DFT, did you maybe "
                                "run a DFT code with outputVxc?\n I will take the external "
                                "implementation"
                                << std::flush;
                    }
                    _doVxc = false;
                }
                CTP_LOG(ctp::logDEBUG, *_pLog)
                        << ctp::TimeStamp()
                        << " Loaded external Vxc matrix"
                        << std::flush;
                vxc_ao = _orbitals.AOVxc();
            } else if (_doVxc) {
                NumericalIntegration numint;
                numint.setXCfunctional(_functional);
                double ScaHFX_temp = numint.getExactExchange(_functional);
                if (ScaHFX_temp != _orbitals.getScaHFX()) {
                    throw std::runtime_error(
                            (boost::format("GWBSE exact exchange a=%s differs from qmpackage "
                            "exact exchange a=%s, probably your functionals are "
                            "inconsistent") %
                            ScaHFX_temp % _orbitals.getScaHFX())
                            .str());
                }
                numint.GridSetup(_grid, _orbitals.QMAtoms(), dftbasis);
                CTP_LOG(ctp::logDEBUG, *_pLog)
                        << ctp::TimeStamp()
                        << " Setup grid for integration with gridsize: "
                        << _grid
                        << " with "
                        << numint.getGridSize()
                        << " points, divided into "
                        << numint.getBoxesSize()
                        << " boxes"
                        << std::flush;
                CTP_LOG(ctp::logDEBUG, *_pLog)
                        << ctp::TimeStamp() << " Integrating Vxc in VOTCA with functional "
                        << _functional
                        << std::flush;
                Eigen::MatrixXd DMAT = _orbitals.DensityMatrixGroundState();
                vxc_ao = numint.IntegrateVXC(DMAT);
                CTP_LOG(ctp::logDEBUG, *_pLog)
                        << ctp::TimeStamp()
                        << " Calculated Vxc in VOTCA"
                        << std::flush;
            } else {
                throw std::runtime_error(
                        "So your DFT data contains no Vxc, if you want to proceed use the "
                        "dovxc option.");
            }
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Set hybrid exchange factor: "
                    << _orbitals.getScaHFX()
                    << std::flush;
            // now get expectation values but only for those in _qpmin:_qpmax range
            Eigen::MatrixXd _mos = _orbitals.MOCoefficients().block(0, _qpmin, _orbitals.MOCoefficients().rows(), _qptotal);
            Eigen::MatrixXd vxc = _mos.transpose() * vxc_ao*_mos;
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Calculated exchange-correlation expectation values "
                    << std::flush;
            return vxc;
        }

        bool GWBSE_Exact::Evaluate() {
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Running GWBSE_Exact :D "
                    << std::flush;

            // ****** Check input ******
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Molecule Coordinates [A] "
                    << std::flush;

            for (QMAtom* atom : _orbitals.QMAtoms()) {

                std::string output = (boost::format("  %1$s"
                        "   %2$+1.4f %3$+1.4f %4$+1.4f")
                        % (atom->getType())
                        % (atom->getPos().getX() * tools::conv::bohr2ang)
                        % (atom->getPos().getY() * tools::conv::bohr2ang)
                        % (atom->getPos().getZ() * tools::conv::bohr2ang)).str();
                
                CTP_LOG(ctp::logDEBUG, *_pLog)
                        << ctp::TimeStamp()
                        << output
                        << std::flush;
                
            }

            // Check which QC program was used for the DFT run
            // -> Implicit info about MO coefficient storage order
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " DFT data was created by "
                    << _orbitals.getQMpackage()
                    << std::flush;
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " homo: "
                    << _homo
                    << std::flush;
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " RPA indices: "
                    << "[" << _rpamin << ", " << _rpamax << "]"
                    << std::flush;
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " QP indices: "
                    << "[" << _qpmin << ", " << _qpmax << "]"
                    << std::flush;
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " BSE indices:"
                    << " v: [" << _bse_vmin << ", " << _bse_vmax << "]"
                    << " c: [" << _bse_cmin << ", " << _bse_cmax << "]"
                    << std::flush;

            // ****** Configure orbitals object ******
            
            _orbitals.setRPAindices(_rpamin, _rpamax);
            _orbitals.setGWAindices(_qpmin, _qpmax);
            _orbitals.setBSEindices(_bse_vmin, _bse_cmax, _bse_maxeigenvectors);
            
            // ****** Set TDA approximation ******

            // ****** Load DFT basis ******
            
            BasisSet dftbs;
            dftbs.LoadBasisSet(_dftbasis_name);
            _orbitals.setDFTbasis(_dftbasis_name);

            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Loaded DFT Basis Set "
                    << _dftbasis_name
                    << std::flush;

            // ****** Fill DFT AO basis ******
            
            const int fragA = -1;

            // Fill DFT AO basis by going through all atoms
            AOBasis dftbasis;
            dftbasis.AOBasisFill(dftbs, _orbitals.QMAtoms(), fragA);
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Filled DFT Basis of size "
                    << dftbasis.AOBasisSize()
                    << std::flush;
            
            if (dftbasis.getAOBasisFragB() > 0 && dftbasis.getAOBasisFragA() > 0) {
                
                CTP_LOG(ctp::logDEBUG, *_pLog)
                        << ctp::TimeStamp() << " FragmentA size "
                        << dftbasis.getAOBasisFragA()
                        << std::flush;
                
                CTP_LOG(ctp::logDEBUG, *_pLog)
                        << ctp::TimeStamp()
                        << " FragmentB size "
                        << dftbasis.getAOBasisFragB()
                        << std::flush;
            }
            
            // ****** Load aux. basis ******

            // Load auxiliary basis set (element-wise information) from xml file
            BasisSet auxbs;
            auxbs.LoadBasisSet(_auxbasis_name);
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Loaded Auxbasis Set "
                    << _auxbasis_name
                    << std::flush;

            // Fill auxiliary AO basis by going through all atoms
            AOBasis auxbasis;
            auxbasis.AOBasisFill(auxbs, _orbitals.QMAtoms());
            _orbitals.setAuxbasis(_auxbasis_name);
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Filled Auxbasis of size "
                    << auxbasis.AOBasisSize()
                    << std::flush;
            
            // ****** Fill Exchange-Correlation Potential (VXC) matrix ******
            
            Eigen::MatrixXd vxc=CalculateVXC(dftbasis);

            // ****** Fill overlap matrix ******

            /*
             * For the representation of 2-point functions with the help of the
             * auxiliary basis, its AO overlap matrix is required.
             * cf. M. Rohlfing, PhD thesis, ch. 3
             */
            AOOverlap auxoverlap;
            auxoverlap.Fill(auxbasis);

            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Filled Aux Overlap matrix of dimension: "
                    << auxoverlap.Matrix().rows()
                    << std::flush;
            
            // ****** Fill Coulomb matrix ******

            /*
             *  For the calculation of Coulomb and exchange term in the self
             *  energy and electron-hole interaction, the Coulomb interaction
             *  is represented using the auxiliary basis set.
             *  Here, we need to prepare the Coulomb matrix expressed in
             *  the AOs of the auxbasis
             */

            // Get Coulomb matrix as AOCoulomb
            AOCoulomb auxcoulomb;
            auxcoulomb.Fill(auxbasis);
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Filled Aux Coulomb matrix of dimension: "
                    << auxcoulomb.Matrix().rows()
                    << std::flush;
            
            Eigen::MatrixXd Coulomb_sqrtInv = auxcoulomb.Pseudo_InvSqrt_GWBSE(auxoverlap, 5e-7);
            
            auxoverlap.FreeMatrix();
            auxcoulomb.FreeMatrix();
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Calculated Matrix Sqrt of Aux Coulomb Matrix"
                    << std::flush;
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Removed "
                    << auxcoulomb.Removedfunctions()
                    << " functions from Aux Coulomb matrix to avoid near linear dependencies"
                    << std::flush;
            
            // ****** Prepare 3-center integral object ******

            // Container => M_mn
            // Prepare 3-center integral object
            TCMatrix_gwbse Mmn;
            Mmn.Initialize(auxbasis.AOBasisSize(), _rpamin, _qpmax, _rpamin, _rpamax); //rpamin here, because RPA needs till rpamin
            Mmn.Fill(auxbasis, dftbasis, _orbitals.MOCoefficients());
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Calculated Mmn_beta (3-center-repulsion x orbitals)  "
                    << std::flush;

            // Make _Mmn symmetric
            Mmn.MultiplyRightWithAuxMatrix(Coulomb_sqrtInv);
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Prepared Mmn "
                    << std::flush;
            
            // ****** gwbse.cc ******

            PPM ppm;
            RPA rpa;
            rpa.configure(_homo, _rpamin, _rpamax);
            Eigen::VectorXd screen_r = Eigen::VectorXd::Zero(1);
            screen_r(0) = ppm.getScreening_r();
            Eigen::VectorXd screen_i = Eigen::VectorXd::Zero(1);
            screen_i(0) = ppm.getScreening_i();
            rpa.setScreening(screen_r, screen_i);
            Sigma sigma=Sigma(_pLog);
            sigma.configure(_homo,_qpmin,_qpmax,0,0);
            sigma.setDFTdata(_orbitals.getScaHFX(),&vxc,&_orbitals.MOEnergies());
            Eigen::VectorXd gwa_energies = Eigen::VectorXd::Zero(_orbitals.getNumberOfLevels());
            double _shift = 0; // No shift?
            for (int i = 0; i < gwa_energies.size(); ++i) {
              gwa_energies(i) = _orbitals.MOEnergies()(i);
              if (i > int(_homo)) {
                gwa_energies(i) += _shift;
              }
            }
            sigma.setGWAEnergies(gwa_energies);
            rpa.calculate_epsilon(gwa_energies,Mmn);
            ppm.PPM_construct_parameters(rpa);
            Mmn.MultiplyRightWithAuxMatrix(ppm.getPpm_phi());
            sigma.CalcdiagElements(Mmn,ppm);
            sigma.CalcOffDiagElements(Mmn,ppm);
            Eigen::MatrixXd Hqp = sigma.SetupFullQPHamiltonian();

            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Computed sigma_c using PPM "
                    << std::flush;
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " sigma_c:\n"
                    << sigma.get_sigma_c()
                    << std::flush;
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Hqp:\n"
                    << Hqp
                    << std::flush;
            
            // ****** Prepare RPA ******
            
            RPA_Spectral rpa_spectral = RPA_Spectral();
            rpa_spectral.configure_bse(_homo, _bse_vmin, _bse_cmax);
            rpa_spectral.configure_qp(_homo, _qpmin, _qpmax);
            rpa_spectral.set_GWAEnergies(_orbitals.MOEnergies());
            rpa_spectral.prepare_decomp(Mmn);
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Prepared RPA Spectral "
                    << std::flush;
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Eigenvalues:\n"
                    << rpa_spectral.get_Omega().transpose()
                    << std::flush;
            
            // ****** Prepare Sigma ******

            Sigma_Spectral sigma_spectral = Sigma_Spectral();
            sigma_spectral.configure_bse(_homo, _bse_vmin, _bse_cmax);
            sigma_spectral.configure_qp(_homo, _qpmin, _qpmax);
            sigma_spectral.configure_g_iter(40, 1e-5);
            sigma_spectral.set_GWAEnergies(_orbitals.MOEnergies());
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Prepared Sigma Spectral "
                    << std::flush;
            
            // ****** Compute Sigma_c ******

            sigma_spectral.setHedin(false);
            sigma_spectral.compute_sigma(Mmn, rpa_spectral, _orbitals.getScaHFX());
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Computed sigma_c exactly "
                    << std::flush;
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " sigma_c:\n"
                    << sigma_spectral.get_sigma_c()
                    << std::flush;
            
            /*
            sigma_spectral.setHedin(true);
            sigma_spectral.compute_sigma(Mmn, rpa_spectral, _orbitals.getScaHFX());
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Computed sigma_c with Hedin's approx. "
                    << std::flush;
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " sigma_c:\n"
                    << sigma_spectral.get_sigma_c()
                    << std::flush;
            */
            
            // ****** Set-up Hamiltonian
            
            Eigen::MatrixXd Hqp_spectral = sigma_spectral.SetupFullQPHamiltonian(vxc);
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Hqp:\n"
                    << Hqp_spectral
                    << std::flush;

            CTP_LOG(ctp::logDEBUG, *_pLog) << std::flush;
            return false;
        }

    }
}