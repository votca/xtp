/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Standard includes
#include <vector>

// Third party includes
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/aomatrix3d.h"
#include "votca/xtp/aopotential.h"
#include "votca/xtp/density_integration.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/mmregion.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/qmmolecule.h"

using boost::format;
using namespace boost::filesystem;
using namespace std;
using std::flush;
using namespace votca::tools;

namespace votca {
namespace xtp {

void DFTEngine::Initialize(Property& options) {

  string key = "package";
  const string key_xtpdft = "package.xtpdft";
  _dftbasis_name = options.get(key + ".basisset").as<string>();

  if (options.get(key + ".use_auxbasisset").as<bool>()) {
    _auxbasis_name = options.get(key + ".auxbasisset").as<string>();
  }

  _four_center_method =
      options.get(key_xtpdft + ".four_center_method").as<std::string>();

  if (_four_center_method != "RI") {
    _with_screening = options.get(key_xtpdft + ".with_screening").as<bool>();
    _screening_eps = options.get(key_xtpdft + ".screening_eps").as<double>();
  }

  if (options.get(key + ".use_ecp").as<bool>()) {
    _ecp_name = options.get(key + ".ecp").as<string>();
    _with_ecp = true;
  } else {
    _with_ecp = false;
  }
  _with_guess = options.get(key + ".read_guess").as<bool>();
  _initial_guess = options.get(key_xtpdft + ".initial_guess").as<string>();

  _grid_name = options.get(key_xtpdft + ".integration_grid").as<string>();
  _xc_functional_name = options.get(key + ".functional").as<string>();

  if (options.get(key_xtpdft + ".use_external_density").as<bool>()) {
    _integrate_ext_density = true;
    _orbfilename = options.ifExistsReturnElseThrowRuntimeError<string>(
        key_xtpdft + ".externaldensity.orbfile");
    _gridquality = options.ifExistsReturnElseThrowRuntimeError<string>(
        key_xtpdft + ".externaldensity.gridquality");
    _state = options.ifExistsReturnElseThrowRuntimeError<string>(
        key_xtpdft + ".externaldensity.state");
  }

  if (options.get(key_xtpdft + ".use_external_field").as<bool>()) {
    _integrate_ext_field = true;

    _extfield = options.ifExistsReturnElseThrowRuntimeError<Eigen::Vector3d>(
        key_xtpdft + ".externalfield");
  }

  _conv_opt.Econverged =
      options.get(key_xtpdft + ".convergence.energy").as<double>();
  _conv_opt.error_converged =
      options.get(key_xtpdft + ".convergence.error").as<double>();
  _max_iter =
      options.get(key_xtpdft + ".convergence.max_iterations").as<Index>();

  string method = options.get(key_xtpdft + ".convergence.method").as<string>();
  if (method == "DIIS") {
    _conv_opt.usediis = true;
  } else if (method == "mixing") {
    _conv_opt.usediis = false;
  }
  if (!_conv_opt.usediis) {
    _conv_opt.histlength = 1;
    _conv_opt.maxout = false;
  }
  _conv_opt.mixingparameter =
      options.get(key_xtpdft + ".convergence.mixing").as<double>();
  _conv_opt.levelshift =
      options.get(key_xtpdft + ".convergence.levelshift").as<double>();
  _conv_opt.levelshiftend =
      options.get(key_xtpdft + ".convergence.levelshift_end").as<double>();
  _conv_opt.maxout =
      options.get(key_xtpdft + ".convergence.DIIS_maxout").as<bool>();
  _conv_opt.histlength =
      options.get(key_xtpdft + ".convergence.DIIS_length").as<Index>();
  _conv_opt.diis_start =
      options.get(key_xtpdft + ".convergence.DIIS_start").as<double>();
  _conv_opt.adiis_start =
      options.get(key_xtpdft + ".convergence.ADIIS_start").as<double>();

  return;
}

void DFTEngine::PrintMOs(const Eigen::VectorXd& MOEnergies, Log::Level level) {
  XTP_LOG(level, *_pLog) << "  Orbital energies: " << flush;
  XTP_LOG(level, *_pLog) << "  index occupation energy(Hartree) " << flush;
  for (Index i = 0; i < MOEnergies.size(); i++) {
    Index occupancy = 0;
    if (i < _numofelectrons / 2) {
      occupancy = 2;
    }
    XTP_LOG(level, *_pLog) << (boost::format(" %1$5d      %2$1d   %3$+1.10f") %
                               i % occupancy % MOEnergies(i))
                                  .str()
                           << flush;
  }
  return;
}

void DFTEngine::CalcElDipole(const Orbitals& orb) const {
  QMState state = QMState("n");
  Eigen::Vector3d result = orb.CalcElDipole(state);
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Electric Dipole is[e*bohr]:\n\t\t dx=" << result[0]
      << "\n\t\t dy=" << result[1] << "\n\t\t dz=" << result[2] << flush;
  return;
}

Mat_p_Energy DFTEngine::CalcEXXs(const Eigen::MatrixXd& MOCoeff,
                                 const Eigen::MatrixXd& Dmat) const {
  if (_four_center_method == "RI") {
    if (_conv_accelerator.getUseMixing() || MOCoeff.rows() == 0) {
      return _ERIs.CalculateEXX(Dmat);
    } else {
      Eigen::MatrixXd occblock =
          MOCoeff.block(0, 0, MOCoeff.rows(), _numofelectrons / 2);
      return _ERIs.CalculateEXX(occblock, Dmat);
    }
  } else {
    if (_four_center_method == "direct") {
      throw std::runtime_error(
          "direct 4c method only works with LDA and GGA functionals.");
    }
    return _ERIs.CalculateEXX_4c_small_molecule(Dmat);
  }
}

tools::EigenSystem DFTEngine::IndependentElectronGuess(
    const Mat_p_Energy& H0) const {
  return _conv_accelerator.SolveFockmatrix(H0.matrix());
}

tools::EigenSystem DFTEngine::ModelPotentialGuess(
    const Mat_p_Energy& H0, const QMMolecule& mol,
    const Vxc_Potential<Vxc_Grid>& vxcpotential) const {
  Eigen::MatrixXd Dmat = AtomicGuess(mol);
  Mat_p_Energy ERIs = CalculateERIs(Dmat);
  Mat_p_Energy e_vxc = vxcpotential.IntegrateVXC(Dmat);
  XTP_LOG(Log::info, *_pLog)
      << TimeStamp() << " Filled DFT Vxc matrix " << flush;

  Eigen::MatrixXd H = H0.matrix() + ERIs.matrix() + e_vxc.matrix();
  if (_ScaHFX > 0) {
    Mat_p_Energy EXXs = CalcEXXs(Eigen::MatrixXd::Zero(0, 0), Dmat);
    H -= 0.5 * _ScaHFX * EXXs.matrix();
  }
  return _conv_accelerator.SolveFockmatrix(H);
}

bool DFTEngine::Evaluate(Orbitals& orb) {
  Prepare(orb.QMAtoms());
  Mat_p_Energy H0 = SetupH0(orb.QMAtoms());
  tools::EigenSystem MOs;
  MOs.eigenvalues() = Eigen::VectorXd::Zero(H0.cols());
  MOs.eigenvectors() = Eigen::MatrixXd::Zero(H0.rows(), H0.cols());
  Vxc_Potential<Vxc_Grid> vxcpotential = SetupVxc(orb.QMAtoms());
  ConfigOrbfile(orb);

  if (_with_guess) {
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Reading guess from orbitals object/file" << flush;
    MOs = orb.MOs();
    MOs.eigenvectors() = OrthogonalizeGuess(MOs.eigenvectors());
  } else {
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Setup Initial Guess using: " << _initial_guess
        << flush;
    if (_initial_guess == "independent") {
      MOs = IndependentElectronGuess(H0);
    } else if (_initial_guess == "atom") {
      MOs = ModelPotentialGuess(H0, orb.QMAtoms(), vxcpotential);
    } else {
      throw std::runtime_error("Initial guess method not known/implemented");
    }
  }

  Eigen::MatrixXd Dmat = _conv_accelerator.DensityMatrix(MOs);
  XTP_LOG(Log::info, *_pLog)
      << TimeStamp() << " Guess Matrix gives N=" << std::setprecision(9)
      << Dmat.cwiseProduct(_dftAOoverlap.Matrix()).sum() << " electrons."
      << flush;

  XTP_LOG(Log::error, *_pLog) << TimeStamp() << " STARTING SCF cycle" << flush;
  XTP_LOG(Log::error, *_pLog)
      << " ----------------------------------------------"
         "----------------------------"
      << flush;

  for (Index this_iter = 0; this_iter < _max_iter; this_iter++) {
    XTP_LOG(Log::error, *_pLog) << flush;
    XTP_LOG(Log::error, *_pLog) << TimeStamp() << " Iteration " << this_iter + 1
                                << " of " << _max_iter << flush;

    Mat_p_Energy e_vxc = vxcpotential.IntegrateVXC(Dmat);
    XTP_LOG(Log::info, *_pLog)
        << TimeStamp() << " Filled DFT Vxc matrix " << flush;

    Mat_p_Energy ERIs = CalculateERIs(Dmat);
    Eigen::MatrixXd H = H0.matrix() + ERIs.matrix() + e_vxc.matrix();
    double Eone = Dmat.cwiseProduct(H0.matrix()).sum();
    double Etwo = 0.5 * ERIs.energy() + e_vxc.energy();
    double exx = 0.0;
    if (_ScaHFX > 0) {
      Mat_p_Energy EXXs = CalcEXXs(MOs.eigenvectors(), Dmat);
      XTP_LOG(Log::info, *_pLog)
          << TimeStamp() << " Filled DFT Electron exchange matrix" << flush;
      H -= 0.5 * _ScaHFX * EXXs.matrix();
      exx = -_ScaHFX / 4 * EXXs.energy();
    }
    Etwo += exx;
    double totenergy = Eone + H0.energy() + Etwo;
    XTP_LOG(Log::info, *_pLog) << TimeStamp() << " Single particle energy "
                               << std::setprecision(12) << Eone << flush;
    XTP_LOG(Log::info, *_pLog) << TimeStamp() << " Two particle energy "
                               << std::setprecision(12) << Etwo << flush;
    XTP_LOG(Log::info, *_pLog)
        << TimeStamp() << std::setprecision(12) << " Local Exc contribution "
        << e_vxc.energy() << flush;
    if (_ScaHFX > 0) {
      XTP_LOG(Log::info, *_pLog)
          << TimeStamp() << std::setprecision(12)
          << " Non local Ex contribution " << exx << flush;
    }
    XTP_LOG(Log::error, *_pLog) << TimeStamp() << " Total Energy "
                                << std::setprecision(12) << totenergy << flush;

    Dmat = _conv_accelerator.Iterate(Dmat, H, MOs, totenergy);

    PrintMOs(MOs.eigenvalues(), Log::info);

    XTP_LOG(Log::info, *_pLog) << "\t\tGAP "
                               << MOs.eigenvalues()(_numofelectrons / 2) -
                                      MOs.eigenvalues()(_numofelectrons / 2 - 1)
                               << flush;

    if (_conv_accelerator.isConverged()) {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Total Energy has converged to "
          << std::setprecision(9) << _conv_accelerator.getDeltaE()
          << "[Ha] after " << this_iter + 1
          << " iterations. DIIS error is converged up to "
          << _conv_accelerator.getDIIsError() << flush;
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Final Single Point Energy "
          << std::setprecision(12) << totenergy << " Ha" << flush;
      XTP_LOG(Log::error, *_pLog) << TimeStamp() << std::setprecision(12)
                                  << " Final Local Exc contribution "
                                  << e_vxc.energy() << " Ha" << flush;
      if (_ScaHFX > 0) {
        XTP_LOG(Log::error, *_pLog)
            << TimeStamp() << std::setprecision(12)
            << " Final Non Local Ex contribution " << exx << " Ha" << flush;
      }

      Mat_p_Energy EXXs = CalcEXXs(MOs.eigenvectors(), Dmat);
      exx = -1.0 / 4.0 * EXXs.energy();
      XTP_LOG(Log::error, *_pLog) << TimeStamp() << std::setprecision(12)
                                  << " EXX energy " << exx << " Ha" << flush;

      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << std::setprecision(12) << " Final EXX Total energy "
          << totenergy - e_vxc.energy() + (1.0 - _ScaHFX) * exx << " Ha"
          << flush;

      PrintMOs(MOs.eigenvalues(), Log::error);
      orb.setQMEnergy(totenergy);
      orb.MOs() = MOs;
      CalcElDipole(orb);
      break;
    } else if (this_iter == _max_iter - 1) {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " DFT calculation has not converged after "
          << _max_iter
          << " iterations. Use more iterations or another convergence "
             "acceleration scheme."
          << std::flush;
      return false;
    }
  }
  return true;
}

Mat_p_Energy DFTEngine::SetupH0(const QMMolecule& mol) const {

  AOKinetic dftAOkinetic;

  dftAOkinetic.Fill(_dftbasis);
  XTP_LOG(Log::info, *_pLog)
      << TimeStamp() << " Filled DFT Kinetic energy matrix ." << flush;

  AOMultipole dftAOESP;
  dftAOESP.FillPotential(_dftbasis, mol);
  XTP_LOG(Log::info, *_pLog)
      << TimeStamp() << " Filled DFT nuclear potential matrix." << flush;

  Eigen::MatrixXd H0 = dftAOkinetic.Matrix() + dftAOESP.Matrix();
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Constructed independent particle hamiltonian "
      << flush;
  double E0 = NuclearRepulsion(mol);
  XTP_LOG(Log::error, *_pLog) << TimeStamp() << " Nuclear Repulsion Energy is "
                              << std::setprecision(9) << E0 << flush;

  if (_with_ecp) {
    AOECP dftAOECP;
    dftAOECP.FillPotential(_dftbasis, _ecp);
    H0 += dftAOECP.Matrix();
    XTP_LOG(Log::info, *_pLog)
        << TimeStamp() << " Filled DFT ECP matrix" << flush;
  }

  if (_addexternalsites) {
    XTP_LOG(Log::error, *_pLog) << TimeStamp() << " " << _externalsites->size()
                                << " External sites" << flush;
    if (_externalsites->size() < 200) {
      XTP_LOG(Log::error, *_pLog)
          << " Name      Coordinates[a0]     charge[e]         dipole[e*a0]    "
             "              quadrupole[e*a0^2]         "
          << flush;

      for (const std::unique_ptr<StaticSite>& site : *_externalsites) {
        std::string output =
            (boost::format("  %1$s"
                           "   %2$+1.4f %3$+1.4f %4$+1.4f"
                           "   %5$+1.4f") %
             site->getElement() % site->getPos()[0] % site->getPos()[1] %
             site->getPos()[2] % site->getCharge())
                .str();
        const Eigen::Vector3d& dipole = site->getDipole();
        output += (boost::format("   %1$+1.4f %2$+1.4f %3$+1.4f") % dipole[0] %
                   dipole[1] % dipole[2])
                      .str();
        if (site->getRank() > 1) {
          Eigen::VectorXd quadrupole = site->Q().tail<5>();
          output += (boost::format(
                         "   %1$+1.4f %2$+1.4f %3$+1.4f %4$+1.4f %5$+1.4f") %
                     quadrupole[0] % quadrupole[1] % quadrupole[2] %
                     quadrupole[3] % quadrupole[4])
                        .str();
        }
        XTP_LOG(Log::error, *_pLog) << output << flush;
      }
    }

    Mat_p_Energy ext_multipoles =
        IntegrateExternalMultipoles(mol, *_externalsites);
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Nuclei-external site interaction energy "
        << std::setprecision(9) << ext_multipoles.energy() << flush;
    E0 += ext_multipoles.energy();
    H0 += ext_multipoles.matrix();
  }

  if (_integrate_ext_density) {
    Orbitals extdensity;
    extdensity.ReadFromCpt(_orbfilename);
    Mat_p_Energy extdensity_result = IntegrateExternalDensity(mol, extdensity);
    E0 += extdensity_result.energy();
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Nuclei- external density interaction energy "
        << std::setprecision(9) << extdensity_result.energy() << flush;
    H0 += extdensity_result.matrix();
  }

  if (_integrate_ext_field) {

    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Integrating external electric field with F[Hrt]="
        << _extfield.transpose() << flush;
    H0 += IntegrateExternalField(mol);
  }

  return Mat_p_Energy(E0, H0);
}

void DFTEngine::SetupInvariantMatrices() {

  _dftAOoverlap.Fill(_dftbasis);

  XTP_LOG(Log::info, *_pLog)
      << TimeStamp() << " Filled DFT Overlap matrix." << flush;

  _conv_opt.numberofelectrons = _numofelectrons;
  _conv_accelerator.Configure(_conv_opt);
  _conv_accelerator.setLogger(_pLog);
  _conv_accelerator.setOverlap(_dftAOoverlap, 1e-8);
  _conv_accelerator.PrintConfigOptions();

  if (_four_center_method == "RI") {
    // prepare invariant part of electron repulsion integrals
    _ERIs.Initialize(_dftbasis, _auxbasis);
    XTP_LOG(Log::info, *_pLog)
        << TimeStamp() << " Inverted AUX Coulomb matrix, removed "
        << _ERIs.Removedfunctions() << " functions from aux basis" << flush;
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp()
        << " Setup invariant parts of Electron Repulsion integrals " << flush;
  } else {

    if (_four_center_method == "cache") {

      XTP_LOG(Log::info, *_pLog)
          << TimeStamp() << " Calculating 4c integrals. " << flush;
      _ERIs.Initialize_4c_small_molecule(_dftbasis);
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Calculated 4c integrals. " << flush;
    }

    if (_with_screening && _four_center_method == "direct") {
      XTP_LOG(Log::info, *_pLog)
          << TimeStamp() << " Calculating 4c diagonals. " << flush;
      _ERIs.Initialize_4c_screening(_dftbasis, _screening_eps);
      XTP_LOG(Log::info, *_pLog)
          << TimeStamp() << " Calculated 4c diagonals. " << flush;
    }
  }

  return;
}

Eigen::MatrixXd DFTEngine::RunAtomicDFT_unrestricted(
    const QMAtom& uniqueAtom) const {
  bool with_ecp = _with_ecp;
  if (uniqueAtom.getElement() == "H" || uniqueAtom.getElement() == "He") {
    with_ecp = false;
  }

  QMMolecule atom = QMMolecule("individual_atom", 0);
  atom.push_back(uniqueAtom);

  BasisSet basisset;
  basisset.Load(_dftbasis_name);
  AOBasis dftbasis;
  dftbasis.Fill(basisset, atom);
  Vxc_Grid grid;
  grid.GridSetup(_grid_name, atom, dftbasis);
  Vxc_Potential<Vxc_Grid> gridIntegration(grid);
  gridIntegration.setXCfunctional(_xc_functional_name);

  ECPAOBasis ecp;
  if (with_ecp) {
    ECPBasisSet ecps;
    ecps.Load(_ecp_name);
    ecp.Fill(ecps, atom);
  }

  Index numofelectrons = uniqueAtom.getNuccharge();
  Index alpha_e = 0;
  Index beta_e = 0;

  if ((numofelectrons % 2) != 0) {
    alpha_e = numofelectrons / 2 + numofelectrons % 2;
    beta_e = numofelectrons / 2;
  } else {
    alpha_e = numofelectrons / 2;
    beta_e = alpha_e;
  }

  AOOverlap dftAOoverlap;
  AOKinetic dftAOkinetic;
  AOMultipole dftAOESP;
  AOECP dftAOECP;
  ERIs ERIs_atom;

  // DFT AOOverlap matrix

  dftAOoverlap.Fill(dftbasis);
  dftAOkinetic.Fill(dftbasis);

  dftAOESP.FillPotential(dftbasis, atom);
  ERIs_atom.Initialize_4c_small_molecule(dftbasis);

  ConvergenceAcc Convergence_alpha;
  ConvergenceAcc Convergence_beta;
  ConvergenceAcc::options opt_alpha = _conv_opt;
  opt_alpha.mode = ConvergenceAcc::KSmode::open;
  opt_alpha.histlength = 20;
  opt_alpha.levelshift = 0.1;
  opt_alpha.levelshiftend = 0.0;
  opt_alpha.usediis = true;
  opt_alpha.adiis_start = 0.0;
  opt_alpha.diis_start = 0.0;
  opt_alpha.numberofelectrons = alpha_e;

  ConvergenceAcc::options opt_beta = opt_alpha;
  opt_beta.numberofelectrons = beta_e;

  Logger log;
  Convergence_alpha.Configure(opt_alpha);
  Convergence_alpha.setLogger(&log);
  Convergence_alpha.setOverlap(dftAOoverlap, 1e-8);
  Convergence_beta.Configure(opt_beta);
  Convergence_beta.setLogger(&log);
  Convergence_beta.setOverlap(dftAOoverlap, 1e-8);
  /**** Construct initial density  ****/

  Eigen::MatrixXd H0 = dftAOkinetic.Matrix() + dftAOESP.Matrix();
  if (with_ecp) {
    dftAOECP.FillPotential(dftbasis, ecp);
    H0 += dftAOECP.Matrix();
  }
  tools::EigenSystem MOs_alpha = Convergence_alpha.SolveFockmatrix(H0);

  Eigen::MatrixXd dftAOdmat_alpha = Convergence_alpha.DensityMatrix(MOs_alpha);
  if (uniqueAtom.getElement() == "H") {
    return dftAOdmat_alpha;
  }
  tools::EigenSystem MOs_beta = Convergence_beta.SolveFockmatrix(H0);
  Eigen::MatrixXd dftAOdmat_beta = Convergence_beta.DensityMatrix(MOs_beta);

  Index maxiter = 80;
  for (Index this_iter = 0; this_iter < maxiter; this_iter++) {
    Mat_p_Energy ERIs = ERIs_atom.CalculateERIs_4c_small_molecule(
        dftAOdmat_alpha + dftAOdmat_beta);
    double E_two_alpha = ERIs.matrix().cwiseProduct(dftAOdmat_alpha).sum();
    double E_two_beta = ERIs.matrix().cwiseProduct(dftAOdmat_beta).sum();
    Eigen::MatrixXd H_alpha = H0 + ERIs.matrix();
    Eigen::MatrixXd H_beta = H0 + ERIs.matrix();

    Mat_p_Energy e_vxc = gridIntegration.IntegrateVXC(dftAOdmat_alpha);
    Eigen::MatrixXd AOVxc_alpha = e_vxc.matrix();
    double E_vxc_alpha = e_vxc.energy();
    H_alpha += AOVxc_alpha;
    E_two_alpha += E_vxc_alpha;

    e_vxc = gridIntegration.IntegrateVXC(dftAOdmat_beta);
    Eigen::MatrixXd AOVxc_beta = e_vxc.matrix();
    double E_vxc_beta = e_vxc.energy();
    H_beta += AOVxc_beta;
    E_two_beta += E_vxc_beta;

    if (_ScaHFX > 0) {
      Mat_p_Energy EXXs =
          ERIs_atom.CalculateEXX_4c_small_molecule(dftAOdmat_alpha);
      double E_exx_alpha =
          -0.5 * _ScaHFX * EXXs.matrix().cwiseProduct(dftAOdmat_alpha).sum();
      H_alpha -= _ScaHFX * EXXs.matrix();
      E_two_alpha += E_exx_alpha;
      ERIs_atom.CalculateEXX_4c_small_molecule(dftAOdmat_beta);
      double E_exx_beta =
          -0.5 * _ScaHFX * EXXs.matrix().cwiseProduct(dftAOdmat_beta).sum();
      H_beta -= _ScaHFX * EXXs.matrix();
      E_two_beta += E_exx_beta;
    }

    double E_one_alpha = dftAOdmat_alpha.cwiseProduct(H0).sum();
    double E_one_beta = dftAOdmat_beta.cwiseProduct(H0).sum();
    double E_alpha = E_one_alpha + E_two_alpha;
    double E_beta = E_one_beta + E_two_beta;
    double totenergy = E_alpha + E_beta;
    // evolve alpha
    dftAOdmat_alpha =
        Convergence_alpha.Iterate(dftAOdmat_alpha, H_alpha, MOs_alpha, E_alpha);
    // evolve beta
    dftAOdmat_beta =
        Convergence_beta.Iterate(dftAOdmat_beta, H_beta, MOs_beta, E_beta);

    XTP_LOG(Log::debug, *_pLog)
        << TimeStamp() << " Iter " << this_iter << " of " << maxiter << " Etot "
        << totenergy << " diise_a " << Convergence_alpha.getDIIsError()
        << " diise_b " << Convergence_beta.getDIIsError() << "\n\t\t a_gap "
        << MOs_alpha.eigenvalues()(alpha_e) -
               MOs_alpha.eigenvalues()(alpha_e - 1)
        << " b_gap "
        << MOs_beta.eigenvalues()(beta_e) - MOs_beta.eigenvalues()(beta_e - 1)
        << " Nalpha="
        << dftAOoverlap.Matrix().cwiseProduct(dftAOdmat_alpha).sum()
        << " Nbeta=" << dftAOoverlap.Matrix().cwiseProduct(dftAOdmat_beta).sum()
        << flush;

    bool converged =
        Convergence_alpha.isConverged() && Convergence_beta.isConverged();
    if (converged || this_iter == maxiter - 1) {

      if (converged) {
        XTP_LOG(Log::info, *_pLog) << TimeStamp() << " Converged after "
                                   << this_iter + 1 << " iterations" << flush;
      } else {
        XTP_LOG(Log::info, *_pLog)
            << TimeStamp() << " Not converged after " << this_iter + 1
            << " iterations. Unconverged density.\n\t\t\t"
            << " DIIsError_alpha=" << Convergence_alpha.getDIIsError()
            << " DIIsError_beta=" << Convergence_beta.getDIIsError() << flush;
      }
      break;
    }
  }
  Eigen::MatrixXd avgmatrix =
      SphericalAverageShells(dftAOdmat_alpha + dftAOdmat_beta, dftbasis);
  XTP_LOG(Log::info, *_pLog)
      << TimeStamp() << " Atomic density Matrix for " << uniqueAtom.getElement()
      << " gives N=" << std::setprecision(9)
      << avgmatrix.cwiseProduct(dftAOoverlap.Matrix()).sum() << " electrons."
      << flush;
  return avgmatrix;
}

Eigen::MatrixXd DFTEngine::AtomicGuess(const QMMolecule& mol) const {

  std::vector<std::string> elements = mol.FindUniqueElements();
  XTP_LOG(Log::info, *_pLog) << TimeStamp() << " Scanning molecule of size "
                             << mol.size() << " for unique elements" << flush;
  QMMolecule uniqueelements = QMMolecule("uniqueelements", 0);
  for (auto element : elements) {
    uniqueelements.push_back(QMAtom(0, element, Eigen::Vector3d::Zero()));
  }

  XTP_LOG(Log::info, *_pLog) << TimeStamp() << " " << uniqueelements.size()
                             << " unique elements found" << flush;
  std::vector<Eigen::MatrixXd> uniqueatom_guesses;
  for (QMAtom& unique_atom : uniqueelements) {
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Calculating atom density for "
        << unique_atom.getElement() << flush;
    Eigen::MatrixXd dmat_unrestricted = RunAtomicDFT_unrestricted(unique_atom);
    uniqueatom_guesses.push_back(dmat_unrestricted);
  }

  Eigen::MatrixXd guess =
      Eigen::MatrixXd::Zero(_dftbasis.AOBasisSize(), _dftbasis.AOBasisSize());
  Index start = 0;
  for (const QMAtom& atom : mol) {
    Index index = 0;
    for (; index < uniqueelements.size(); index++) {
      if (atom.getElement() == uniqueelements[index].getElement()) {
        break;
      }
    }
    Eigen::MatrixXd& dmat_unrestricted = uniqueatom_guesses[index];
    guess.block(start, start, dmat_unrestricted.rows(),
                dmat_unrestricted.cols()) = dmat_unrestricted;
    start += dmat_unrestricted.rows();
  }

  return guess;
}

void DFTEngine::ConfigOrbfile(Orbitals& orb) {
  if (_with_guess) {

    if (orb.hasDFTbasisName()) {
      if (orb.getDFTbasisName() != _dftbasis_name) {
        throw runtime_error(
            (boost::format("Basisset Name in guess orb file "
                           "and in dftengine option file differ %1% vs %2%") %
             orb.getDFTbasisName() % _dftbasis_name)
                .str());
      }
    } else {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp()
          << " WARNING: "
             "Orbital file has no basisset information,"
             "using it as a guess might work or not for calculation with "
          << _dftbasis_name << flush;
    }
  }
  orb.setDFTbasisName(_dftbasis_name);
  orb.setBasisSetSize(_dftbasis.AOBasisSize());
  orb.setXCFunctionalName(_xc_functional_name);
  orb.setScaHFX(_ScaHFX);
  if (_with_ecp) {
    orb.setECPName(_ecp_name);
  }
  if (_four_center_method == "RI") {
    orb.setAuxbasisName(_auxbasis_name);
  }

  if (_with_guess) {
    if (orb.hasECPName() || _with_ecp) {
      if (orb.getECPName() != _ecp_name) {
        throw runtime_error(
            (boost::format("ECPs in orb file: %1% and options %2% differ") %
             orb.getECPName() % _ecp_name)
                .str());
      }
    }
    if (orb.getNumberOfAlphaElectrons() != _numofelectrons / 2) {
      throw runtime_error(
          (boost::format("Number of electron in guess orb file: %1% and in "
                         "dftengine: %2% differ.") %
           orb.getNumberOfAlphaElectrons() % (_numofelectrons / 2))
              .str());
    }
    if (orb.getBasisSetSize() != _dftbasis.AOBasisSize()) {
      throw runtime_error((boost::format("Number of levels in guess orb file: "
                                         "%1% and in dftengine: %2% differ.") %
                           orb.getBasisSetSize() % _dftbasis.AOBasisSize())
                              .str());
    }
  } else {
    orb.setNumberOfAlphaElectrons(_numofelectrons / 2);
    orb.setNumberOfOccupiedLevels(_numofelectrons / 2);
  }
  return;
}

void DFTEngine::Prepare(QMMolecule& mol) {
  XTP_LOG(Log::error, *_pLog) << TimeStamp() << " Using "
                              << OPENMP::getMaxThreads() << " threads" << flush;

  if (XTP_HAS_MKL_OVERLOAD()) {
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Using MKL overload for Eigen " << flush;
  } else {
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp()
        << " Using native Eigen implementation, no BLAS overload " << flush;
  }

  XTP_LOG(Log::error, *_pLog) << " Molecule Coordinates [A] " << flush;
  for (const QMAtom& atom : mol) {
    const Eigen::Vector3d pos = atom.getPos() * tools::conv::bohr2ang;
    std::string output = (boost::format("  %1$s"
                                        "   %2$+1.4f %3$+1.4f %4$+1.4f") %
                          atom.getElement() % pos[0] % pos[1] % pos[2])
                             .str();

    XTP_LOG(Log::error, *_pLog) << output << flush;
  }
  BasisSet dftbasisset;
  dftbasisset.Load(_dftbasis_name);

  _dftbasis.Fill(dftbasisset, mol);
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Loaded DFT Basis Set " << _dftbasis_name << " with "
      << _dftbasis.AOBasisSize() << " functions" << flush;

  if (_four_center_method == "RI") {
    BasisSet auxbasisset;
    auxbasisset.Load(_auxbasis_name);
    _auxbasis.Fill(auxbasisset, mol);
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Loaded AUX Basis Set " << _auxbasis_name << " with "
        << _auxbasis.AOBasisSize() << " functions" << flush;
  }
  if (_with_ecp) {
    ECPBasisSet ecpbasisset;
    ecpbasisset.Load(_ecp_name);
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Loaded ECP library " << _ecp_name << flush;

    std::vector<std::string> results = _ecp.Fill(ecpbasisset, mol);
    XTP_LOG(Log::info, *_pLog) << TimeStamp() << " Filled ECP Basis of size "
                               << _ecp.ECPAOBasisSize() << flush;
    if (results.size() > 0) {
      std::string message = "";
      for (const std::string& element : results) {
        message += " " + element;
      }
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Found no ECPs for elements" << message << flush;
    }
  }

  for (const QMAtom& atom : mol) {
    _numofelectrons += atom.getNuccharge();
  }

  // here number of electrons is actually the total number, everywhere else in
  // votca it is just alpha_electrons
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Total number of electrons: " << _numofelectrons
      << flush;

  SetupInvariantMatrices();
  return;
}

Vxc_Potential<Vxc_Grid> DFTEngine::SetupVxc(const QMMolecule& mol) {
  _ScaHFX = Vxc_Potential<Vxc_Grid>::getExactExchange(_xc_functional_name);
  if (_ScaHFX > 0) {
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Using hybrid functional with alpha=" << _ScaHFX
        << flush;
  }
  Vxc_Grid grid;
  grid.GridSetup(_grid_name, mol, _dftbasis);
  Vxc_Potential<Vxc_Grid> vxc(grid);
  vxc.setXCfunctional(_xc_functional_name);
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Setup numerical integration grid " << _grid_name
      << " for vxc functional " << _xc_functional_name << flush;
  XTP_LOG(Log::info, *_pLog)
      << "\t\t "
      << " with " << grid.getGridSize() << " points"
      << " divided into " << grid.getBoxesSize() << " boxes" << flush;
  return vxc;
}

double DFTEngine::NuclearRepulsion(const QMMolecule& mol) const {
  double E_nucnuc = 0.0;

  for (Index i = 0; i < mol.size(); i++) {
    const Eigen::Vector3d& r1 = mol[i].getPos();
    double charge1 = double(mol[i].getNuccharge());
    for (Index j = 0; j < i; j++) {
      const Eigen::Vector3d& r2 = mol[j].getPos();
      double charge2 = double(mol[j].getNuccharge());
      E_nucnuc += charge1 * charge2 / (r1 - r2).norm();
    }
  }
  return E_nucnuc;
}

// average atom densities matrices, for SP and other combined shells average
// each subshell separately.
Eigen::MatrixXd DFTEngine::SphericalAverageShells(
    const Eigen::MatrixXd& dmat, const AOBasis& dftbasis) const {
  Eigen::MatrixXd avdmat = Eigen::MatrixXd::Zero(dmat.rows(), dmat.cols());
  Index start = 0.0;
  std::vector<Index> starts;
  std::vector<Index> ends;
  for (const AOShell& shell : dftbasis) {
    Index end = shell.getNumFunc() + start;
    starts.push_back(start);
    ends.push_back(end);
    start = end;
  }
  for (Index k = 0; k < Index(starts.size()); k++) {
    Index s1 = starts[k];
    Index e1 = ends[k];
    Index len1 = e1 - s1;
    for (Index l = 0; l < Index(starts.size()); l++) {
      Index s2 = starts[l];
      Index e2 = ends[l];
      Index len2 = e2 - s2;
      double diag = 0.0;
      double offdiag = 0.0;
      for (Index i = 0; i < len1; ++i) {
        for (Index j = 0; j < len2; ++j) {
          if (i == j) {
            diag += dmat(s1 + i, s2 + j);
          } else {
            offdiag += dmat(s1 + i, s2 + j);
          }
        }
      }
      if (len1 == len2) {
        diag = diag / double(len1);
        offdiag = offdiag / double(len1 * (len1 - 1));
      } else {
        double avg = (diag + offdiag) / double(len1 * len2);
        diag = avg;
        offdiag = avg;
      }
      for (Index i = 0; i < len1; ++i) {
        for (Index j = 0; j < len2; ++j) {
          if (i == j) {
            avdmat(s1 + i, s2 + j) = diag;
          } else {
            avdmat(s1 + i, s2 + j) = offdiag;
          }
        }
      }
    }
  }
  return avdmat;
}

double DFTEngine::ExternalRepulsion(
    const QMMolecule& mol,
    const std::vector<std::unique_ptr<StaticSite> >& multipoles) const {

  if (multipoles.size() == 0) {
    return 0;
  }

  double E_ext = 0;
  eeInteractor interactor;
  for (const QMAtom& atom : mol) {
    StaticSite nucleus = StaticSite(atom, double(atom.getNuccharge()));
    for (const std::unique_ptr<StaticSite>& site : *_externalsites) {
      if ((site->getPos() - nucleus.getPos()).norm() < 1e-7) {
        XTP_LOG(Log::error, *_pLog) << TimeStamp()
                                    << " External site sits on nucleus, "
                                       "interaction between them is ignored."
                                    << flush;
        continue;
      }
      E_ext += interactor.CalcStaticEnergy_site(*site, nucleus);
    }
  }
  return E_ext;
}

Eigen::MatrixXd DFTEngine::IntegrateExternalField(const QMMolecule& mol) const {

  AODipole dipole;
  dipole.setCenter(mol.getPos());
  dipole.Fill(_dftbasis);
  Eigen::MatrixXd result =
      Eigen::MatrixXd::Zero(dipole.Dimension(), dipole.Dimension());
  for (Index i = 0; i < 3; i++) {
    result -= dipole.Matrix()[i] * _extfield[i];
  }
  return result;
}

Mat_p_Energy DFTEngine::IntegrateExternalMultipoles(
    const QMMolecule& mol,
    const std::vector<std::unique_ptr<StaticSite> >& multipoles) const {

  Mat_p_Energy result(_dftbasis.AOBasisSize(), _dftbasis.AOBasisSize());
  AOMultipole dftAOESP;

  dftAOESP.FillPotential(_dftbasis, multipoles);
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Filled DFT external multipole potential matrix"
      << flush;
  result.matrix() = dftAOESP.Matrix();
  result.energy() = ExternalRepulsion(mol, multipoles);

  return result;
}

Mat_p_Energy DFTEngine::IntegrateExternalDensity(
    const QMMolecule& mol, const Orbitals& extdensity) const {
  BasisSet basis;
  basis.Load(extdensity.getDFTbasisName());
  AOBasis aobasis;
  aobasis.Fill(basis, extdensity.QMAtoms());
  Vxc_Grid grid;
  grid.GridSetup(_gridquality, extdensity.QMAtoms(), aobasis);
  DensityIntegration<Vxc_Grid> numint(grid);
  Eigen::MatrixXd dmat = extdensity.DensityMatrixFull(_state);

  numint.IntegrateDensity(dmat);
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Calculated external density" << flush;
  Eigen::MatrixXd e_contrib = numint.IntegratePotential(_dftbasis);
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Calculated potential from electron density" << flush;
  AOMultipole esp;
  esp.FillPotential(_dftbasis, extdensity.QMAtoms());

  double nuc_energy = 0.0;
  for (const QMAtom& atom : mol) {
    nuc_energy +=
        numint.IntegratePotential(atom.getPos()) * double(atom.getNuccharge());
    for (const QMAtom& extatom : extdensity.QMAtoms()) {
      const double dist = (atom.getPos() - extatom.getPos()).norm();
      nuc_energy +=
          double(atom.getNuccharge()) * double(extatom.getNuccharge()) / dist;
    }
  }
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Calculated potential from nuclei" << flush;
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Electrostatic: " << nuc_energy << flush;
  return Mat_p_Energy(nuc_energy, e_contrib + esp.Matrix());
}

Mat_p_Energy DFTEngine::CalculateERIs(const Eigen::MatrixXd& DMAT) const {
  if (_four_center_method == "RI") {
    return _ERIs.CalculateERIs(DMAT);
  } else if (_four_center_method == "cache") {
    return _ERIs.CalculateERIs_4c_small_molecule(DMAT);
  } else if (_four_center_method == "direct") {
    return _ERIs.CalculateERIs_4c_direct(_dftbasis, DMAT);
  } else {
    throw std::runtime_error("ERI method not known.");
  }
}

Eigen::MatrixXd DFTEngine::OrthogonalizeGuess(
    const Eigen::MatrixXd& GuessMOs) const {
  Eigen::MatrixXd nonortho =
      GuessMOs.transpose() * _dftAOoverlap.Matrix() * GuessMOs;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(nonortho);
  Eigen::MatrixXd result = GuessMOs * es.operatorInverseSqrt();
  return result;
}

}  // namespace xtp
}  // namespace votca
