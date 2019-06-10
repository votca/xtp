/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <votca/ctp/logger.h>
#include <votca/tools/constants.h>
#include <votca/xtp/bse.h>
#include <votca/xtp/gw.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/numerical_integrations.h>
#include <votca/xtp/orbitals.h>

using boost::format;
using namespace boost::filesystem;
using std::flush;
namespace votca {
namespace xtp {

int GWBSE::CountCoreLevels() {
  int ignored_corelevels = 0;
  if (!_orbitals.hasECPName()) {
    BasisSet basis;
    basis.LoadPseudopotentialSet("corelevels");
    int coreElectrons = 0;
    for (const auto& atom : _orbitals.QMAtoms()) {
      coreElectrons += basis.getElement(atom->getType()).getNcore();
    }
    ignored_corelevels = coreElectrons / 2;
  }
  return ignored_corelevels;
}

void GWBSE::Initialize(tools::Property& options) {

  std::string key = Identify();

  // getting level ranges
  int rpamax = 0;
  int rpamin = 0;  // never changes
  int qpmin = 0;
  int qpmax = 0;
  int bse_vmin = 0;
  int bse_cmax = 0;

  int homo = _orbitals.getHomo();  // indexed from 0
  int num_of_levels = _orbitals.getBasisSetSize();
  int num_of_occlevels = _orbitals.getNumberOfAlphaElectrons();

  std::vector<std::string> range_choice = {"default", "factor", "explicit",
                                           "full"};
  std::string ranges =
      options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(
          key + ".ranges", range_choice);

  // now check validity, and get rpa, qp, and bse level ranges accordingly

  if (ranges == "factor") {

    double rpamaxfactor = options.get(key + ".rpamax").as<double>();
    rpamax = int(rpamaxfactor * num_of_levels) - 1;  // total number of levels

    double qpminfactor = options.get(key + ".qpmin").as<double>();
    qpmin = num_of_occlevels - int(qpminfactor * num_of_occlevels) - 1;

    double qpmaxfactor = options.get(key + ".qpmax").as<double>();
    qpmax = num_of_occlevels + int(qpmaxfactor * num_of_occlevels) - 1;

    double bseminfactor = options.get(key + ".bsemin").as<double>();
    bse_vmin = num_of_occlevels - int(bseminfactor * num_of_occlevels) - 1;

    double bsemaxfactor = options.get(key + ".bsemax").as<double>();
    bse_cmax = num_of_occlevels + int(bsemaxfactor * num_of_occlevels) - 1;

  } else if (ranges == "explicit") {
    // get explicit numbers
    rpamax = options.get(key + ".rpamax").as<int>();
    qpmin = options.get(key + ".qpmin").as<int>();
    qpmax = options.get(key + ".qpmax").as<int>();
    bse_vmin = options.get(key + ".bsemin").as<int>();
    bse_cmax = options.get(key + ".bsemax").as<int>();
  } else if (ranges == "default") {
    rpamax = num_of_levels - 1;
    qpmin = 0;
    qpmax = 2 * homo + 1;
    bse_vmin = 0;
    bse_cmax = 2 * homo + 1;
  } else if (ranges == "full") {
    rpamax = num_of_levels - 1;
    qpmin = 0;
    qpmax = num_of_levels - 1;
    bse_vmin = 0;
    bse_cmax = num_of_levels - 1;
  }
  std::string ignore_corelevels =
      options.ifExistsReturnElseReturnDefault<std::string>(
          key + ".ignore_corelevels", "none");

  if (ignore_corelevels == "RPA" || ignore_corelevels == "GW" ||
      ignore_corelevels == "BSE") {
    int ignored_corelevels = CountCoreLevels();
    if (ignore_corelevels == "RPA") {
      rpamin = ignored_corelevels;
    }
    if (ignore_corelevels == "GW" || ignore_corelevels == "RPA") {
      if (qpmin < ignored_corelevels) {
        qpmin = ignored_corelevels;
      }
    }
    if (ignore_corelevels == "GW" || ignore_corelevels == "RPA" ||
        ignore_corelevels == "BSE") {
      if (bse_vmin < ignored_corelevels) {
        bse_vmin = ignored_corelevels;
      }
    }

    CTP_LOG(ctp::logDEBUG, *_pLog)
        << ctp::TimeStamp() << " Ignoring " << ignored_corelevels
        << " core levels for " << ignore_corelevels << " and beyond." << flush;
  }

  // check maximum and minimum sizes
  if (rpamax > num_of_levels) rpamax = num_of_levels - 1;
  if (qpmax > num_of_levels) qpmax = num_of_levels - 1;
  if (bse_cmax > num_of_levels) bse_cmax = num_of_levels - 1;
  if (bse_vmin < 0) bse_vmin = 0;
  if (qpmin < 0) qpmin = 0;

  // some QP - BSE consistency checks are required
  if (bse_vmin < qpmin) qpmin = bse_vmin;
  if (bse_cmax > qpmax) qpmax = bse_cmax;

  _gwopt.homo = homo;
  _gwopt.qpmin = qpmin;
  _gwopt.qpmax = qpmax;
  _gwopt.rpamin = rpamin;
  _gwopt.rpamax = rpamax;

  _bseopt.vmin = bse_vmin;
  _bseopt.cmax = bse_cmax;
  _bseopt.homo = homo;
  _bseopt.qpmin = qpmin;
  _bseopt.rpamin = rpamin;
  _bseopt.rpamax = rpamax;

  _orbitals.setRPAindices(rpamin, rpamax);
  _orbitals.setGWindices(qpmin, qpmax);
  _orbitals.setBSEindices(bse_vmin, bse_cmax);

  int bse_vmax = homo;
  int bse_cmin = homo + 1;
  int bse_vtotal = bse_vmax - bse_vmin + 1;
  int bse_ctotal = bse_cmax - bse_cmin + 1;
  int bse_size = bse_vtotal * bse_ctotal;

  CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Set RPA level range ["
                                 << rpamin << ":" << rpamax << "]" << flush;
  CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Set GW  level range ["
                                 << qpmin << ":" << qpmax << "]" << flush;
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp() << " Set BSE level range occ[" << bse_vmin << ":"
      << bse_vmax << "]  virt[" << bse_cmin << ":" << bse_cmax << "]" << flush;
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp() << " BSE Hamiltonian has size " << bse_size << "x"
      << bse_size << flush;

  _gwopt.reset_3c = options.ifExistsReturnElseReturnDefault<int>(
      key + ".rebuild_threecenter_freq", _gwopt.reset_3c);

  _bseopt.nmax = options.ifExistsReturnElseReturnDefault<int>(key + ".exctotal",
                                                              _bseopt.nmax);
  if (_bseopt.nmax > bse_size || _bseopt.nmax < 0) _bseopt.nmax = bse_size;

  // eigensolver options
  if (options.exists(key + ".eigensolver")) {
    _bseopt.davidson = options.ifExistsReturnElseReturnDefault<bool>(
        key + ".eigensolver.dodavidson", _bseopt.davidson);

    if (_bseopt.davidson) {

      _bseopt.matrixfree = options.ifExistsReturnElseReturnDefault<bool>(
          key + ".eigensolver.domatrixfree", _bseopt.matrixfree);

      _bseopt.davidson_correction =
          options.ifExistsReturnElseReturnDefault<std::string>(
              key + ".eigensolver.davidson_correction",
              _bseopt.davidson_correction);

      _bseopt.davidson_ortho =
          options.ifExistsReturnElseReturnDefault<std::string>(
              key + ".eigensolver.davidson_ortho", _bseopt.davidson_ortho);

      _bseopt.davidson_tolerance =
          options.ifExistsReturnElseReturnDefault<std::string>(
              key + ".eigensolver.davidson_tolerance",
              _bseopt.davidson_tolerance);

      _bseopt.davidson_update =
          options.ifExistsReturnElseReturnDefault<std::string>(
              key + ".eigensolver.davidson_update", _bseopt.davidson_update);

      _bseopt.davidson_maxiter = options.ifExistsReturnElseReturnDefault<int>(
          key + ".eigensolver.davidson_maxiter", _bseopt.davidson_maxiter);

      std::vector<std::string> _dcorr = {"DPR", "OLSEN"};
      options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(
          key + ".eigensolver.davidson_correction", _dcorr);

      std::vector<std::string> _dortho = {"GS", "QR"};
      options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(
          key + ".eigensolver.davidson_ortho", _dortho);

      std::vector<std::string> _dtol = {"strict", "normal", "loose"};
      options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(
          key + ".eigensolver.davidson_tolerance", _dtol);

      std::vector<std::string> _dup = {"min", "safe", "max"};
      options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(
          key + ".eigensolver.davidson_update", _dup);

      // check size
      if (_bseopt.nmax > bse_size / 4) {
        CTP_LOG(ctp::logDEBUG, *_pLog)
            << ctp::TimeStamp()
            << " Warning : Too many eigenvalues required for Davidson. Default "
               "to Lapack diagonalization"
            << flush;
        _bseopt.davidson = false;
      }
    }
  }

  _fragA = options.ifExistsReturnElseReturnDefault<int>(key + ".fragment", -1);

  _bseopt.useTDA = options.ifExistsReturnElseReturnDefault<bool>(
      key + ".useTDA", _bseopt.useTDA);
  _orbitals.setTDAApprox(_bseopt.useTDA);
  if (!_bseopt.useTDA) {
    CTP_LOG(ctp::logDEBUG, *_pLog) << " BSE type: full" << flush;
  } else {
    CTP_LOG(ctp::logDEBUG, *_pLog) << " BSE type: TDA" << flush;
  }

  _openmp_threads =
      options.ifExistsReturnElseReturnDefault<int>(key + ".openmp", 0);

  if (options.exists(key + ".vxc")) {
    _doVxc =
        options.ifExistsReturnElseThrowRuntimeError<bool>(key + ".vxc.dovxc");
    if (_doVxc) {
      _functional = options.ifExistsReturnElseThrowRuntimeError<std::string>(
          key + ".vxc.functional");
      _grid = options.ifExistsReturnElseReturnDefault<std::string>(
          key + ".vxc.grid", "medium");
    }
  }

  if (options.exists(key + ".gwbasis")) {
    _auxbasis_name = options.ifExistsReturnElseThrowRuntimeError<std::string>(
        key + ".gwbasis");
    if (options.exists(key + ".auxbasis")) {
      throw std::runtime_error(
          std::string() +
          "Cannot use the options \"gwbasis\" and \"auxbasis\" "
          "simultaneously. " +
          " Use option \"auxbasis\" to specify the auxiliary basis instead.");
    } else {
      CTP_LOG(ctp::logDEBUG, *_pLog)
          << " Warning: The option \"gwbasis\" is outdated."
          << " Use option \"auxbasis\" to specify the auxiliary basis instead. "
          << flush;
    }
  } else {
    _auxbasis_name = options.ifExistsReturnElseThrowRuntimeError<std::string>(
        key + ".auxbasis");
  }

  _dftbasis_name = options.ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".dftbasis");
  if (_dftbasis_name != _orbitals.getDFTbasisName()) {
    throw std::runtime_error(
        "Name of the Basisset from .orb file: " + _orbitals.getDFTbasisName() +
        " and from GWBSE optionfile " + _dftbasis_name + " do not agree.");
  }

  std::vector<std::string> choices = {"G0W0", "evGW"};
  std::string mode =
      options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(
          key + ".mode", choices);
  if (mode == "G0W0") {
    _gwopt.gw_sc_max_iterations = 1;
  }
  CTP_LOG(ctp::logDEBUG, *_pLog) << " Running GW as: " << mode << flush;
  _gwopt.ScaHFX = _orbitals.getScaHFX();

  _gwopt.shift = options.ifExistsReturnElseReturnDefault<double>(key + ".shift",
                                                                 _gwopt.shift);
  _gwopt.g_sc_limit = options.ifExistsReturnElseReturnDefault<double>(
      key + ".g_sc_limit",
      _gwopt.g_sc_limit);  // convergence criteria for qp iteration [Hartree]]
  _gwopt.g_sc_max_iterations = options.ifExistsReturnElseReturnDefault<int>(
      key + ".g_sc_max_iterations",
      _gwopt.g_sc_max_iterations);  // convergence criteria for qp iteration
                                    // [Hartree]]

  _gwopt.gw_sc_max_iterations = options.ifExistsReturnElseReturnDefault<int>(
      key + ".gw_sc_max_iterations",
      _gwopt.gw_sc_max_iterations);  // convergence criteria for qp iteration
                                     // [Hartree]]

  _gwopt.gw_sc_limit = options.ifExistsReturnElseReturnDefault<double>(
      key + ".gw_sc_limit",
      _gwopt.gw_sc_limit);  // convergence criteria for shift it
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << " g_sc_limit [Hartree]: " << _gwopt.g_sc_limit << flush;
  if (_gwopt.gw_sc_max_iterations > 1) {
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << " gw_sc_limit [Hartree]: " << _gwopt.gw_sc_limit << flush;
  }
  _bseopt.min_print_weight = options.ifExistsReturnElseReturnDefault<double>(
      key + ".bse_print_weight", _bseopt.min_print_weight);
  // print exciton WF composition weight larger than minimum

  _gwopt.sigma_integration =
      options.ifExistsReturnElseReturnDefault<std::string>(
          key + ".sigma_integrator", _gwopt.sigma_integration);
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << " Sigma Integration: " << _gwopt.sigma_integration << flush;
  if (_gwopt.sigma_integration == "exact") {
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << " RPA size: " << (homo + 1 - rpamin) * (rpamax - homo) << flush;
  }
  
  _gwopt.sigma_eta = options.ifExistsReturnElseReturnDefault<double>(
      key + ".sigma_eta",
      CustomOpts::SigmaSpectralEta());
  if (_gwopt.sigma_integration == "exact") {
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << " eta: " << _gwopt.sigma_eta << flush;
  }
  
  // TODO: Allow string input
  std::vector<std::string> root_finder_choice =
      {"fixed point", "bisection", "regula falsi", "grid"};
  _gwopt.gw_sc_root_finder_method =
      options.ifExistsReturnElseReturnDefault<int>(
          key + ".gw_sc_root_finder_method", _gwopt.gw_sc_root_finder_method);
  if (_gwopt.gw_sc_root_finder_method < 0 || _gwopt.gw_sc_root_finder_method >= root_finder_choice.size()) {
    throw std::runtime_error(
        (boost::format("Root finder must be within [0, %d]") % (root_finder_choice.size() - 1)).str());
  }
  _gwopt.gw_sc_root_finder_range = options.ifExistsReturnElseReturnDefault<double>(
      key + ".gw_sc_root_finder_range",
      _gwopt.gw_sc_root_finder_range); 
  _gwopt.gw_sc_root_finder_steps = options.ifExistsReturnElseReturnDefault<int>(
      key + ".gw_sc_root_finder_steps",
      _gwopt.gw_sc_root_finder_steps);
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << " Root finder method: " << root_finder_choice[_gwopt.gw_sc_root_finder_method] << flush;
  if (_gwopt.gw_sc_root_finder_method == 3) {
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << " Root finder range: " << _gwopt.gw_sc_root_finder_range << flush;
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << " Root finder steps: " << _gwopt.gw_sc_root_finder_steps << flush;
    _gwopt.g_sc_max_iterations = 1; // Grid method does not support mutliple g iterations
  }

  // setting some defaults
  _do_bse_singlets = false;
  _do_bse_triplets = false;

  // possible tasks
  std::string tasks_string =
      options.ifExistsReturnElseThrowRuntimeError<std::string>(key + ".tasks");
  if (tasks_string.find("all") != std::string::npos) {
    _do_gw = true;
    _do_bse_singlets = true;
    _do_bse_triplets = true;
  }
  if (tasks_string.find("gw") != std::string::npos) _do_gw = true;
  if (tasks_string.find("singlets") != std::string::npos)
    _do_bse_singlets = true;
  if (tasks_string.find("triplets") != std::string::npos)
    _do_bse_triplets = true;
  // special construction for ibse mode

  std::string store_string =
      options.ifExistsReturnElseThrowRuntimeError<std::string>(key + ".store");
  if ((store_string.find("all") != std::string::npos) ||
      (store_string.find("") != std::string::npos)) {
    // store according to tasks choice
    if (_do_bse_singlets) _store_bse_singlets = true;
    if (_do_bse_triplets) _store_bse_triplets = true;
  }
  if (store_string.find("singlets") != std::string::npos)
    _store_bse_singlets = true;
  if (store_string.find("triplets") != std::string::npos)
    _store_bse_triplets = true;

  CTP_LOG(ctp::logDEBUG, *_pLog) << " Tasks: " << flush;
  if (_do_gw) {
    CTP_LOG(ctp::logDEBUG, *_pLog) << " GW " << flush;
  }
  if (_do_bse_singlets) {
    CTP_LOG(ctp::logDEBUG, *_pLog) << " singlets " << flush;
  }
  if (_do_bse_triplets) {
    CTP_LOG(ctp::logDEBUG, *_pLog) << " triplets " << flush;
  }
  CTP_LOG(ctp::logDEBUG, *_pLog) << " Store: " << flush;
  if (_do_gw) {
    CTP_LOG(ctp::logDEBUG, *_pLog) << " GW " << flush;
  }
  if (_store_bse_singlets) {
    CTP_LOG(ctp::logDEBUG, *_pLog) << " singlets " << flush;
  }
  if (_store_bse_triplets) {
    CTP_LOG(ctp::logDEBUG, *_pLog) << " triplets " << flush;
  }

  return;
}

void GWBSE::addoutput(tools::Property& summary) {

  const double hrt2ev = tools::conv::hrt2ev;
  tools::Property& gwbse_summary = summary.add("GWBSE", "");
  if (_do_gw) {
    gwbse_summary.setAttribute("units", "eV");
    gwbse_summary.setAttribute(
        "DFTEnergy",
        (format("%1$+1.6f ") % (_orbitals.getQMEnergy() * hrt2ev)).str());

    tools::Property& dft_summary = gwbse_summary.add("dft", "");
    dft_summary.setAttribute("HOMO", _gwopt.homo);
    dft_summary.setAttribute("LUMO", _gwopt.homo + 1);

    for (int state = 0; state < _gwopt.qpmax + 1 - _gwopt.qpmin; state++) {
      tools::Property& level_summary = dft_summary.add("level", "");
      level_summary.setAttribute("number", state + _gwopt.qpmin);
      level_summary.add("dft_energy",
                        (format("%1$+1.6f ") %
                         (_orbitals.QPpertEnergies().col(0)(state) * hrt2ev))
                            .str());
      level_summary.add("gw_energy",
                        (format("%1$+1.6f ") %
                         (_orbitals.QPpertEnergies().col(4)(state) * hrt2ev))
                            .str());

      level_summary.add(
          "qp_energy",
          (format("%1$+1.6f ") % (_orbitals.QPdiagEnergies()(state) * hrt2ev))
              .str());
    }
  }
  if (_do_bse_singlets) {
    tools::Property& singlet_summary = gwbse_summary.add("singlets", "");
    for (int state = 0; state < _bseopt.nmax; ++state) {
      tools::Property& level_summary = singlet_summary.add("level", "");
      level_summary.setAttribute("number", state + 1);
      level_summary.add("omega",
                        (format("%1$+1.6f ") %
                         (_orbitals.BSESingletEnergies()(state) * hrt2ev))
                            .str());
      if (_orbitals.hasTransitionDipoles()) {

        const tools::vec& dipoles = (_orbitals.TransitionDipoles())[state];
        double f =
            2 * dipoles * dipoles * _orbitals.BSESingletEnergies()(state) / 3.0;

        level_summary.add("f", (format("%1$+1.6f ") % f).str());
        tools::Property& dipol_summary = level_summary.add(
            "Trdipole", (format("%1$+1.4f %2$+1.4f %3$+1.4f") % dipoles.getX() %
                         dipoles.getY() % dipoles.getZ())
                            .str());
        dipol_summary.setAttribute("unit", "e*bohr");
        dipol_summary.setAttribute("gauge", "length");
      }
    }
  }
  if (_do_bse_triplets) {
    tools::Property& triplet_summary = gwbse_summary.add("triplets", "");
    for (int state = 0; state < _bseopt.nmax; ++state) {
      tools::Property& level_summary = triplet_summary.add("level", "");
      level_summary.setAttribute("number", state + 1);
      level_summary.add("omega",
                        (format("%1$+1.6f ") %
                         (_orbitals.BSETripletEnergies()(state) * hrt2ev))
                            .str());
    }
  }
  return;
}

/*
 *    Many-body Green's fuctions theory implementation
 *
 *  data required from orbitals file
 *  - atomic coordinates
 *  - DFT molecular orbitals (energies and coeffcients)
 *  - DFT exchange-correlation potential matrix in atomic orbitals
 *  - number of electrons, number of levels
 */

Eigen::MatrixXd GWBSE::CalculateVXC(const AOBasis& dftbasis) {

  Eigen::MatrixXd vxc_ao;
  if (_orbitals.hasAOVxc()) {
    if (_doVxc) {
      if (_orbitals.getQMpackage() == "xtp") {
        CTP_LOG(ctp::logDEBUG, *_pLog)
            << ctp::TimeStamp() << " Taking VXC from xtp DFT run." << flush;
      } else {
        CTP_LOG(ctp::logDEBUG, *_pLog)
            << ctp::TimeStamp()
            << " There is already a Vxc matrix loaded from DFT, did you maybe "
               "run a DFT code with outputVxc?\n I will take the external "
               "implementation"
            << flush;
      }
      _doVxc = false;
    }
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << ctp::TimeStamp() << " Loaded external Vxc matrix" << flush;
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
        << " Setup grid for integration with gridsize: " << _grid << " with "
        << numint.getGridSize() << " points, divided into "
        << numint.getBoxesSize() << " boxes" << flush;
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << ctp::TimeStamp() << " Integrating Vxc in VOTCA with functional "
        << _functional << flush;
    Eigen::MatrixXd DMAT = _orbitals.DensityMatrixGroundState();

    vxc_ao = numint.IntegrateVXC(DMAT);
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << ctp::TimeStamp() << " Calculated Vxc in VOTCA" << flush;

  } else {
    throw std::runtime_error(
        "So your DFT data contains no Vxc, if you want to proceed use the "
        "dovxc option.");
  }

  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp()
      << " Set hybrid exchange factor: " << _orbitals.getScaHFX() << flush;
  int qptotal = _gwopt.qpmax - _gwopt.qpmin + 1;
  int basissize = _orbitals.MOCoefficients().rows();
  Eigen::MatrixXd mos =
      _orbitals.MOCoefficients().block(0, _gwopt.qpmin, basissize, qptotal);

  Eigen::MatrixXd vxc = mos.transpose() * vxc_ao * mos;
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp()
      << " Calculated exchange-correlation expectation values " << flush;

  return vxc;
}

bool GWBSE::Evaluate() {

// set the parallelization
#ifdef _OPENMP
  if (_openmp_threads > 0) {
    omp_set_num_threads(_openmp_threads);
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << ctp::TimeStamp() << " Using " << omp_get_max_threads() << " threads"
        << flush;
  }
#endif

  if (tools::globals::VOTCA_MKL) {
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << ctp::TimeStamp() << " Using MKL overload for Eigen " << flush;
  } else {
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << ctp::TimeStamp()
        << " Using native Eigen implementation, no BLAS overload " << flush;
  }

  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp() << " Molecule Coordinates [A] " << flush;
  for (QMAtom* atom : _orbitals.QMAtoms()) {
    std::string output =
        (boost::format("  %1$s"
                       "   %2$+1.4f %3$+1.4f %4$+1.4f") %
         atom->getType() % (atom->getPos().getX() * tools::conv::bohr2ang) %
         (atom->getPos().getY() * tools::conv::bohr2ang) %
         (atom->getPos().getZ() * tools::conv::bohr2ang))
            .str();

    CTP_LOG(ctp::logDEBUG, *_pLog) << output << flush;
  }

  std::string dft_package = _orbitals.getQMpackage();
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp() << " DFT data was created by " << dft_package
      << flush;

  BasisSet dftbs;
  dftbs.LoadBasisSet(_dftbasis_name);

  CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Loaded DFT Basis Set "
                                 << _dftbasis_name << flush;

  // fill DFT AO basis by going through all atoms
  AOBasis dftbasis;
  dftbasis.AOBasisFill(dftbs, _orbitals.QMAtoms(), _fragA);
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp() << " Filled DFT Basis of size "
      << dftbasis.AOBasisSize() << flush;
  if (dftbasis.getAOBasisFragB() > 0 && dftbasis.getAOBasisFragA() > 0) {
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " FragmentA size "
                                   << dftbasis.getAOBasisFragA() << flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " FragmentB size "
                                   << dftbasis.getAOBasisFragB() << flush;
  }

  // load auxiliary basis set (element-wise information) from xml file
  BasisSet auxbs;
  auxbs.LoadBasisSet(_auxbasis_name);
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp() << " Loaded Auxbasis Set " << _auxbasis_name << flush;

  // fill auxiliary AO basis by going through all atoms
  AOBasis auxbasis;
  auxbasis.AOBasisFill(auxbs, _orbitals.QMAtoms());
  _orbitals.setAuxbasisName(_auxbasis_name);
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp() << " Filled Auxbasis of size "
      << auxbasis.AOBasisSize() << flush;

  TCMatrix_gwbse Mmn;
  // rpamin here, because RPA needs till rpamin
  Mmn.Initialize(auxbasis.AOBasisSize(), _gwopt.rpamin, _gwopt.qpmax,
                 _gwopt.rpamin, _gwopt.rpamax);
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp()
      << " Calculating Mmn_beta (3-center-repulsion x orbitals)  " << flush;
  Mmn.Fill(auxbasis, dftbasis, _orbitals.MOCoefficients());
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp() << " Removed " << Mmn.Removedfunctions()
      << " functions from Aux Coulomb matrix to avoid near linear dependencies"
      << flush;
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp()
      << " Calculated Mmn_beta (3-center-repulsion x orbitals)  " << flush;

  Eigen::MatrixXd Hqp;

  if (_do_gw) {
    Eigen::MatrixXd vxc = CalculateVXC(dftbasis);
    GW gw = GW(*_pLog, Mmn, vxc, _orbitals.MOEnergies());
    gw.configure(_gwopt);
    gw.CalculateGWPerturbation();

    // store perturbative QP energy data in orbitals object (DFT, S_x,S_c, V_xc,
    // E_qp)
    _orbitals.QPpertEnergies() = gw.getGWAResults();

    CTP_LOG(ctp::logDEBUG, *_pLog)
        << ctp::TimeStamp() << " Calculating offdiagonal part of Sigma  "
        << flush;
    gw.CalculateHQP();
    CTP_LOG(ctp::logDEBUG, *_pLog)
        << ctp::TimeStamp() << " Calculated offdiagonal part of Sigma  "
        << flush;
    Hqp = gw.getHQP();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es =
        gw.DiagonalizeQPHamiltonian();
    if (es.info() == Eigen::ComputationInfo::Success) {
      CTP_LOG(ctp::logDEBUG, *_pLog)
          << ctp::TimeStamp() << " Diagonalized QP Hamiltonian  " << flush;
    }

    _orbitals.QPdiagCoefficients() = es.eigenvectors();
    _orbitals.QPdiagEnergies() = es.eigenvalues();
  } else {
    if (_orbitals.hasQPdiag()) {
      const Eigen::MatrixXd& qpcoeff = _orbitals.QPdiagCoefficients();
      Hqp = qpcoeff * _orbitals.QPdiagEnergies().asDiagonal() *
            qpcoeff.transpose();
    } else {
      throw std::runtime_error("orb file has no QPcoefficients");
    }
  }

  // proceed only if BSE requested
  if (_do_bse_singlets || _do_bse_triplets) {
    BSE bse = BSE(_orbitals, *_pLog, Mmn, Hqp);
    bse.configure(_bseopt);

    if (_do_bse_triplets) {
      bse.Solve_triplets();
      CTP_LOG(ctp::logDEBUG, *_pLog)
          << ctp::TimeStamp() << " Solved BSE for triplets " << flush;
      bse.Analyze_triplets(dftbasis);
      if (!_store_bse_triplets) {
        bse.FreeTriplets();
      }
    }

    if (_do_bse_singlets) {
      bse.Solve_singlets();
      CTP_LOG(ctp::logDEBUG, *_pLog)
          << ctp::TimeStamp() << " Solved BSE for singlets " << flush;
      bse.Analyze_singlets(dftbasis);
      if (!_store_bse_singlets) {
        bse.FreeSinglets();
      }
    }
  }
  CTP_LOG(ctp::logDEBUG, *_pLog)
      << ctp::TimeStamp() << " GWBSE calculation finished " << flush;
  return true;
}
}  // namespace xtp
};  // namespace votca
