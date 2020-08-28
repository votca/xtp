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
#include <cstdio>
#include <iomanip>

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/xtp/basisset.h"
#include "votca/xtp/ecpaobasis.h"
#include "votca/xtp/orbitals.h"

// Local private VOTCA includes
#include "orca.h"

namespace votca {
namespace xtp {
using namespace std;

void Orca::Initialize(const tools::Property& options) {

  // good luck

  // Orca file names
  const std::string& fileName =
      options.ifExistsReturnElseReturnDefault<std::string>("job_name",
                                                           "system");

  _input_file_name = fileName + ".inp";
  _log_file_name = fileName + ".log";
  _shell_file_name = fileName + ".sh";
  _mo_file_name = fileName + ".gbw";

  ParseCommonOptions(options);
}

/* Custom basis sets are written on a per-element basis to
 * the system.bas/aux file(s), which are then included in the
 * Orca input file using GTOName = "system.bas/aux"
 */
void Orca::WriteBasisset(const QMMolecule& qmatoms, std::string& bs_name,
                         std::string& el_file_name) {

  std::vector<std::string> UniqueElements = qmatoms.FindUniqueElements();

  tools::Elements elementInfo;
  BasisSet bs;
  bs.Load(bs_name);
  XTP_LOG(Log::error, *_pLog) << "Loaded Basis Set " << bs_name << flush;
  ofstream el_file;

  el_file.open(el_file_name);
  el_file << "$DATA" << endl;

  for (const std::string& element_name : UniqueElements) {
    const Element& element = bs.getElement(element_name);
    el_file << elementInfo.getEleFull(element_name) << endl;
    for (const Shell& shell : element) {
      el_file << xtp::EnumToString(shell.getL()) << " " << shell.getSize()
              << endl;
      Index sh_idx = 0;
      for (const GaussianPrimitive& gaussian : shell) {
        sh_idx++;
        el_file << " " << sh_idx << " " << indent(gaussian.decay());
        el_file << " " << indent(gaussian.contraction());
        el_file << endl;
      }
    }
  }
  el_file << "STOP\n";
  el_file.close();

  return;
}

/* Coordinates are written in standard Element,x,y,z format to the
 * input file.
 */
void Orca::WriteCoordinates(std::ofstream& inp_file,
                            const QMMolecule& qmatoms) {

  for (const QMAtom& atom : qmatoms) {
    Eigen::Vector3d pos = atom.getPos() * tools::conv::bohr2ang;
    inp_file << setw(3) << atom.getElement() << setw(12)
             << setiosflags(ios::fixed) << setprecision(5) << pos.x()
             << setw(12) << setiosflags(ios::fixed) << setprecision(5)
             << pos.y() << setw(12) << setiosflags(ios::fixed)
             << setprecision(5) << pos.z() << endl;
  }
  inp_file << "* \n" << endl;
  return;
}

/* If custom ECPs are used, they need to be specified in the input file
 * in a section following the basis set includes.
 */
void Orca::WriteECP(std::ofstream& inp_file, const QMMolecule& qmatoms) {

  inp_file << endl;
  std::vector<std::string> UniqueElements = qmatoms.FindUniqueElements();

  ECPBasisSet ecp;
  ecp.Load(_settings.get("ecp"));

  XTP_LOG(Log::error, *_pLog)
      << "Loaded Pseudopotentials " << _settings.get("ecp") << flush;

  for (const std::string& element_name : UniqueElements) {
    try {
      ecp.getElement(element_name);
    } catch (std::runtime_error& error) {
      XTP_LOG(Log::error, *_pLog)
          << "No pseudopotential for " << element_name << " available" << flush;
      continue;
    }
    const ECPElement& element = ecp.getElement(element_name);

    inp_file << "\n"
             << "NewECP"
             << " " << element_name << endl;
    inp_file << "N_core"
             << " " << element.getNcore() << endl;
    inp_file << "lmax"
             << " " << EnumToString(element.getLmax()) << endl;
    // For Orca the order doesn't matter but let's write it in ascending order
    // write remaining shells in ascending order s,p,d...
    for (Index i = 0; i <= Index(element.getLmax()); i++) {
      for (const ECPShell& shell : element) {
        if (Index(shell.getL()) == i) {
          // shell type, number primitives, scale factor
          inp_file << xtp::EnumToString(shell.getL()) << " " << shell.getSize()
                   << endl;
          Index sh_idx = 0;
          for (const ECPGaussianPrimitive& gaussian : shell) {
            sh_idx++;
            inp_file << sh_idx << " " << gaussian._decay << " "
                     << gaussian._contraction << " " << gaussian._power << endl;
          }
        }
      }
    }
    inp_file << "end\n "
             << "\n"
             << endl;
  }
  return;
}

void Orca::WriteChargeOption() {
  this->_settings.add("orca.pointcharges", "\"background.crg\"");
}

/* For QM/MM the molecules in the MM environment are represented by
 * their atomic partial charge distributions. ORCA expects them in
 * q,x,y,z format in a separate file "background.crg"
 */
void Orca::WriteBackgroundCharges() {

  std::ofstream crg_file;
  std::string _crg_file_name_full = _run_dir + "/background.crg";
  crg_file.open(_crg_file_name_full);
  Index total_background = 0;

  for (const std::unique_ptr<StaticSite>& site : _externalsites) {
    if (site->getCharge() != 0.0) {
      total_background++;
    }
    std::vector<MinimalMMCharge> split_multipoles = SplitMultipoles(*site);
    total_background += split_multipoles.size();
  }  // counting only

  crg_file << total_background << endl;
  boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
  // now write
  for (const std::unique_ptr<StaticSite>& site : _externalsites) {
    Eigen::Vector3d pos = site->getPos() * tools::conv::bohr2ang;
    string sitestring =
        boost::str(fmt % site->getCharge() % pos.x() % pos.y() % pos.z());
    if (site->getCharge() != 0.0) {
      crg_file << sitestring << endl;
    }
    std::vector<MinimalMMCharge> split_multipoles = SplitMultipoles(*site);
    for (const auto& mpoles : split_multipoles) {
      Eigen::Vector3d pos2 = mpoles._pos * tools::conv::bohr2ang;
      string multipole =
          boost::str(fmt % mpoles._q % pos2.x() % pos2.y() % pos2.z());
      crg_file << multipole << endl;
    }
  }

  return;
}

/**
 * Prepares the *.inp file from a vector of segments
 * Appends a guess constructed from monomer orbitals if supplied, Not
 * implemented yet
 */
bool Orca::WriteInputFile(const Orbitals& orbitals) {

  std::vector<std::string> results;
  std::string temp_suffix = "/id";
  std::string scratch_dir_backup = _scratch_dir;
  std::ofstream inp_file;
  std::string inp_file_name_full = _run_dir + "/" + _input_file_name;
  inp_file.open(inp_file_name_full);
  // header
  inp_file << "* xyz  " << _charge << " " << _spin << endl;
  Index threads = OPENMP::getMaxThreads();
  const QMMolecule& qmatoms = orbitals.QMAtoms();
  // put coordinates
  WriteCoordinates(inp_file, qmatoms);
  // add parallelization info
  inp_file << "%pal\n"
           << "nprocs " << threads << "\nend"
           << "\n"
           << endl;
  // basis set info
  std::string el_file_name = _run_dir + "/" + "system.bas";
  WriteBasisset(qmatoms, _basisset_name, el_file_name);
  inp_file << "%basis\n";
  inp_file << "GTOName"
           << " "
           << "="
           << "\"system.bas\";" << endl;
  if (_settings.has_key("auxbasisset")) {
    std::string aux_file_name = _run_dir + "/" + "system.aux";
    std::string auxbasisset_name = _settings.get("auxbasisset");
    WriteBasisset(qmatoms, auxbasisset_name, aux_file_name);
    inp_file << "GTOAuxName"
             << " "
             << "="
             << "\"system.aux\";" << endl;
  }  // write_auxbasis set

  // ECPs
  if (_settings.has_key("ecp")) {
    WriteECP(inp_file, qmatoms);
  }
  inp_file << "end\n "
           << "\n"
           << endl;  // This end is for the basis set block
  if (_settings.get<bool>("write_charges")) {
    WriteBackgroundCharges();
  }

  // Write Orca section specified by the user
  for (const auto& prop : this->_settings.property("orca")) {
    const std::string& prop_name = prop.name();
    if (prop_name == "pointcharges") {
      _options += this->CreateInputSection("orca.pointcharges", true);
    } else if (prop_name != "method") {
      _options += this->CreateInputSection("orca." + prop_name);
    }
  }
  // Write main DFT method
  _options += this->WriteMethod();
  inp_file << _options;
  inp_file.close();
  // and now generate a shell script to run both jobs, if neccessary

  XTP_LOG(Log::info, *_pLog)
      << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;
  _scratch_dir = scratch_dir_backup + temp_suffix;
  WriteShellScript();
  _scratch_dir = scratch_dir_backup;
  return true;
}

bool Orca::WriteShellScript() {
  ofstream shell_file;
  std::string shell_file_name_full = _run_dir + "/" + _shell_file_name;
  shell_file.open(shell_file_name_full);
  shell_file << "#!/bin/bash" << endl;
  shell_file << "mkdir -p " << _scratch_dir << endl;

  if (_settings.get<bool>("read_guess")) {
    if (!(boost::filesystem::exists(_run_dir + "/molA.gbw") &&
          boost::filesystem::exists(_run_dir + "/molB.gbw"))) {
      throw runtime_error(
          "Using guess relies on a molA.gbw and a molB.gbw file being in the "
          "directory.");
    }
    shell_file << _settings.get("executable")
               << "_mergefrag molA.gbw molB.gbw dimer.gbw > merge.log" << endl;
  }
  shell_file << _settings.get("executable") << " " << _input_file_name << " > "
             << _log_file_name << endl;  //" 2> run.error" << endl;
  shell_file.close();
  return true;
}

/**
 * Runs the Orca job.
 */
bool Orca::Run() {

  XTP_LOG(Log::error, *_pLog) << "Running Orca job" << flush;

  if (std::system(nullptr)) {

    std::string command = "cd " + _run_dir + "; sh " + _shell_file_name;
    Index check = std::system(command.c_str());
    if (check == -1) {
      XTP_LOG(Log::error, *_pLog)
          << _input_file_name << " failed to start" << flush;
      return false;
    }
    if (CheckLogFile()) {
      XTP_LOG(Log::error, *_pLog) << "Finished Orca job" << flush;
      return true;
    } else {
      XTP_LOG(Log::error, *_pLog) << "Orca job failed" << flush;
    }
  } else {
    XTP_LOG(Log::error, *_pLog)
        << _input_file_name << " failed to start" << flush;
    return false;
  }

  return true;
}

/**
 * Cleans up after the Orca job
 */
void Orca::CleanUp() {

  if (_settings.get<bool>("read_guess")) {
    remove((_run_dir + "/" + "molA.gbw").c_str());
    remove((_run_dir + "/" + "molB.gbw").c_str());
    remove((_run_dir + "/" + "dimer.gbw").c_str());
  }
  // cleaning up the generated files
  if (_cleanup.size() != 0) {
    tools::Tokenizer tok_cleanup(_cleanup, ",");
    std::vector<std::string> cleanup_info;
    tok_cleanup.ToVector(cleanup_info);
    for (const std::string& substring : cleanup_info) {
      if (substring == "inp") {
        std::string file_name = _run_dir + "/" + _input_file_name;
        remove(file_name.c_str());
      }

      if (substring == "bas") {
        std::string file_name = _run_dir + "/system.bas";
        remove(file_name.c_str());
      }

      if (substring == "log") {
        std::string file_name = _run_dir + "/" + _log_file_name;
        remove(file_name.c_str());
      }

      if (substring == "gbw") {
        std::string file_name = _run_dir + "/" + _mo_file_name;
        remove(file_name.c_str());
      }

      if (substring == "ges") {
        std::string file_name = _run_dir + "/system.ges";
        remove(file_name.c_str());
      }
      if (substring == "prop") {
        std::string file_name = _run_dir + "/system.prop";
        remove(file_name.c_str());
      }
    }
  }
  return;
}

StaticSegment Orca::GetCharges() const {

  StaticSegment result("charges", 0);

  XTP_LOG(Log::error, *_pLog) << "Parsing " << _log_file_name << flush;
  std::string log_file_name_full = _run_dir + "/" + _log_file_name;
  std::string line;

  std::ifstream input_file(log_file_name_full);
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);
    GetCoordinates(result, line, input_file);

    std::string::size_type charge_pos = line.find("CHELPG Charges");

    if (charge_pos != std::string::npos) {
      XTP_LOG(Log::error, *_pLog) << "Getting charges" << flush;
      getline(input_file, line);
      std::vector<std::string> row = GetLineAndSplit(input_file, "\t ");
      Index nfields = Index(row.size());
      bool hasAtoms = result.size() > 0;
      while (nfields == 4) {
        Index atom_id = boost::lexical_cast<Index>(row.at(0));
        std::string atom_type = row.at(1);
        double atom_charge = boost::lexical_cast<double>(row.at(3));
        row = GetLineAndSplit(input_file, "\t ");
        nfields = Index(row.size());
        if (hasAtoms) {
          StaticSite& temp = result.at(atom_id);
          if (temp.getElement() != atom_type) {
            throw std::runtime_error(
                "Getting charges failed. Mismatch in elemts:" +
                temp.getElement() + " vs " + atom_type);
          }
          temp.setCharge(atom_charge);
        } else {
          StaticSite temp =
              StaticSite(atom_id, atom_type, Eigen::Vector3d::Zero());
          temp.setCharge(atom_charge);
          result.push_back(temp);
        }
      }
    }
  }
  return result;
}

Eigen::Matrix3d Orca::GetPolarizability() const {
  std::string line;
  ifstream input_file((_run_dir + "/" + _log_file_name));
  bool has_pol = false;

  Eigen::Matrix3d pol = Eigen::Matrix3d::Zero();
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);

    std::string::size_type pol_pos = line.find("THE POLARIZABILITY TENSOR");
    if (pol_pos != std::string::npos) {
      XTP_LOG(Log::error, *_pLog) << "Getting polarizability" << flush;
      getline(input_file, line);
      getline(input_file, line);
      getline(input_file, line);

      if (line.find("The raw cartesian tensor (atomic units)") ==
          std::string::npos) {
        throw std::runtime_error(
            "Could not find cartesian polarization tensor");
      }

      for (Index i = 0; i < 3; i++) {
        getline(input_file, line);
        tools::Tokenizer tok2(line, " ");
        std::vector<std::string> values = tok2.ToVector();
        if (values.size() != 3) {
          throw std::runtime_error("polarization line " + line +
                                   " cannot be parsed");
        }
        Eigen::Vector3d row;
        row << std::stod(values[0]), std::stod(values[1]), std::stod(values[2]);
        pol.row(i) = row;
      }

      has_pol = true;
    }
  }
  if (!has_pol) {
    throw std::runtime_error("Could not find polarization in logfile");
  }
  return pol;
}

bool Orca::ParseLogFile(Orbitals& orbitals) {
  bool found_success = false;
  orbitals.setQMpackage(getPackageName());
  orbitals.setDFTbasisName(_basisset_name);
  if (_settings.has_key("ecp")) {
    orbitals.setECPName(_settings.get("ecp"));
  }

  XTP_LOG(Log::error, *_pLog) << "Parsing " << _log_file_name << flush;
  std::string log_file_name_full = _run_dir + "/" + _log_file_name;
  // check if LOG file is complete
  if (!CheckLogFile()) {
    return false;
  }
  std::map<Index, double> energies;
  std::map<Index, double> occupancy;

  std::string line;
  Index levels = 0;
  Index number_of_electrons = 0;
  std::vector<std::string> results;

  std::ifstream input_file(log_file_name_full);

  if (input_file.fail()) {
    XTP_LOG(Log::error, *_pLog)
        << "File " << log_file_name_full << " not found " << flush;
    return false;
  } else {
    XTP_LOG(Log::error, *_pLog)
        << "Reading Coordinates and occupationnumbers and energies from "
        << log_file_name_full << flush;
  }
  // Coordinates of the final configuration depending on whether it is an
  // optimization or not

  QMMolecule& mol = orbitals.QMAtoms();
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);

    GetCoordinates(mol, line, input_file);

    std::string::size_type energy_pos = line.find("FINAL SINGLE");
    if (energy_pos != std::string::npos) {

      boost::algorithm::split(results, line, boost::is_any_of(" "),
                              boost::algorithm::token_compress_on);
      std::string energy = results[4];
      boost::trim(energy);
      orbitals.setQMEnergy(boost::lexical_cast<double>(energy));
      XTP_LOG(Log::error, *_pLog) << (boost::format("QM energy[Hrt]: %4.6f ") %
                                      orbitals.getDFTTotalEnergy())
                                         .str()
                                  << flush;
    }

    std::string::size_type HFX_pos = line.find("Fraction HF Exchange ScalHFX");
    if (HFX_pos != std::string::npos) {
      boost::algorithm::split(results, line, boost::is_any_of(" "),
                              boost::algorithm::token_compress_on);
      double ScaHFX = boost::lexical_cast<double>(results.back());
      orbitals.setScaHFX(ScaHFX);
      XTP_LOG(Log::error, *_pLog)
          << "DFT with " << ScaHFX << " of HF exchange!" << flush;
    }

    std::string::size_type dim_pos = line.find("Basis Dimension");
    if (dim_pos != std::string::npos) {
      boost::algorithm::split(results, line, boost::is_any_of(" "),
                              boost::algorithm::token_compress_on);
      std::string dim =
          results[4];  // The 4th element of results vector is the Basis Dim
      boost::trim(dim);
      levels = boost::lexical_cast<Index>(dim);
      XTP_LOG(Log::info, *_pLog) << "Basis Dimension: " << levels << flush;
      XTP_LOG(Log::info, *_pLog) << "Energy levels: " << levels << flush;
    }

    std::string::size_type OE_pos = line.find("ORBITAL ENERGIES");
    if (OE_pos != std::string::npos) {

      number_of_electrons = 0;
      getline(input_file, line);
      getline(input_file, line);
      getline(input_file, line);
      if (line.find("E(Eh)") == std::string::npos) {
        XTP_LOG(Log::error, *_pLog)
            << "Warning: Orbital Energies not found in log file" << flush;
      }
      for (Index i = 0; i < levels; i++) {
        results = GetLineAndSplit(input_file, " ");
        std::string no = results[0];
        boost::trim(no);
        Index levelnumber = boost::lexical_cast<Index>(no);
        if (levelnumber != i) {
          XTP_LOG(Log::error, *_pLog) << "Have a look at the orbital energies "
                                         "something weird is going on"
                                      << flush;
        }
        std::string oc = results[1];
        boost::trim(oc);
        double occ = boost::lexical_cast<double>(oc);
        // We only count alpha electrons, each orbital must be empty or doubly
        // occupied
        if (occ == 2 || occ == 1) {
          number_of_electrons++;
          occupancy[i] = occ;
          if (occ == 1) {
            XTP_LOG(Log::error, *_pLog)
                << "Watch out! No distinction between alpha and beta "
                   "electrons. Check if occ = 1 is suitable for your "
                   "calculation "
                << flush;
          }
        } else if (occ == 0) {
          occupancy[i] = occ;
        } else {
          throw runtime_error(
              "Only empty or doubly occupied orbitals are allowed not "
              "running the right kind of DFT calculation");
        }
        std::string e = results[2];
        boost::trim(e);
        energies[i] = boost::lexical_cast<double>(e);
      }
    }

    std::string::size_type success =
        line.find("*                     SUCCESS                       *");
    if (success != std::string::npos) {
      found_success = true;
    }
  }

  XTP_LOG(Log::info, *_pLog)
      << "Alpha electrons: " << number_of_electrons << flush;
  Index occupied_levels = number_of_electrons;
  Index unoccupied_levels = levels - occupied_levels;
  XTP_LOG(Log::info, *_pLog) << "Occupied levels: " << occupied_levels << flush;
  XTP_LOG(Log::info, *_pLog)
      << "Unoccupied levels: " << unoccupied_levels << flush;

  /************************************************************/

  // copying information to the orbitals object

  orbitals.setBasisSetSize(levels);
  orbitals.setNumberOfAlphaElectrons(number_of_electrons);
  orbitals.setNumberOfOccupiedLevels(occupied_levels);

  // copying energies to a vector
  orbitals.MOs().eigenvalues().resize(levels);
  //_level = 1;
  for (Index i = 0; i < levels; i++) {
    orbitals.MOs().eigenvalues()[i] = energies[i];
  }

  XTP_LOG(Log::error, *_pLog) << "Done reading Log file" << flush;

  return found_success;
}
template <class T>
void Orca::GetCoordinates(T& mol, string& line, ifstream& input_file) const {
  std::string::size_type coordinates_pos =
      line.find("CARTESIAN COORDINATES (ANGSTROEM)");

  using Atom = typename std::iterator_traits<typename T::iterator>::value_type;

  if (coordinates_pos != std::string::npos) {
    XTP_LOG(Log::error, *_pLog) << "Getting the coordinates" << flush;
    bool has_QMAtoms = mol.size() > 0;
    // three garbage lines
    getline(input_file, line);
    // now starts the data in format
    // _id type Qnuc x y z
    vector<string> row = GetLineAndSplit(input_file, "\t ");
    Index nfields = Index(row.size());
    Index atom_id = 0;
    while (nfields == 4) {
      string atom_type = row.at(0);
      double x = boost::lexical_cast<double>(row.at(1));
      double y = boost::lexical_cast<double>(row.at(2));
      double z = boost::lexical_cast<double>(row.at(3));
      row = GetLineAndSplit(input_file, "\t ");
      nfields = Index(row.size());
      Eigen::Vector3d pos(x, y, z);
      pos *= tools::conv::ang2bohr;
      if (has_QMAtoms == false) {
        mol.push_back(Atom(atom_id, atom_type, pos));
      } else {
        Atom& pAtom = mol.at(atom_id);
        pAtom.setPos(pos);
      }
      atom_id++;
    }
  }
}

bool Orca::CheckLogFile() {
  // check if the log file exists
  ifstream input_file(_run_dir + "/" + _log_file_name);

  if (input_file.fail()) {
    XTP_LOG(Log::error, *_pLog) << "Orca LOG is not found" << flush;
    return false;
  };

  std::string line;
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);
    std::string::size_type error = line.find("FATAL ERROR ENCOUNTERED");
    if (error != std::string::npos) {
      XTP_LOG(Log::error, *_pLog) << "ORCA encountered a fatal error, maybe a "
                                     "look in the log file may help."
                                  << flush;
      return false;
    }
    error = line.find(
        "mpirun detected that one or more processes exited with non-zero "
        "status");
    if (error != std::string::npos) {
      XTP_LOG(Log::error, *_pLog)
          << "ORCA had an mpi problem, maybe your openmpi version is not good."
          << flush;
      return false;
    }
  }
  return true;
}

// Parses the Orca gbw file and stores data in the Orbitals object

bool Orca::ParseMOsFile(Orbitals& orbitals) {
  if (!CheckLogFile()) {
    return false;
  }
  std::vector<double> coefficients;
  Index basis_size = orbitals.getBasisSetSize();
  if (basis_size == 0) {
    throw runtime_error(
        "Basis size not set, calculator does not parse log file first");
  }

  XTP_LOG(Log::error, *_pLog)
      << "Reading the gbw file, this may or may not work so be careful: "
      << flush;
  ifstream infile;
  infile.open(_run_dir + "/" + _mo_file_name, ios::binary | ios::in);
  if (!infile) {
    throw runtime_error("Could not open " + _mo_file_name + " file");
  }
  infile.seekg(24, ios::beg);
  std::array<char, 8> buffer;
  infile.read(buffer.data(), 8);
  if (!infile) {
    infile.close();
    return false;
  }
  Index offset = *((Index*)buffer.data());

  infile.seekg(offset, ios::beg);
  infile.read(buffer.data(), 4);
  if (!infile) {
    infile.close();
    return false;
  }
  int op_read = *((int*)buffer.data());
  infile.seekg(offset + 4, ios::beg);
  infile.read(buffer.data(), 4);
  if (!infile) {
    infile.close();
    return false;
  }
  int dim_read = *((int*)buffer.data());
  infile.seekg(offset + 8, ios::beg);
  XTP_LOG(Log::info, *_pLog) << "Number of operators: " << op_read
                             << " Basis dimension: " << dim_read << flush;
  Index n = op_read * dim_read * dim_read;
  for (Index i = 0; i < n; i++) {
    infile.read(buffer.data(), 8);
    if (!infile) {
      infile.close();
      return false;
    }
    double mocoeff = *((double*)buffer.data());
    coefficients.push_back(mocoeff);
  }

  infile.close();
  // i -> MO, j -> AO
  orbitals.MOs().eigenvectors().resize(basis_size, basis_size);
  for (Index i = 0; i < basis_size; i++) {
    for (Index j = 0; j < basis_size; j++) {
      orbitals.MOs().eigenvectors()(j, i) = coefficients[j * basis_size + i];
    }
  }
  ReorderOutput(orbitals);
  XTP_LOG(Log::error, *_pLog) << "Done parsing" << flush;
  return true;
}

std::string Orca::indent(const double& number) {
  std::stringstream ssnumber;
  if (number >= 0) {
    ssnumber << "    ";
  } else {
    ssnumber << "   ";
  }
  ssnumber << setiosflags(ios::fixed) << setprecision(15) << std::scientific
           << number;
  std::string snumber = ssnumber.str();
  return snumber;
}

std::string Orca::CreateInputSection(const std::string& key,
                                     bool single_line) const {
  std::stringstream stream;
  std::string section = key.substr(key.find(".") + 1);
  stream << "%" << section;
  if (single_line) {
    stream << " " << _settings.get(key) << "\n";
  } else {
    stream << "\n"
           << this->_settings.get(key) << "\n"
           << "end\n";
  }

  return stream.str();
}

std::string Orca::WriteMethod() const {
  std::stringstream stream;
  std::string opt = (_settings.get<bool>("optimize")) ? "Opt" : "";
  const tools::Property& orca = _settings.property("orca");
  std::string user_method =
      (orca.exists("method")) ? orca.get("method").as<std::string>() : "";
  std::string convergence = "";
  if (!orca.exists("scf")) {
    convergence =
        this->_convergence_map.at(_settings.get("convergence_tightness"));
  }
  stream << "! DFT " << this->GetOrcaFunctionalName() << " " << convergence
         << " " << opt
         << " "
         // additional properties provided by the user
         << user_method << "\n";
  return stream.str();
}

std::string Orca::GetOrcaFunctionalName() const {

  char* votca_share = getenv("VOTCASHARE");
  if (votca_share == nullptr) {
    return _settings.get("functional");
  } else {
    tools::Property all_functionals;

    auto xml_file = std::string(getenv("VOTCASHARE")) +
                    std::string("/xtp/data/orca_functional_names.xml");

    all_functionals.LoadFromXML(xml_file);

    const tools::Property& functional_names =
        all_functionals.get("functionals");

    std::string input_name = _settings.get("functional");
    // Some functionals have a composed named separated by a space
    // In the case just look for the first part
    std::size_t plus = input_name.find(' ');
    if (plus != std::string::npos) {
      input_name = input_name.substr(0, plus);
    }

    if (functional_names.exists(input_name)) {
      return functional_names.get(input_name).as<std::string>();
    } else {
      std::ostringstream oss;
      oss << "The libxc functional \"" << input_name << "\"\n"
          << "doesn't seem to have a corresponding name in Orca.\n"
          << "Check the "
          << "\"${VOTCASHARE}/xtp/data/orca_functional_names.xml\""
          << "file for the whole list of known libxc/orca functionals";
      throw runtime_error(oss.str());
    }
  }
}

}  // namespace xtp
}  // namespace votca
