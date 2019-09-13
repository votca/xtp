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

#include "gaussian.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <iomanip>
#include <stdio.h>
#include <votca/tools/constants.h>
#include <votca/xtp/ecpaobasis.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {
using namespace std;

void Gaussian::Initialize(tools::Property& options) {

  // GAUSSIAN file names
  std::string fileName = "system";
  _input_file_name = fileName + ".com";
  _log_file_name = fileName + ".log";
  _shell_file_name = fileName + ".sh";
  _mo_file_name = "fort.7";
  std::string key = "package";
  ParseCommonOptions(options);

  if (options.exists(key + ".vdWRadii")) {
    _vdWfooter = options.get(key + ".vdWRadii").as<std::string>();
  }

  /* G09 by default deletes functions from the basisset according to some
   * criterion based on, a.o., the contraction coefficients. This can lead
   * to inconsistencies when MOs are later used in VOTCA's GWBSE modules
   * (and other post-processing routines). G09's default can be modified
   * by the keywork int=nobasistransform. This will add this keyword
   * automatically to the _options string for runs with G09.
   */
  if (_executable == "g09") {
    std::string::size_type basistransform_pos =
        (boost::algorithm::to_lower_copy(_options)).find("nobasistransform");
    if (basistransform_pos == std::string::npos) {
      _options = _options + " int=nobasistransform ";
    }
  }
  std::string::size_type iop_pos;
  // check if the guess keyword is present, if yes, append the guess later
  if (_write_guess) {
    iop_pos = _options.find("cards");
    if (iop_pos != std::string::npos) {
      _options = _options + " cards ";
    }
  }

  // check if the basis set is available ("/gen")
  iop_pos = _options.find("gen");
  if (iop_pos != std::string::npos) {
    _write_basis_set = true;
  } else {
    _write_basis_set = false;
  }

  // check if pseudopotentials are required ("pseudo")
  iop_pos = _options.find("pseudo");
  if (iop_pos != std::string::npos) {
    _write_pseudopotentials = true;
  } else {
    _write_pseudopotentials = false;
  }
}

void Gaussian::WriteChargeOption() {
  std::string::size_type iop_pos = _options.find("charge");
  if (iop_pos == std::string::npos) {
    std::string::size_type pos = _options.find('\n');
    if (pos != std::string::npos) {
      _options.insert(pos, " charge");
    } else {
      _options = _options + " charge";
    }
  }
}

/* Custom basis sets are written on a per-element basis to
 * 'elementname'.gbs files, which are then included in the
 * Gaussian input file using @'elementname'.gbs
 */
void Gaussian::WriteBasisset(std::ofstream& com_file,
                             const QMMolecule& qmatoms) {

  std::vector<std::string> UniqueElements = qmatoms.FindUniqueElements();
  BasisSet bs;
  bs.Load(_basisset_name);
  XTP_LOG(logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;

  for (const std::string& element_name : UniqueElements) {

    const Element& element = bs.getElement(element_name);
    /* Write each basis set to a element_name.gbs file
     * and include the gbs file in the com-file via Gaussian's @ function
     */
    std::ofstream el_file;
    std::string el_file_name = _run_dir + "/" + element_name + ".gbs";

    el_file.open(el_file_name);
    // element name, [possibly indeces of centers], zero to indicate the end
    com_file << "@" << element_name << ".gbs" << endl;
    el_file << element_name << " 0" << endl;
    for (const Shell& shell : element) {
      // gaussian can only use S,P,SP,D,F,G shells so we split them up if not SP
      if (shell.getType() == "SP" || !shell.isCombined()) {
        // shell type, number primitives, scale factor
        el_file << shell.getType() << " " << shell.getSize() << " "
                << FortranFormat(shell.getScale()) << endl;
        for (const GaussianPrimitive& gaussian : shell) {
          el_file << FortranFormat(gaussian._decay);
          for (const double& contraction : gaussian._contraction) {
            if (contraction != 0.0) {
              el_file << " " << FortranFormat(contraction);
            }
          }
          el_file << endl;
        }
      } else {
        string type = shell.getType();
        for (unsigned i = 0; i < type.size(); ++i) {
          string subtype = string(type, i, 1);
          el_file << subtype << " " << shell.getSize() << " "
                  << FortranFormat(shell.getScale()) << endl;

          for (const GaussianPrimitive& gaussian : shell) {
            el_file << FortranFormat(gaussian._decay);
            el_file << " "
                    << FortranFormat(gaussian._contraction[FindLmax(subtype)]);
          }
          el_file << endl;
        }
      }
    }

    el_file << "****\n";
    el_file.close();
  }

  com_file << endl;
  return;
}

/* If custom ECPs are used, they need to be specified in the input file
 * in a section following the basis set includes.
 */
void Gaussian::WriteECP(std::ofstream& com_file, const QMMolecule& qmatoms) {
  std::vector<std::string> UniqueElements = qmatoms.FindUniqueElements();

  ECPBasisSet ecp;
  ecp.Load(_ecp_name);

  XTP_LOG(logDEBUG, *_pLog) << "Loaded Pseudopotentials " << _ecp_name << flush;

  for (const std::string& element_name : UniqueElements) {
    try {
      ecp.getElement(element_name);
    } catch (std::runtime_error& error) {
      XTP_LOG(logDEBUG, *_pLog)
          << "No pseudopotential for " << element_name << " available" << flush;
      continue;
    }
    const ECPElement& element = ecp.getElement(element_name);
    // element name, [possibly indeces of centers], zero to indicate the end
    com_file << element_name << " 0\n"
             << _ecp_name << " " << element.getLmax() << " "
             << element.getNcore() << endl;

    for (const ECPShell& shell : element) {
      // shell type, number primitives, scale factor
      com_file << shell.getType() << endl;
      com_file << shell.getSize() << endl;

      for (const ECPGaussianPrimitive& gaussian : shell) {
        com_file << gaussian._power << " " << FortranFormat(gaussian._decay)
                 << " " << FortranFormat(gaussian._contraction) << endl;
      }
    }
  }
  com_file << endl;
  return;
}

/* For QM/MM the molecules in the MM environment are represented by
 * their atomic partial charge distributions. Triggered by the option
 * keyword "charge" Gaussian expects them in x,y,z,q format in the
 * input file. In g03 AFTER basis sets and ECPs, in g09 BEFORE.
 */

void Gaussian::WriteBackgroundCharges(std::ofstream& com_file) {

  boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
  for (const std::unique_ptr<StaticSite>& site : _externalsites) {
    Eigen::Vector3d pos = site->getPos() * tools::conv::bohr2ang;
    string sitestring =
        boost::str(fmt % pos.x() % pos.y() % pos.z() % site->getCharge());
    if (site->getCharge() != 0.0) com_file << sitestring << endl;

    std::vector<MinimalMMCharge> split_multipoles = SplitMultipoles(*site);
    for (const auto& mpoles : split_multipoles) {
      Eigen::Vector3d pos = mpoles._pos * tools::conv::bohr2ang;
      string multipole =
          boost::str(fmt % pos.x() % pos.y() % pos.z() % mpoles._q);
      com_file << multipole << endl;
    }
  }

  com_file << endl;
  return;
}

/* An initial guess for the electron density can be provided by
 * a set of molecular orbital coefficients in the input file,
 * triggered by the 'guess=cards' keyword. This MUST be done in
 * Fortran fixed format 5D15.8. The information about the guess
 * itself is taken from a prepared orbitals object.
 */
void Gaussian::WriteGuess(const Orbitals& orbitals_guess,
                          std::ofstream& com_file) {
  Eigen::MatrixXd MOs = ReorderMOsBack(orbitals_guess);
  com_file << "(5D15.8)" << endl;
  int level = 1;
  int ncolumns = 5;
  for (int i = 0; i < MOs.cols(); ++i) {
    com_file << setw(5) << level << endl;
    Eigen::VectorXd mr = MOs.col(i);
    int column = 1;
    for (unsigned j = 0; j < mr.size(); ++j) {
      com_file << FortranFormat(mr[j]);
      if (column == ncolumns) {
        com_file << std::endl;
        column = 0;
      }
      column++;
    }
    level++;
    if (column != 1) com_file << endl;
  }
  com_file << 0 << endl;
  return;
}

/* Coordinates are written in standard Element,x,y,z format to the
 * input file.
 */
void Gaussian::WriteCoordinates(std::ofstream& com_file,
                                const QMMolecule& qmatoms) {
  for (const QMAtom& atom : qmatoms) {
    Eigen::Vector3d pos = atom.getPos() * tools::conv::bohr2ang;
    com_file << setw(3) << atom.getElement() << setw(12)
             << setiosflags(ios::fixed) << setprecision(5) << pos.x()
             << setw(12) << setiosflags(ios::fixed) << setprecision(5)
             << pos.y() << setw(12) << setiosflags(ios::fixed)
             << setprecision(5) << pos.z() << endl;
  }
  com_file << endl;
  return;
}

/* Standard Gaussian Header is written to the input file, with checkpoint,
 * memory, shared processor request, option string containing all
 * relevant keywords, charge, and spin information.
 */
void Gaussian::WriteHeader(std::ofstream& com_file) {
  if (_memory.size()) com_file << "%mem=" << _memory << endl;

  int threads = OPENMP::getMaxThreads();
  if (threads > 0) com_file << "%nprocshared=" << threads << endl;
  if (_options.size()) com_file << _options << endl;

  com_file << endl;
  com_file << "TITLE ";

  com_file << endl << endl;
  com_file << setw(2) << _charge << setw(2) << _spin << endl;
  return;
}

/**
 * Prepares the com file from a vector of segments
 * Appends a guess constructed from monomer orbitals if supplied
 */
bool Gaussian::WriteInputFile(const Orbitals& orbitals) {

  std::string temp_suffix = "/id";
  std::string scratch_dir_backup = _scratch_dir;

  std::ofstream com_file;
  std::string com_file_name_full = _run_dir + "/" + _input_file_name;
  com_file.open(com_file_name_full);

  // header
  WriteHeader(com_file);

  const QMMolecule& qmatoms = orbitals.QMAtoms();

  WriteCoordinates(com_file, qmatoms);

  /* The order in which the following information has to appear
   * in the Gaussian input file is different from version g03 to
   * version g09. The newest version (g2016) is not supported.
   */
  if (_executable == "g03") {

    // if we need to write basis sets, do it now
    if (_write_basis_set) WriteBasisset(com_file, qmatoms);

    // write ECPs
    if (_write_pseudopotentials) WriteECP(com_file, qmatoms);

    // write the background charges
    if (_write_charges) WriteBackgroundCharges(com_file);

    // write inital guess
    if (_write_guess) {
      WriteGuess(orbitals, com_file);
    }

  } else if (_executable == "g09") {

    // write the background charges
    // if (_write_charges) WriteBackgroundCharges(_com_file, qmatoms);
    if (_write_charges) WriteBackgroundCharges(com_file);
    // if we need to write basis sets, do it now
    if (_write_basis_set) WriteBasisset(com_file, qmatoms);

    // write ECPs
    if (_write_pseudopotentials) WriteECP(com_file, qmatoms);

    // write inital guess
    if (_write_guess) {
      WriteGuess(orbitals, com_file);
    }

  } else {
    throw std::runtime_error(
        "Gaussian executable unknown. Must be either g03 or g09.");
  }

  com_file << _vdWfooter << endl;

  com_file << endl;
  com_file.close();
  // and now generate a shell script to run both jobs
  XTP_LOG(logDEBUG, *_pLog)
      << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;

  _scratch_dir = scratch_dir_backup + temp_suffix;
  WriteShellScript();
  _scratch_dir = scratch_dir_backup;

  return true;
}

/* Gaussian will be executed within a shell in order to set some
 * environment variables for the local SCRATCH directory and
 * (legacy mode) running a second instance for AO matrix of Vxc
 * using patched g03. This function writes the shell script.
 */
bool Gaussian::WriteShellScript() {
  std::ofstream shell_file;

  std::string shell_file_name_full = _run_dir + "/" + _shell_file_name;

  shell_file.open(shell_file_name_full);

  shell_file << "#!/bin/tcsh" << endl;
  shell_file << "mkdir -p " << _scratch_dir << endl;
  shell_file << "setenv GAUSS_SCRDIR " << _scratch_dir << endl;
  shell_file << _executable << " " << _input_file_name << endl;
  shell_file.close();

  return true;
}

/**
 * Runs the Gaussian job.
 */
bool Gaussian::Run() {
  XTP_LOG(logDEBUG, *_pLog) << "GAUSSIAN: running [" << _executable << " "
                            << _input_file_name << "]" << flush;

  if (std::system(NULL)) {
    // if scratch is provided, run the shell script;
    // otherwise run gaussian directly and rely on global variables
    std::string command;
    if (_scratch_dir.size() != 0) {
      command = "cd " + _run_dir + "; tcsh " + _shell_file_name;
      //            _command  = "cd " + _run_dir + "; mkdir -p " + _scratch_dir
      //            +"; " + _executable + " " + _input_file_name;
    } else {
      command = "cd " + _run_dir + "; mkdir -p $GAUSS_SCRDIR; " + _executable +
                " " + _input_file_name;
    }
    int check = std::system(command.c_str());
    if (check == -1) {
      XTP_LOG(logERROR, *_pLog)
          << _input_file_name << " failed to start" << flush;
      return false;
    }
    if (CheckLogFile()) {
      XTP_LOG(logDEBUG, *_pLog) << "GAUSSIAN: finished job" << flush;
      return true;
    } else {
      XTP_LOG(logDEBUG, *_pLog) << "GAUSSIAN: job failed" << flush;
    }
  } else {
    XTP_LOG(logERROR, *_pLog)
        << _input_file_name << " failed to start" << flush;
    return false;
  }
  return true;
}

/**
 * Cleans up after the Gaussian job
 */
void Gaussian::CleanUp() {

  // cleaning up the generated files
  if (_cleanup.size() != 0) {

    XTP_LOG(logDEBUG, *_pLog) << "Removing " << _cleanup << " files" << flush;
    tools::Tokenizer tok_cleanup(_cleanup, ", ");
    std::vector<std::string> cleanup_info;
    tok_cleanup.ToVector(cleanup_info);
    for (const std::string& substring : cleanup_info) {

      if (substring == "com") {
        std::string file_name = _run_dir + "/" + _input_file_name;
        remove(file_name.c_str());
      }

      if (substring == "sh") {
        std::string file_name = _run_dir + "/" + _shell_file_name;
        remove(file_name.c_str());
      }

      if (substring == "log") {
        std::string file_name = _run_dir + "/" + _log_file_name;
        remove(file_name.c_str());
      }

      if (substring == "fort.7") {
        std::string file_name = _run_dir + "/" + substring;
        remove(file_name.c_str());
      }

      if (substring == "gbs" && _write_basis_set) {
        std::vector<std::string> fileswithfileending;
        boost::filesystem::recursive_directory_iterator fit(_run_dir);
        boost::filesystem::recursive_directory_iterator endit;
        while (fit != endit) {
          if (boost::filesystem::is_regular_file(*fit) &&
              fit->path().extension() == substring)
            fileswithfileending.push_back(fit->path().filename().string());
          ++fit;
        }
        for (const auto filename : fileswithfileending) {
          std::string file_name = _run_dir + "/" + filename;
          remove(file_name.c_str());
        }
      }
    }
  }
  return;
}

/**
 * Reads in the MO coefficients from a GAUSSIAN fort.7 file
 */
bool Gaussian::ParseMOsFile(Orbitals& orbitals) {
  std::map<int, std::vector<double> > coefficients;
  std::map<int, double> energies;

  std::string line;
  unsigned levels = 0;
  unsigned level = 0;
  unsigned basis_size = 0;

  std::string orb_file_name_full = _mo_file_name;
  if (_run_dir != "") orb_file_name_full = _run_dir + "/" + _mo_file_name;
  std::ifstream input_file(orb_file_name_full);

  if (input_file.fail()) {
    XTP_LOG(logERROR, *_pLog)
        << "File " << _mo_file_name << " with molecular orbitals is not found "
        << flush;
    return false;
  } else {
    XTP_LOG(logDEBUG, *_pLog) << "Reading MOs from " << _mo_file_name << flush;
  }

  // number of coefficients per line is  in the first line of the file (5D15.8)
  getline(input_file, line);
  std::vector<std::string> strs;
  boost::algorithm::split(strs, line, boost::is_any_of("(D)"));
  std::string format = strs.at(2);

  while (input_file) {

    getline(input_file, line);
    // if a line has an equality sign, must be energy
    std::string::size_type energy_pos = line.find("=");

    if (energy_pos != std::string::npos) {

      boost::trim(line);
      tools::Tokenizer tok(line, "\t =");
      std::vector<std::string> results = tok.ToVector();
      level = boost::lexical_cast<int>(results.front());
      boost::replace_first(results.back(), "D", "e");
      energies[level] = boost::lexical_cast<double>(results.back());
      levels++;

    } else {

      while (line.size() > 1) {
        std::string coefficient;
        coefficient.assign(line, 0, 15);
        boost::trim(coefficient);
        boost::replace_first(coefficient, "D", "e");
        double coef = boost::lexical_cast<double>(coefficient);
        coefficients[level].push_back(coef);
        line.erase(0, 15);
      }
    }
  }

  // some sanity checks
  XTP_LOG(logDEBUG, *_pLog) << "Energy levels: " << levels << flush;
  std::map<int, std::vector<double> >::iterator iter = coefficients.begin();
  basis_size = iter->second.size();

  for (iter = coefficients.begin()++; iter != coefficients.end(); iter++) {
    if (iter->second.size() != basis_size) {
      XTP_LOG(logERROR, *_pLog)
          << "Error reading " << _mo_file_name
          << ". Basis set size change from level to level." << flush;
      return false;
    }
  }

  XTP_LOG(logDEBUG, *_pLog) << "Basis set size: " << basis_size << flush;

  // copying information to the orbitals object
  orbitals.setBasisSetSize(basis_size);  // = _basis_size;

  // copying energies to the orbitals object
  Eigen::VectorXd& mo_energies = orbitals.MOs().eigenvalues();
  mo_energies.resize(levels);
  for (int i = 0; i < mo_energies.size(); i++) mo_energies[i] = energies[i + 1];

  // copying mo coefficients to the orbitals object
  Eigen::MatrixXd& mo_coefficients = orbitals.MOs().eigenvectors();
  mo_coefficients.resize(levels, basis_size);
  for (int i = 0; i < mo_coefficients.rows(); i++) {
    for (int j = 0; j < mo_coefficients.cols(); j++) {
      mo_coefficients(j, i) = coefficients[i + 1][j];
    }
  }

  ReorderOutput(orbitals);
  XTP_LOG(logDEBUG, *_pLog) << "GAUSSIAN: done reading MOs" << flush;

  return true;
}

bool Gaussian::CheckLogFile() const {

  // check if the log file exists
  boost::filesystem::path arg_path;
  char ch;

  std::string full_name = (arg_path / _run_dir / _log_file_name).c_str();
  ifstream input_file(full_name);

  if (input_file.fail()) {
    XTP_LOG(logERROR, *_pLog)
        << "GAUSSIAN: " << full_name << " is not found" << flush;
    return false;
  };

  input_file.seekg(0, ios_base::end);  // go to the EOF

  // get empty lines and end of lines out of the way
  do {
    input_file.seekg(-2, ios_base::cur);
    input_file.get(ch);
  } while (ch == '\n' || ch == ' ' || ch == '\t' ||
           (int)input_file.tellg() == -1);

  // get the beginning of the line or the file
  do {
    input_file.seekg(-2, ios_base::cur);
    input_file.get(ch);
  } while (ch != '\n' && (int)input_file.tellg() != -1);

  std::string line;
  getline(input_file, line);
  input_file.close();

  std::string::size_type success = line.find("Normal termination of Gaussian");
  if (success == std::string::npos) {
    XTP_LOG(logERROR, *_pLog)
        << "GAUSSIAN: " << full_name << " is incomplete" << flush;
    return false;
  } else {
    return true;
  }
}

StaticSegment Gaussian::GetCharges() const {
  XTP_LOG(logDEBUG, *_pLog) << "GAUSSIAN: parsing " << _log_file_name << flush;

  StaticSegment result("charges", 0);
  std::string log_file_name_full = _log_file_name;
  if (_run_dir != "") log_file_name_full = _run_dir + "/" + _log_file_name;

  // check if LOG file is complete
  if (!CheckLogFile()) throw std::runtime_error("logfile is not complete");
  std::string line;
  ifstream input_file(log_file_name_full);
  bool has_charges = false;
  std::vector<std::string> archive;
  while (input_file) {

    getline(input_file, line);
    boost::trim(line);
    GetArchive(archive, line, input_file);

    std::string::size_type charge_pos = line.find("Charges from ESP fit, RMS");
    if (charge_pos != std::string::npos) {
      has_charges = true;
      XTP_LOG(logDEBUG, *_pLog) << "Getting charges" << flush;
      getline(input_file, line);
      getline(input_file, line);

      std::vector<std::string> row = GetLineAndSplit(input_file, "\t ");
      int nfields = row.size();

      while (nfields == 3) {
        int atom_id = boost::lexical_cast<int>(row.at(0)) - 1;
        std::string atom_type = row.at(1);
        double atom_charge = boost::lexical_cast<double>(row.at(2));
        row = GetLineAndSplit(input_file, "\t ");
        nfields = row.size();

        StaticSite temp =
            StaticSite(atom_id, atom_type, Eigen::Vector3d::Zero());
        temp.setCharge(atom_charge);
        result.push_back(temp);
      }
    }
  }
  if (archive.empty()) {
    throw std::runtime_error(
        "Gaussian log file is missing the archive at the end.");
  }

  GetCoordinates(result, archive);

  if (!has_charges) {
    throw std::runtime_error("Charges not found in logfile.");
  }
  return result;
}

Eigen::Matrix3d Gaussian::GetPolarizability() const {

  if (!CheckLogFile()) {
    throw std::runtime_error("logfile not correctly formatted");
  }
  std::string line;
  ifstream input_file((_run_dir + "/" + _log_file_name));
  bool has_pol = false;

  std::vector<double> polar_coeff;

  Eigen::Matrix3d pol = Eigen::Matrix3d::Zero();
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);

    std::string::size_type pol_pos =
        line.find("Dipole polarizability, Alpha (input orientation)");
    if (pol_pos != std::string::npos) {
      XTP_LOG(logDEBUG, *_pLog) << "Getting polarizability" << flush;
      getline(input_file, line);
      getline(input_file, line);
      getline(input_file, line);
      getline(input_file, line);
      getline(input_file, line);
      for (int i = 0; i < 6; i++) {
        getline(input_file, line);
        tools::Tokenizer tok2(line, " ");
        std::vector<std::string> values = tok2.ToVector();
        if (values.size() != 4) {
          throw std::runtime_error("Polarisation line " + line +
                                   " cannot be parsed");
        }
        std::string value = values[1];
        boost::replace_first(value, "D", "e");
        polar_coeff.push_back(std::stod(value));  // xx xy yy zx zy zz
      }
      pol << polar_coeff[0], polar_coeff[1], polar_coeff[3], polar_coeff[1],
          polar_coeff[2], polar_coeff[4], polar_coeff[3], polar_coeff[4],
          polar_coeff[5];

      has_pol = true;
    }
  }
  if (!has_pol) {
    throw std::runtime_error("Could not find polarisation in logfile");
  }
  return pol;
}

void Gaussian::GetArchive(std::vector<std::string>& archive, std::string& line,
                          std::ifstream& input_file) const {
  std::string::size_type endseg_pos = line.find("\\");
  if (endseg_pos != std::string::npos) {
    boost::trim(line);
    std::string sum = line;
    while (line.size() != 0) {
      getline(input_file, line);
      boost::trim(line);
      sum += line;
    }
    std::vector<std::string> strings;
    boost::iter_split(strings, sum, boost::first_finder("\\\\"));
    archive = strings;
  }
}

double Gaussian::GetQMEnergy(const std::vector<std::string>& archive) const {
  double qm_energy = 0.0;
  tools::Tokenizer tok3(archive[4], "\\");
  std::vector<std::string> blocks = tok3.ToVector();
  map<std::string, std::string> properties;
  for (const std::string& block : blocks) {
    tools::Tokenizer tok4(block, "=");
    std::vector<std::string> property = tok4.ToVector();
    properties[property[0]] = property[1];
  }
  if (properties.count("HF") > 0) {
    qm_energy = boost::lexical_cast<double>(properties["HF"]);

    XTP_LOG(logDEBUG, *_pLog)
        << (boost::format("QM energy[Hrt]: %4.6f ") % qm_energy).str() << flush;
  } else {
    throw std::runtime_error("ERROR No energy in archive");
  }
  return qm_energy;
}

template <class T>
void Gaussian::GetCoordinates(T& mol,
                              const std::vector<std::string>& archive) const {
  typedef typename std::iterator_traits<typename T::iterator>::value_type Atom;
  XTP_LOG(logDEBUG, *_pLog) << "Getting the coordinates" << flush;
  bool has_atoms = mol.size() > 0;
  tools::Tokenizer tok(archive[3], "\\");
  std::vector<std::string> atom_block = tok.ToVector();
  for (unsigned i = 0; i < atom_block.size() - 1; i++) {
    tools::Tokenizer tok2(atom_block[i + 1], ",");
    std::vector<std::string> atom = tok2.ToVector();
    std::string atom_type = atom[0];
    int endindex = atom.size() - 1;
    double x = boost::lexical_cast<double>(atom[endindex - 2]);
    double y = boost::lexical_cast<double>(atom[endindex - 1]);
    double z = boost::lexical_cast<double>(atom[endindex]);
    Eigen::Vector3d pos(x, y, z);
    pos *= tools::conv::ang2bohr;
    if (has_atoms == false) {
      mol.push_back(Atom(i, atom_type, pos));
    } else {
      Atom& pAtom = mol.at(i);
      pAtom.setPos(pos);
    }
  }
}

/**
 * Parses the Gaussian Log file and stores data in the Orbitals object
 */
bool Gaussian::ParseLogFile(Orbitals& orbitals) {
  std::string line;
  std::vector<std::string> results;
  bool has_occupied_levels = false;
  bool has_unoccupied_levels = false;

  int occupied_levels = 0;
  int unoccupied_levels = 0;
  int number_of_electrons = 0;
  int basis_set_size = 0;

  XTP_LOG(logDEBUG, *_pLog) << "GAUSSIAN: parsing " << _log_file_name << flush;

  std::string log_file_name_full = _log_file_name;
  if (_run_dir != "") log_file_name_full = _run_dir + "/" + _log_file_name;

  // check if LOG file is complete
  if (!CheckLogFile()) return false;

  // save qmpackage name
  orbitals.setQMpackage(getPackageName());
  orbitals.setDFTbasisName(_basisset_name);

  if (_write_pseudopotentials) {
    orbitals.setECPName(_ecp_name);
  }

  double self_energy = 0.0;
  bool ScaHFX_found = false;
  // Start parsing the file line by line
  ifstream input_file(log_file_name_full);
  std::vector<std::string> archive;

  while (input_file) {

    getline(input_file, line);
    boost::trim(line);

    /* Check for ScaHFX = factor of HF exchange included in functional */
    std::string::size_type HFX_pos = line.find("ScaHFX=");
    if (HFX_pos != std::string::npos) {
      tools::Tokenizer tok(line, "\t ");
      std::vector<std::string> line_split = tok.ToVector();
      double ScaHFX = boost::lexical_cast<double>(line_split.back());
      orbitals.setScaHFX(ScaHFX);
      ScaHFX_found = true;
      XTP_LOG(logDEBUG, *_pLog)
          << "DFT with " << ScaHFX << " of HF exchange!" << flush;
    }

    /*
     * number of occupied and virtual orbitals
     * N alpha electrons      M beta electrons
     */
    std::string::size_type electrons_pos = line.find("alpha electrons");
    if (electrons_pos != std::string::npos) {
      tools::Tokenizer tok(line, "\t ");
      std::vector<std::string> line_split = tok.ToVector();
      number_of_electrons = boost::lexical_cast<int>(line_split.front());
      orbitals.setNumberOfAlphaElectrons(number_of_electrons);
      XTP_LOG(logDEBUG, *_pLog)
          << "Alpha electrons: " << number_of_electrons << flush;
    }

    /*
     * basis set size
     * N basis functions,  M primitive gaussians,   K cartesian basis functions
     */
    std::string::size_type basis_pos = line.find("basis functions,");
    if (basis_pos != std::string::npos) {
      tools::Tokenizer tok(line, "\t ");
      std::vector<std::string> line_split = tok.ToVector();
      basis_set_size = boost::lexical_cast<int>(line_split.front());
      orbitals.setBasisSetSize(basis_set_size);
      XTP_LOG(logDEBUG, *_pLog)
          << "Basis functions: " << basis_set_size << flush;
    }

    /*
     * energies of occupied/unoccupied levels
     * Alpha  occ.(virt.) eigenvalues -- e1 e2 e3 e4 e5
     */
    std::string::size_type eigenvalues_pos = line.find("Alpha");
    if (eigenvalues_pos != std::string::npos) {
      std::list<std::string> stringList;
      while (eigenvalues_pos != std::string::npos && !has_occupied_levels &&
             !has_unoccupied_levels) {

        boost::iter_split(stringList, line, boost::first_finder("--"));
        boost::trim(stringList.back());
        tools::Tokenizer tok(stringList.back(), "\t ");
        std::vector<std::string> energies = tok.ToVector();

        if (stringList.front().find("virt.") != std::string::npos) {
          unoccupied_levels += energies.size();
          energies.clear();
        }
        if (stringList.front().find("occ.") != std::string::npos) {
          occupied_levels += energies.size();
          energies.clear();
        }
        getline(input_file, line);
        eigenvalues_pos = line.find("Alpha");
        boost::trim(line);

        if (eigenvalues_pos == std::string::npos) {
          has_occupied_levels = true;
          has_unoccupied_levels = true;
          orbitals.setNumberOfOccupiedLevels(occupied_levels);
          XTP_LOG(logDEBUG, *_pLog)
              << "Occupied levels: " << occupied_levels << flush;
          XTP_LOG(logDEBUG, *_pLog)
              << "Unoccupied levels: " << unoccupied_levels << flush;
        }
      }  // end of the while loop
    }    // end of the eigenvalue parsing

    /*
     * Coordinates of the final configuration
     * stored in the archive at the end of the file
     */
    GetArchive(archive, line, input_file);

    std::string::size_type self_energy_pos =
        line.find("Self energy of the charges");

    if (self_energy_pos != std::string::npos) {
      XTP_LOG(logDEBUG, *_pLog) << "Getting the self energy\n";
      tools::Tokenizer tok(line, "=");
      std::vector<std::string> block = tok.ToVector();
      tools::Tokenizer tok2(block[1], "\t ");
      std::vector<std::string> energy = tok2.ToVector();
      self_energy = boost::lexical_cast<double>(energy[0]);
      XTP_LOG(logDEBUG, *_pLog) << "Self energy " << self_energy << flush;
    }

  }  // end of reading the file line-by-line

  if (archive.empty()) {
    throw std::runtime_error(
        "Gaussian log file is missing the archive at the end.");
  }

  QMMolecule& mol = orbitals.QMAtoms();
  GetCoordinates(mol, archive);
  double qm_energy = GetQMEnergy(archive);
  orbitals.setQMEnergy(qm_energy - self_energy);
  XTP_LOG(logDEBUG, *_pLog) << "Done parsing" << flush;
  input_file.close();

  if (!ScaHFX_found) {
    XTP_LOG(logDEBUG, *_pLog)
        << "WARNING === WARNING \n, could not find ScaHFX= entry in log."
           "\n probably you forgt #P in the beginning of the input file.\n"
           " If you are running a hybrid functional calculation redo it! Now! "
           "Please!\n ===WARNING=== \n"
        << flush;
    orbitals.setScaHFX(0.0);
  }

  return true;
}

std::string Gaussian::FortranFormat(double number) {
  std::stringstream ssnumber;
  if (number >= 0) ssnumber << " ";
  ssnumber << setiosflags(ios::fixed) << setprecision(8) << std::scientific
           << number;
  std::string snumber = ssnumber.str();
  boost::replace_first(snumber, "e", "D");
  return snumber;
}

}  // namespace xtp
}  // namespace votca
