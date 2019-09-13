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

#pragma once
#ifndef VOTCA_XTP_GWBSEENGINE_H
#define VOTCA_XTP_GWBSEENGINE_H

#include <votca/tools/property.h>
#include <votca/xtp/logger.h>

namespace votca {
namespace xtp {
class QMPackage;
class Orbitals;

/**
 * \brief Electronic Excitations via Density-Functional Theory
 *
 * Evaluates electronic ground state in molecular systems based on
 * density functional theory with Gaussian Orbitals.
 *
 */

class GWBSEEngine {
 public:
  std::string Identify() { return "gwbse_engine"; }

  void Initialize(tools::Property& options, std::string archive_filename);
  void ExcitationEnergies(Orbitals& orbitals);

  void setLog(Logger* pLog) { _pLog = pLog; }

  void setQMPackage(QMPackage* qmpackage) { _qmpackage = qmpackage; }

  std::string GetDFTLog() const { return _dftlog_file; };

  void setLoggerFile(std::string logger_file) { _logger_file = logger_file; };

  void setRedirectLogger(bool redirect_logger) {
    _redirect_logger = redirect_logger;
  };

  const tools::Property& ReportSummary() const { return _summary; };

 private:
  QMPackage* _qmpackage;

  Logger* _pLog;

  // task options
  bool _do_guess = false;
  bool _do_dft_input = false;
  bool _do_dft_run = false;
  bool _do_dft_parse = false;
  bool _do_gwbse = false;
  bool _redirect_logger = false;

  // DFT log and MO file names
  std::string _MO_file;      // file containing the MOs from qmpackage...
  std::string _dftlog_file;  // file containing the Energies etc... from
                             // qmpackage...
  std::string _logger_file;
  std::string _archive_file;
  std::string _guess_archiveA;
  std::string _guess_archiveB;

  // Options for GWBSE module
  tools::Property _gwbse_options;
  tools::Property _summary;

  void WriteLoggerToFile(Logger* pLog);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GWBSEENGINE_H
