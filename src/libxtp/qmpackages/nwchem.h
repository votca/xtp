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
#ifndef __VOTCA_XTP_NWCHEM_H
#define __VOTCA_XTP_NWCHEM_H

#include <votca/xtp/qmpackage.h>

namespace votca {
namespace xtp {
/**
    \brief Wrapper for the Gaussian program

    The Gaussian class executes the Gaussian package
    and extracts information from its log and io files

*/
class Orbitals;
class NWChem : public QMPackage {
 public:
  std::string getPackageName() const { return "nwchem"; }

  void Initialize(tools::Property& options);

  bool WriteInputFile(const Orbitals& orbitals);

  bool Run();

  void CleanUp();

  bool ParseLogFile(Orbitals& orbitals);

  bool ParseMOsFile(Orbitals& orbitals);

  StaticSegment GetCharges() const;

  Eigen::Matrix3d GetPolarizability() const;

 private:
  std::string ascii_mo_file_name;
  bool CheckLogFile() const;
  bool WriteShellScript();
  bool WriteGuess(const Orbitals& orbitals);

  void WriteBasisset(std::ofstream& nw_file, const QMMolecule& qmatoms);
  void WriteECP(std::ofstream& nw_file, const QMMolecule& qmatoms);

  std::string FortranFormat(const double& number);
  int WriteBackgroundCharges(std::ofstream& nw_file);
  void WriteChargeOption();
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_NWCHEM_H */
