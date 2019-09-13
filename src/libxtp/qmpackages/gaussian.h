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
#ifndef __VOTCA_XTP_GAUSSIAN_H
#define __VOTCA_XTP_GAUSSIAN_H

#include <votca/xtp/qmpackage.h>

namespace votca {
namespace xtp {
/**
    \brief Wrapper for the Gaussian program

    The Gaussian class executes the Gaussian package
    and extracts information from its log and io files

*/

class Orbitals;
class Gaussian : public QMPackage {
 public:
  std::string getPackageName() const { return "gaussian"; }

  void Initialize(tools::Property& options);

  bool WriteInputFile(const Orbitals& orbitals);

  bool Run();

  void CleanUp();

  bool ParseLogFile(Orbitals& orbitals);

  bool ParseMOsFile(Orbitals& orbitals);

  StaticSegment GetCharges() const;

  Eigen::Matrix3d GetPolarizability() const;

 private:
  bool WriteShellScript();

  bool CheckLogFile() const;

  std::string _vdWfooter = "";

  std::string FortranFormat(double number);
  void GetArchive(std::vector<std::string>& archive, std::string& line,
                  std::ifstream& input_file) const;

  template <class T>
  void GetCoordinates(T& mol, const std::vector<std::string>& archive) const;

  double GetQMEnergy(const std::vector<std::string>& archive) const;

  void WriteBasisset(std::ofstream& com_file, const QMMolecule& qmatoms);
  void WriteECP(std::ofstream& com_file, const QMMolecule& qmatoms);
  void WriteBackgroundCharges(std::ofstream& com_file);
  void WriteGuess(const Orbitals& orbitals_guess, std::ofstream& com_file);
  void WriteCoordinates(std::ofstream& com_file, const QMMolecule& qmatoms);
  void WriteHeader(std::ofstream& com_file);

  void WriteChargeOption();
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_GAUSSAIN_H */
