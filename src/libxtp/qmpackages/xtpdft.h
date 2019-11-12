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
#ifndef __VOTCA_XTP_XTPDFT_H
#define __VOTCA_XTP_XTPDFT_H

#include <string>
#include <votca/xtp/dftengine.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/polarsite.h>
#include <votca/xtp/qmpackage.h>

namespace votca {
namespace xtp {

/**
    \brief Wrapper for the internal XTP DFT engine


 */

class XTPDFT : public QMPackage {
 public:
  std::string getPackageName() const override { return "xtp"; }

  void Initialize(tools::Property& options) override;

  bool WriteInputFile(const Orbitals& orbitals) override;

  bool Run() override;

  void CleanUp() override;

  bool CheckLogFile();

  bool ParseLogFile(Orbitals& orbitals) override;

  bool ParseMOsFile(Orbitals& orbitals) override;

  StaticSegment GetCharges() const override {
    throw std::runtime_error(
        "If you want partial charges just run the 'partialcharges' calculator");
  }

  Eigen::Matrix3d GetPolarizability() const override {
    throw std::runtime_error(
        "GetPolarizability() is not implemented for xtpdft");
  }

 private:
  void WriteChargeOption() override { return; }
  tools::Property _xtpdft_options;

  Orbitals _orbitals;
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_XTPDFT_H */
