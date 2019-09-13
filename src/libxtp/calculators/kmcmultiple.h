/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef __VOTCA_KMC_MULTIPLE_H
#define __VOTCA_KMC_MULTIPLE_H

#include <fstream>
#include <votca/xtp/kmccalculator.h>
namespace votca {
namespace xtp {

class KMCMultiple : public KMCCalculator {
 public:
  KMCMultiple(){};
  ~KMCMultiple(){};
  bool WriteToStateFile() const { return false; }
  std::string Identify() { return "kmcmultiple"; }
  void Initialize(tools::Property& options);
  bool EvaluateFrame(Topology& top);

 private:
  void RunVSSM();
  void PrintChargeVelocity(double simtime) const;

  void PrintDiagDandMu(const Eigen::Matrix3d& avgdiffusiontensor,
                       double simtime, unsigned long step) const;

  void WriteToEnergyFile(std::fstream& tfile, double simtime,
                         unsigned long step) const;

  void WriteToTrajectory(std::fstream& traj,
                         std::vector<Eigen::Vector3d>& startposition,
                         double simtime, unsigned long step) const;

  void PrintDiffandMu(const Eigen::Matrix3d& avgdiffusiontensor, double simtime,
                      unsigned long step) const;

  double _runtime;
  double _outputtime;
  std::string _timefile = "";
  int _intermediateoutput_frequency = 10000;
  unsigned long _diffusionresolution = 1000;
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_KMC_MULTIPLE_H */
