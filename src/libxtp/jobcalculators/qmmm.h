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
#ifndef VOTCA_XTP_QMMM_H
#define VOTCA_XTP_QMMM_H

#include <votca/xtp/parallelxjobcalc.h>

namespace votca {
namespace xtp {

/**
 * \brief QM/MM with different regions around
 *
 * Calculates properties of different regions inside a multiregion calculation
 *
 * Callname: qmmm
 */

class QMMM : public ParallelXJobCalc<std::vector<Job> > {
 public:
  void Initialize(tools::Property& options) override;
  std::string Identify() override { return "qmmm"; }
  Job::JobResult EvalJob(const Topology& top, Job& job,
                         QMThread& Thread) override;
  void WriteJobFile(const Topology& top) override;
  void ReadJobFile(Topology& top) override;

 private:
  bool hasQMRegion() const;
  tools::Property _regions_def;

  Index _max_iterations = 100;
  bool _print_regions_pdb = false;

  bool _write_parse = false;
  std::vector<QMState> _states;
};

}  // namespace xtp
}  // namespace votca
#endif
