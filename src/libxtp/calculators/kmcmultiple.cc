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

#include "kmcmultiple.h"
#include "votca/xtp/qmstate.h"
#include <boost/format.hpp>
#include <locale>
#include <votca/tools/constants.h>
#include <votca/tools/property.h>
#include <votca/xtp/gnode.h>
#include <votca/xtp/topology.h>

using namespace std;

namespace votca {
namespace xtp {
void KMCMultiple::Initialize(tools::Property& options) {
  std::string key = "options." + Identify();
  ParseCommonOptions(options);
  _runtime =
      options.ifExistsReturnElseThrowRuntimeError<double>(key + ".runtime");
  _field = options.ifExistsReturnElseReturnDefault<Eigen::Vector3d>(
      key + ".field", Eigen::Vector3d::Zero());
  double mtobohr = 1E9 * tools::conv::nm2bohr;
  _field *=
      (tools::conv::ev2hrt / mtobohr);  // Converting from V/m to Hartree/bohr

  _outputtime =
      options.ifExistsReturnElseReturnDefault<double>(key + ".outputtime", 0);
  _timefile = options.ifExistsReturnElseReturnDefault<std::string>(
      key + ".timefile", _timefile);

  std::string carriertype =
      options.ifExistsReturnElseReturnDefault<std::string>(key + ".carriertype",
                                                           "e");
  _carriertype = QMStateType(carriertype);
  if (!_carriertype.isKMCState()) {
    throw runtime_error("KMC cannot be run for state:" +
                        _carriertype.ToLongString());
  }

  return;
}

void KMCMultiple::PrintDiffandMu(const Eigen::Matrix3d& avgdiffusiontensor,
                                 double simtime, unsigned long step) const {
  double absolute_field = _field.norm();

  if (absolute_field == 0) {
    unsigned long diffusionsteps = step / _diffusionresolution;
    Eigen::Matrix3d result =
        avgdiffusiontensor / (diffusionsteps * 2 * simtime * _numberofcarriers);
    cout << endl
         << "Step: " << step
         << " Diffusion tensor averaged over all carriers (nm^2/s):" << endl
         << result * tools::conv::bohr2nm * tools::conv::bohr2nm << endl;
  } else {
    double average_mobility = 0;
    double bohr2Hrts_to_nm2Vs =
        tools::conv::bohr2nm * tools::conv::bohr2nm / tools::conv::hrt2ev;
    cout << endl << "Mobilities (nm^2/Vs): " << endl;
    for (int i = 0; i < _numberofcarriers; i++) {
      Eigen::Vector3d velocity = _carriers[i].get_dRtravelled() / simtime;
      double mobility =
          velocity.dot(_field) / (absolute_field * absolute_field);
      cout << std::scientific << "    carrier " << i + 1
           << ": mu=" << mobility * bohr2Hrts_to_nm2Vs << endl;
      average_mobility +=
          velocity.dot(_field) / (absolute_field * absolute_field);
    }
    average_mobility /= _numberofcarriers;
    cout << std::scientific
         << "  Overall average mobility in field direction <mu>="
         << average_mobility * bohr2Hrts_to_nm2Vs << " nm^2/Vs  " << endl;
  }
}

void KMCMultiple::WriteToTrajectory(std::fstream& traj,
                                    vector<Eigen::Vector3d>& startposition,
                                    double simtime, unsigned long step) const {
  traj << simtime << "\t";
  traj << step << "\t";
  for (int i = 0; i < _numberofcarriers; i++) {
    Eigen::Vector3d pos = startposition[i] + _carriers[i].get_dRtravelled();
    traj << pos.x() * tools::conv::bohr2nm << "\t";
    traj << pos.y() * tools::conv::bohr2nm << "\t";
    traj << pos.z() * tools::conv::bohr2nm;
    if (i < _numberofcarriers - 1) {
      traj << "\t";
    } else {
      traj << endl;
    }
  }
}

void KMCMultiple::WriteToEnergyFile(std::fstream& tfile, double simtime,
                                    unsigned long step) const {
  double absolute_field = _field.norm();
  double currentenergy = 0;
  double currentmobility = 0;
  Eigen::Vector3d dr_travelled_current = Eigen::Vector3d::Zero();
  double dr_travelled_field = 0.0;
  Eigen::Vector3d avgvelocity_current = Eigen::Vector3d::Zero();
  if (absolute_field != 0) {
    for (const auto& carrier : _carriers) {
      dr_travelled_current += carrier.get_dRtravelled();
      currentenergy += carrier.getCurrentEnergy();
    }
    dr_travelled_current /= _numberofcarriers;
    currentenergy /= _numberofcarriers;
    avgvelocity_current = dr_travelled_current / simtime;
    currentmobility =
        avgvelocity_current.dot(_field) / (absolute_field * absolute_field);
    dr_travelled_field = dr_travelled_current.dot(_field) / absolute_field;
  }
  double bohr2Hrts_to_nm2Vs =
      tools::conv::bohr2nm * tools::conv::bohr2nm / tools::conv::hrt2ev;
  tfile << simtime << "\t" << step << "\t"
        << currentenergy * tools::conv::hrt2ev << "\t"
        << currentmobility * bohr2Hrts_to_nm2Vs << "\t"
        << dr_travelled_field * tools::conv::bohr2nm << "\t"
        << dr_travelled_current.norm() * tools::conv::bohr2nm << "\t" << endl;
}

void KMCMultiple::PrintDiagDandMu(const Eigen::Matrix3d& avgdiffusiontensor,
                                  double simtime, unsigned long step) const {
  unsigned long diffusionsteps = step / _diffusionresolution;
  Eigen::Matrix3d result =
      avgdiffusiontensor / (diffusionsteps * 2 * simtime * _numberofcarriers);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
  es.computeDirect(result);
  double bohr2_nm2 = tools::conv::bohr2nm * tools::conv::bohr2nm;
  cout << endl << "Eigenvalues: " << endl << endl;
  for (int i = 0; i < 3; i++) {
    cout << "Eigenvalue: " << es.eigenvalues()(i) * bohr2_nm2 << endl
         << "Eigenvector: ";
    cout << es.eigenvectors().col(i).x() << "   ";
    cout << es.eigenvectors().col(i).y() << "   ";
    cout << es.eigenvectors().col(i).z() << endl << endl;
  }
  double bohr2Hrts_to_nm2Vs =
      tools::conv::bohr2nm * tools::conv::bohr2nm / tools::conv::hrt2ev;
  // calculate average mobility from the Einstein relation
  if (_field.norm() == 0) {
    cout << "The following value is calculated using the Einstein relation "
            "and assuming an isotropic medium"
         << endl;
    double avgD = 1. / 3. * es.eigenvalues().sum();
    double average_mobility = std::abs(avgD / _temperature);
    cout << std::scientific << "  Overall average mobility <mu>="
         << average_mobility * bohr2Hrts_to_nm2Vs << " nm^2/Vs " << endl;
  }
}

void KMCMultiple::PrintChargeVelocity(double simtime) const {
  Eigen::Vector3d avg_dr_travelled = Eigen::Vector3d::Zero();
  for (int i = 0; i < _numberofcarriers; i++) {
    cout << std::scientific << "    carrier " << i + 1 << ": "
         << _carriers[i].get_dRtravelled().transpose() / simtime *
                tools::conv::bohr2nm
         << endl;
    avg_dr_travelled += _carriers[i].get_dRtravelled();
  }
  avg_dr_travelled /= _numberofcarriers;

  Eigen::Vector3d avgvelocity = avg_dr_travelled / simtime;
  cout << std::scientific << "  Overall average velocity (nm/s): "
       << avgvelocity.transpose() * tools::conv::bohr2nm << endl;
}

void KMCMultiple::RunVSSM() {

  int realtime_start = time(NULL);
  cout << endl << "Algorithm: VSSM for Multiple Charges" << endl;
  cout << "number of carriers: " << _numberofcarriers << endl;
  cout << "number of nodes: " << _nodes.size() << endl;

  bool checkifoutput = (_outputtime != 0);
  double nexttrajoutput = 0;
  unsigned long maxsteps = _runtime;
  unsigned long outputstep = _outputtime;
  bool stopontime = false;

  if (_runtime > 100) {
    cout << "stop condition: " << maxsteps << " steps." << endl;

    if (checkifoutput) {
      cout << "output frequency: ";
      cout << "every " << outputstep << " steps." << endl;
    }
  } else {
    stopontime = true;
    cout << "stop condition: " << _runtime << " seconds runtime." << endl;

    if (checkifoutput) {
      cout << "output frequency: ";
      cout << "every " << _outputtime << " seconds." << endl;
    }
  }
  cout << "(If you specify runtimes larger than 100 kmcmultiple assumes that "
          "you are specifying the number of steps for both runtime and "
          "outputtime.)"
       << endl;

  if (!stopontime && _outputtime != 0 && floor(_outputtime) != _outputtime) {
    throw runtime_error(
        "ERROR in kmcmultiple: runtime was specified in steps (>100) and "
        "outputtime in seconds (not an integer). Please use the same units for "
        "both input parameters.");
  }

  if (_numberofcarriers > int(_nodes.size())) {
    throw runtime_error(
        "ERROR in kmcmultiple: specified number of carriers is greater than "
        "the "
        "number of nodes. This conflicts with single occupation.");
  }

  std::fstream traj;
  std::fstream tfile;

  if (checkifoutput) {

    cout << "Writing trajectory to " << _trajectoryfile << "." << endl;
    traj.open(_trajectoryfile, std::fstream::out);

    traj << "time[s]\tsteps\t";
    for (int i = 0; i < _numberofcarriers; i++) {
      traj << "carrier" << i + 1 << "_x\t";
      traj << "carrier" << i + 1 << "_y\t";
      traj << "carrier" << i + 1 << "_z";
      if (i < _numberofcarriers - 1) {
        traj << "'\t";
      }
    }
    traj << endl;
    if (!_timefile.empty()) {
      cout << "Writing time dependence of energy and mobility to " << _timefile
           << "." << endl;
      tfile.open(_timefile, std::fstream::out);
      tfile << "time[s]\t "
               "steps\tenergy_per_carrier[eV]\tmobility[nm**2/"
               "Vs]\tdistance_fielddirection[nm]\tdistance_absolute[nm]"
            << endl;
    }
  }
  RandomlyCreateCharges();
  vector<Eigen::Vector3d> startposition(_numberofcarriers,
                                        Eigen::Vector3d::Zero());
  for (int i = 0; i < _numberofcarriers; i++) {
    startposition[i] = _carriers[i].getCurrentPosition();
  }

  traj << 0 << "\t";
  traj << 0 << "\t";
  for (int i = 0; i < _numberofcarriers; i++) {
    traj << startposition[i].x() * tools::conv::bohr2nm << "\t";
    traj << startposition[i].y() * tools::conv::bohr2nm << "\t";
    traj << startposition[i].z() * tools::conv::bohr2nm;
    if (i < _numberofcarriers - 1) {
      traj << "\t";
    } else {
      traj << endl;
    }
  }

  vector<GNode*> forbiddennodes;
  vector<GNode*> forbiddendests;

  Eigen::Matrix3d avgdiffusiontensor = Eigen::Matrix3d::Zero();

  double simtime = 0.0;
  unsigned long step = 0;

  while (((stopontime && simtime < _runtime) ||
          (!stopontime && step < maxsteps))) {

    if ((time(NULL) - realtime_start) > _maxrealtime * 60. * 60.) {
      cout << endl
           << "Real time limit of " << _maxrealtime << " hours ("
           << int(_maxrealtime * 60 * 60 + 0.5)
           << " seconds) has been reached. Stopping here." << endl
           << endl;
      break;
    }

    double cumulated_rate = 0;
    for (const auto& carrier : _carriers) {
      cumulated_rate += carrier.getCurrentEscapeRate();
    }
    if (cumulated_rate <= 0) {  // this should not happen: no possible jumps
                                // defined for a node
      throw runtime_error(
          "ERROR in kmcmultiple: Incorrect rates. All "
          "the escape rates for the current setting are 0.");
    }

    double dt = Promotetime(cumulated_rate);

    simtime += dt;
    step++;

    for (auto& carrier : _carriers) {
      carrier.updateOccupationtime(dt);
    }

    ResetForbiddenlist(forbiddennodes);
    bool level1step = true;
    while (level1step) {

      // determine which electron will escape
      GNode* newnode = NULL;
      Chargecarrier* affectedcarrier = ChooseAffectedCarrier(cumulated_rate);

      if (CheckForbidden(affectedcarrier->getCurrentNode(), forbiddennodes)) {
        continue;
      }
      ResetForbiddenlist(forbiddendests);
      while (true) {
        // LEVEL 2

        const GLink& event =
            ChooseHoppingDest(affectedcarrier->getCurrentNode());
        newnode = event.getDestination();

        if (newnode == NULL) {
          AddtoForbiddenlist(affectedcarrier->getCurrentNode(), forbiddennodes);
          break;  // select new escape node (ends level 2 but without setting
                  // level1step to 1)
        }

        // check after the event if this was allowed
        if (CheckForbidden(*newnode, forbiddendests)) {
          continue;
        }

        // if the new segment is unoccupied: jump; if not: add to forbidden
        // list and choose new hopping destination
        if (newnode->isOccupied()) {
          if (CheckSurrounded(affectedcarrier->getCurrentNode(),
                              forbiddendests)) {
            AddtoForbiddenlist(affectedcarrier->getCurrentNode(),
                               forbiddennodes);
            break;  // select new escape node (ends level 2 but without
                    // setting level1step to 1)
          }
          AddtoForbiddenlist(*newnode, forbiddendests);
          continue;  // select new destination
        } else {
          affectedcarrier->jumpAccordingEvent(event);
          level1step = false;
          break;  // this ends LEVEL 2 , so that the time is updated and the
                  // next MC step started
        }
        // END LEVEL 2
      }
      // END LEVEL 1
    }

    if (step % _diffusionresolution == 0) {
      for (const auto& carrier : _carriers) {
        avgdiffusiontensor += (carrier.get_dRtravelled()) *
                              (carrier.get_dRtravelled()).transpose();
      }
    }

    if (step != 0 && step % _intermediateoutput_frequency == 0) {
      PrintDiffandMu(avgdiffusiontensor, simtime, step);
    }

    if (checkifoutput) {
      bool outputsteps = (!stopontime && step % outputstep == 0);
      bool outputtime = (stopontime && simtime > nexttrajoutput);
      if (outputsteps || outputtime) {
        // write to trajectory file
        nexttrajoutput = simtime + _outputtime;
        WriteToTrajectory(traj, startposition, simtime, step);
        if (!_timefile.empty()) {
          WriteToEnergyFile(tfile, simtime, step);
        }
      }
    }
  }  // KMC

  if (checkifoutput) {
    traj.close();
    if (!_timefile.empty()) {
      tfile.close();
    }
  }

  WriteOccupationtoFile(simtime, _occfile);

  cout << endl << "finished KMC simulation after " << step << " steps." << endl;
  cout << "simulated time " << simtime << " seconds." << endl;
  cout << "runtime: ";
  cout << endl << endl;

  PrintChargeVelocity(simtime);

  cout << endl << "Distances travelled (nm): " << endl;
  for (int i = 0; i < _numberofcarriers; i++) {
    cout << std::scientific << "    carrier " << i + 1 << ": "
         << _carriers[i].get_dRtravelled().transpose() * tools::conv::bohr2nm
         << endl;
  }

  PrintDiffandMu(avgdiffusiontensor, simtime, step);
  PrintDiagDandMu(avgdiffusiontensor, simtime, step);

  return;
}

bool KMCMultiple::EvaluateFrame(Topology& top) {
  std::cout << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  std::cout << "      KMC FOR MULTIPLE CHARGES" << std::endl;
  std::cout << "-----------------------------------" << std::endl << std::endl;

  // Initialise random number generator
  if (tools::globals::verbose) {
    cout << endl << "Initialising random number generator" << endl;
  }
  srand(_seed);  // srand expects any integer in order to initialise the random
                 // number generator
  _RandomVariable = tools::Random2();
  _RandomVariable.init(rand(), rand(), rand(), rand());

  LoadGraph(top);
  RunVSSM();

  return true;
}

}  // namespace xtp
}  // namespace votca
