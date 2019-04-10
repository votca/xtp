/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_CUSTOM_OPTS_H
#define _VOTCA_XTP_CUSTOM_OPTS_H

#include <votca/tools/property.h>
#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {

class CustomTools {
 public:
  // TODO: How could I merge ExportMat and ExportVec?
  static void ExportMat(std::string filename, const Eigen::MatrixXd& mat);
  static void ExportVec(std::string filename, const Eigen::VectorXd& vec);
  static void AppendRow(std::string filename, const Eigen::VectorXd& row);
};

class CustomOpts {
 private:
  static CustomOpts _instance;
  CustomOpts() {}
  void Parse(tools::Property& options);
  void Report();

  bool _hedin = false;
  bool _gsc_export = false;
  double _gsc_alpha = 0.0;
  int _sigma_export_range = 0;
  double _sigma_export_delta = 1.0;
  double _sigma_spectral_eta = 1e-4;
  bool _rpa_spectrum_export = false;

 public:
  static void Load();

  static bool Hedin() { return _instance._hedin; }
  static bool GSCExport() { return _instance._gsc_export; }
  static double GSCAlpha() { return _instance._gsc_alpha; }
  static int SigmaExportRange() { return _instance._sigma_export_range; }
  static double SigmaExportDelta() { return _instance._sigma_export_delta; }
  static double SigmaSpectralEta() { return _instance._sigma_spectral_eta; }
  static double RPASpectrumExport() { return _instance._rpa_spectrum_export; }
};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_CUSTOM_TOOLS_H */
