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

#include <fstream>
#include <iostream>
#include <string.h>
#include <votca/xtp/customtools.h>
namespace votca {
namespace xtp {

/* Custom Tools */

void CustomTools::ExportMat(std::string filename, const Eigen::MatrixXd& mat) {
  Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", "\n", "", "");
  std::ofstream file;
  file.open(filename, std::ios_base::trunc);
  file << mat.format(fmt) << std::endl;
  file.close();
}

void CustomTools::ExportVec(std::string filename, const Eigen::VectorXd& vec) {
  Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", "\n", "", "");
  std::ofstream file;
  file.open(filename, std::ios_base::trunc);
  file << vec.format(fmt) << std::endl;
  file.close();
}

void CustomTools::AppendRow(std::string filename, const Eigen::VectorXd& row) {
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ",
                      "", "", "", "");
  std::ofstream file;
  file.open(filename, std::ios_base::app);
  file << row.format(fmt) << std::endl;
  file.close();
}

/* Custom Options */

CustomOpts CustomOpts::_instance;

void CustomOpts::Load() {
  tools::Property options;
  try {
    options.LoadFromXML("customopts.xml");
  } catch (...) {
    return;
  }
  CustomOpts::_instance.Parse(options);
  std::cout << std::endl << "************************************************";
  std::cout << std::endl << "Loaded custom options";
  CustomOpts::_instance.Report();
  std::cout << std::endl << "************************************************";
}

void CustomOpts::Parse(tools::Property& options) {
  _hedin =
      options.ifExistsReturnElseReturnDefault<bool>("customopts.hedin", _hedin);
  _gsc_export = options.ifExistsReturnElseReturnDefault<bool>(
      "customopts.gsc_export", _gsc_export);
  _gsc_alpha = options.ifExistsReturnElseReturnDefault<double>(
      "customopts.gsc_alpha", _gsc_alpha);
  _sigma_export_range = options.ifExistsReturnElseReturnDefault<int>(
      "customopts.sigma_export_range", _sigma_export_range);
  _sigma_export_delta = options.ifExistsReturnElseReturnDefault<double>(
      "customopts.sigma_export_delta", _sigma_export_delta);
  _sigma_spectral_eta = options.ifExistsReturnElseReturnDefault<double>(
      "customopts.sigma_spectral_eta", _sigma_spectral_eta);
  _rpa_spectrum_export = options.ifExistsReturnElseReturnDefault<bool>(
      "customopts.rpa_spectrum_export", _rpa_spectrum_export);
}

void CustomOpts::Report() {
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ",
                      "", "", "", "");
  std::cout << std::endl << "Hedin:        " << _hedin;
  std::cout << std::endl << "GSC Export:   " << _gsc_export;
  std::cout << std::endl << "GSC Alpha:    " << _gsc_alpha;
  std::cout << std::endl
            << "Sigma export: " << _sigma_export_range << ", "
            << _sigma_export_delta;
  std::cout << std::endl << "Sigma Eta:    " << _sigma_spectral_eta;
  std::cout << std::endl << "Spec. export: " << _rpa_spectrum_export;
}

}  // namespace xtp
}  // namespace votca
