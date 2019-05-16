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

void CustomTools::ExportMatBinary(std::string filename, const Eigen::MatrixXd& mat) {
  std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
  Eigen::MatrixXd::Index rows = mat.rows(), cols = mat.cols();
  out.write((char*) (&rows), sizeof(Eigen::MatrixXd::Index));
  out.write((char*) (&cols), sizeof(Eigen::MatrixXd::Index));
  out.write((char*) mat.data(), rows * cols * sizeof(Eigen::MatrixXd::Scalar));
  out.close();
}

/* Custom Options */

CustomOpts CustomOpts::_instance;

void CustomOpts::Load() {
  tools::Property options;
  try {
    load_property_from_xml(options, "customopts.xml");
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
  _sigma_export_converged = options.ifExistsReturnElseReturnDefault<bool>(
      "customopts.sigma_export_converged", _sigma_export_converged);
  _sigma_export_binary = options.ifExistsReturnElseReturnDefault<bool>(
      "customopts.sigma_export_binary", _sigma_export_binary);
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
            << "Sigma export: "
            << "range: "  << _sigma_export_range << ", "
            << "delta: "  << _sigma_export_delta << ", "
            << "binary: " << _sigma_export_binary;
  std::cout << std::endl << "Sigma Eta:    " << _sigma_spectral_eta;
  std::cout << std::endl << "Spec. export: " << _rpa_spectrum_export;
}

}  // namespace xtp
}  // namespace votca
