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

void CustomTools::ExportMatBinary(std::string filename, const Eigen::MatrixXd& mat) {
  std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
  Eigen::MatrixXd::Index rows = mat.rows(), cols = mat.cols();
  out.write((char*) (&rows), sizeof(Eigen::MatrixXd::Index));
  out.write((char*) (&cols), sizeof(Eigen::MatrixXd::Index));
  out.write((char*) mat.data(), rows * cols * sizeof(Eigen::MatrixXd::Scalar));
  out.close();
}

Eigen::MatrixXd CustomTools::ImportMatBinary(std::string filename) {
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  Eigen::MatrixXd::Index rows = 0, cols = 0;
  in.read((char*) (&rows), sizeof(typename Eigen::MatrixXd::Index));
  in.read((char*) (&cols), sizeof(typename Eigen::MatrixXd::Index));
  Eigen::MatrixXd mat;
  mat.resize(rows, cols);
  in.read((char*) mat.data(), rows * cols * sizeof(Eigen::MatrixXd::Scalar));
  in.close();
  return mat;
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
  _sigma_spectral_eta = options.ifExistsReturnElseReturnDefault<double>(
      "customopts.sigma_spectral_eta", _sigma_spectral_eta);
  _gsc_export = options.ifExistsReturnElseReturnDefault<bool>(
      "customopts.gsc_export", _gsc_export);
  _gsc_alpha = options.ifExistsReturnElseReturnDefault<double>(
      "customopts.gsc_alpha", _gsc_alpha);
  _rpa_spectrum_export = options.ifExistsReturnElseReturnDefault<bool>(
      "customopts.rpa_spectrum_export", _rpa_spectrum_export);
  _sigma_export_range = options.ifExistsReturnElseReturnDefault<int>(
      "customopts.sigma_export_range", _sigma_export_range);
  _sigma_export_delta = options.ifExistsReturnElseReturnDefault<double>(
      "customopts.sigma_export_delta", _sigma_export_delta);
  _sigma_export_converged = options.ifExistsReturnElseReturnDefault<bool>(
      "customopts.sigma_export_converged", _sigma_export_converged);
  _sigma_export_binary = options.ifExistsReturnElseReturnDefault<bool>(
      "customopts.sigma_export_binary", _sigma_export_binary);
  _sigma_matrix_export = options.ifExistsReturnElseReturnDefault<bool>(
      "customopts.sigma_matrix_export", _sigma_matrix_export);
  _gw_dft_shift = options.ifExistsReturnElseReturnDefault<double>(
      "customopts.gw_dft_shift", _gw_dft_shift);
  _rpa_energies_import = options.ifExistsReturnElseReturnDefault<bool>(
      "customopts.rpa_energies_import", _rpa_energies_import);
  _rpa_energies_export = options.ifExistsReturnElseReturnDefault<bool>(
      "customopts.rpa_energies_export", _rpa_energies_export);
}

void CustomOpts::Report() {
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ",
                      "", "", "", "");
  std::cout << std::endl << "Sigma Eta: " << _sigma_spectral_eta;
  std::cout << std::endl << "GSC Export: " << _gsc_export;
  std::cout << std::endl << "GSC Alpha: " << _gsc_alpha;
  std::cout << std::endl << "RPA Spect. export: " << _rpa_spectrum_export;
  std::cout << std::endl << "Sigma diag. export: "
            << "range: "  << _sigma_export_range << ", "
            << "delta: "  << _sigma_export_delta << ", "
            << "binary: " << _sigma_export_binary;
  std::cout << std::endl << "Sigm. mat. export: " << _sigma_matrix_export;
  std::cout << std::endl << "GW DFT shift: " << _gw_dft_shift;
  std::cout << std::endl << "RPA energies import: " << _rpa_energies_import;
  std::cout << std::endl << "RPA energies export: " << _rpa_energies_export;
}

/* GW Self Consistency Logger */

GWSelfConsistencyLogger GWSelfConsistencyLogger::_instance;

void GWSelfConsistencyLogger::SelfInitialize(int qp_total, int g_max) {
  _qp_total = qp_total;
  _g_max = g_max;
  _gw_iter = 0;
  _g_iter = 0;
  _log = Eigen::MatrixXd::Zero(_qp_total, _g_max);
  // BEGIN Create file
  int zero = 0;
  std::ofstream out("gwsc.bin", std::ios::out | std::ios::binary | std::ios::trunc);
  out.write((char*) (&zero), sizeof(int)); // version
  out.write((char*) (&zero), sizeof(int)); // gw iters
  out.write((char*) (&zero), sizeof(int)); // conv
  out.close();
  // END Create file
}

void GWSelfConsistencyLogger::SelfLogFrequencies(const Eigen::VectorXd& frequencies) {
  if (_g_iter < _g_max) { _log.col(_g_iter++) = frequencies; }
}

void GWSelfConsistencyLogger::SelfWriteGWIter(bool conv) {
  if (_g_iter == 0) { return; }
  // BEGIN Write frequencies
  int iconv = conv ? 1 : 0;
  Eigen::MatrixXd mat = _log.leftCols(_g_iter).transpose();
  std::ofstream out("gwsc.bin", std::ios::out | std::ios::binary | std::ios::app);
  Eigen::MatrixXd::Index rows = mat.rows(), cols = mat.cols();
  out.write((char*) (&iconv), sizeof(int));
  out.write((char*) (&rows), sizeof(Eigen::MatrixXd::Index));
  out.write((char*) (&cols), sizeof(Eigen::MatrixXd::Index));
  out.write((char*) mat.data(), rows * cols * sizeof(Eigen::MatrixXd::Scalar));
  out.close();
  // END Write frequencies
  _gw_iter++;
  _g_iter = 0;
  _log = Eigen::MatrixXd::Zero(_qp_total, _g_max);
}

void GWSelfConsistencyLogger::SelfWriteCount(bool conv) {
  // BEGIN Write count
  int iconv = conv ? 1 : 0;
  std::ofstream out("gwsc.bin", std::ios::in | std::ios::out | std::ios::binary);
  out.seekp(sizeof(int), std::ios::beg);
  out.write((char*) (&_gw_iter), sizeof(int));
  out.write((char*) (&iconv), sizeof(int));
  out.close();
  // END Write count
}

}  // namespace xtp
}  // namespace votca
