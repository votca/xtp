#ifndef _VOTCA_XTP_CUSTOM_OPTS_H
#define _VOTCA_XTP_CUSTOM_OPTS_H

#include <votca/tools/property.h>
#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {

class CustomTools {
 public:
  static void ExportMat(std::string filename, const Eigen::MatrixXd& mat);
  static void ExportMatBinary(std::string filename, const Eigen::MatrixXd& mat);
  static Eigen::MatrixXd ImportMatBinary(std::string filename);
};

class CustomOpts {
 private:
  static CustomOpts _instance;
  CustomOpts() {}

  // Sigma exact options
  double _sigma_spectral_eta = 1e-1;
  bool   _COHSEX             = false;
  // G/GW SC cycle options
  bool   _gsc_export = false;
  double _gsc_alpha  = 0.0;
  // RPA spectrum export
  bool _rpa_spectrum_export = false; // TODO
  // Sigma_c diagonal export
  int    _sigma_export_range     = 0;
  double _sigma_export_delta     = 1.0;
  bool   _sigma_export_converged = false;
  bool   _sigma_export_binary    = false;
  // Sigma_c matrix export
  bool _sigma_matrix_export = false;
  // GW DFT Shift
  double _gw_dft_shift = 0.0;
  // RPA ENergies
  bool _rpa_energies_import = false;
  bool _rpa_energies_export = false;
  
  void Parse(tools::Property& options);
  void Report();

 public:
  static void Load();
  
  static double SigmaSpectralEta() { return _instance._sigma_spectral_eta; }
  static bool COHSEX() { return _instance._COHSEX; }
  static bool GSCExport() { return _instance._gsc_export; }
  static double GSCAlpha() { return _instance._gsc_alpha; }
  static bool RPASpectrumExport() { return _instance._rpa_spectrum_export; }
  static int SigmaExportRange() { return _instance._sigma_export_range; }
  static double SigmaExportDelta() { return _instance._sigma_export_delta; }
  static bool SigmaExportConverged() { return _instance._sigma_export_converged; }
  static bool SigmaExportBinary() { return _instance._sigma_export_binary; }
  static bool SigmaMatrixExport() { return _instance._sigma_matrix_export; }
  static double GWDFTShift() { return _instance._gw_dft_shift; }
  static bool RPAEnergiesImport() { return _instance._rpa_energies_import; }
  static bool RPAEnergiesExport() { return _instance._rpa_energies_export; }
};

class GWSelfConsistencyLogger {
 private:
  static GWSelfConsistencyLogger _instance;
  GWSelfConsistencyLogger() {}
  int _qp_total;
  int _g_max;
  int _gw_iter;
  int _g_iter;
  Eigen::MatrixXd _log;

  void SelfInitialize(int qp_total, int gw_max);
  void SelfLogFrequencies(const Eigen::VectorXd& frequencies);
  void SelfWriteGWIter(bool conv);
  void SelfWriteCount(bool conv);
  
 public:
  static void Initialize(int qp_total, int g_max) {
    _instance.SelfInitialize(qp_total, g_max);
  }
  static void LogFrequencies(const Eigen::VectorXd& frequencies) {
    _instance.SelfLogFrequencies(frequencies);
  }
  static void WriteGWIter(bool conv) {
    _instance.SelfWriteGWIter(conv);
  }
  static void WriteCount(bool conv) {
    _instance.SelfWriteCount(conv);
  }
};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_CUSTOM_TOOLS_H */
