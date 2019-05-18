#ifndef _VOTCA_XTP_CUSTOM_OPTS_H
#define _VOTCA_XTP_CUSTOM_OPTS_H

#include <votca/tools/property.h>
#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {

class CustomTools {
 public:
  static void ExportMat(std::string filename, const Eigen::MatrixXd& mat);
  static void ExportVec(std::string filename, const Eigen::VectorXd& vec);
  static void AppendRow(std::string filename, const Eigen::VectorXd& row);
  static void ExportMatBinary(std::string filename, const Eigen::MatrixXd& mat);
};

class CustomOpts {
 private:
  static CustomOpts _instance;
  CustomOpts() {}
  void Parse(tools::Property& options);
  void Report();

  // Sigma exact options
  double _sigma_spectral_eta = 1e-1;
  bool   _COHSEX             = false;
  // G/GW SC cycle options
  bool   _gsc_export = false;
  double _gsc_alpha  = 0.0;
  // RPA spectrum export
  bool _rpa_spectrum_export = false;
  // Sigma_c diagonal export
  int    _sigma_export_range     = 0;
  double _sigma_export_delta     = 1.0;
  bool   _sigma_export_converged = false;
  bool   _sigma_export_binary    = false;
  // Sigma_c matrix export
  bool _sigma_matrix_export = false;

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
};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_CUSTOM_TOOLS_H */
