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
  bool _sigma_export_converged = false;
  double _sigma_spectral_eta = 1e-4;
  bool _rpa_spectrum_export = false;

 public:
  static void Load();

  static bool Hedin() { return _instance._hedin; }
  static bool GSCExport() { return _instance._gsc_export; }
  static double GSCAlpha() { return _instance._gsc_alpha; }
  static int SigmaExportRange() { return _instance._sigma_export_range; }
  static double SigmaExportDelta() { return _instance._sigma_export_delta; }
  static bool SigmaExportConverged() { return _instance._sigma_export_converged; }
  static double SigmaSpectralEta() { return _instance._sigma_spectral_eta; }
  static double RPASpectrumExport() { return _instance._rpa_spectrum_export; }
};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_CUSTOM_TOOLS_H */
