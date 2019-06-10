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

  // GW SC cycle export
  bool _gw_sc_export = false;
  // RPA spectrum export
  bool _rpa_spectrum_export = false; // TODO
  // Sigma_c diagonal export
  int    _sigma_export_range     = 0;
  double _sigma_export_delta     = 1.0;
  bool   _sigma_export_converged = false;
  // Sigma_c matrix export
  bool _sigma_matrix_export = false;
  // RPA Energies import/export
  bool _rpa_energies_import = false;
  bool _rpa_energies_export = false;
  // Misc.
  bool _export_binary = false;
  
  void Parse(tools::Property& options);
  void Report();

 public:
  static void Load();
  
  // GW SC cycle export
  static bool GWSCExport() { return _instance._gw_sc_export; }
  // RPA spectrum export
  static bool RPASpectrumExport() { return _instance._rpa_spectrum_export; }
  // Sigma_c diagonal export
  static int    SigmaExportRange()     { return _instance._sigma_export_range; }
  static double SigmaExportDelta()     { return _instance._sigma_export_delta; }
  static bool   SigmaExportConverged() { return _instance._sigma_export_converged; }
  // Sigma_c matrix export
  static bool SigmaMatrixExport() { return _instance._sigma_matrix_export; }
  // RPA Energies import/export
  static bool RPAEnergiesImport() { return _instance._rpa_energies_import; }
  static bool RPAEnergiesExport() { return _instance._rpa_energies_export; }
  // Misc.
  static bool ExportBinary() { return _instance._export_binary; }
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
