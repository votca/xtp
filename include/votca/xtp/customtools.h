#ifndef _VOTCA_XTP_CUSTOM_OPTS_H
#define	_VOTCA_XTP_CUSTOM_OPTS_H

#include <votca/tools/property.h>
#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {

class CustomOpts {
  private:
    static CustomOpts _instance;
    CustomOpts() {}
    void Parse(tools::Property& options);
    void Report();
    
    bool   _hedin      = false;
    bool   _gsc_export = false; // Export g self-consistency intermediate results
    double _gsc_alpha  = 0.0;
  
  public:
    static void Load();
    
    static bool Hedin() { return _instance._hedin; }
    static bool GSCExport() { return _instance._gsc_export; }
    static double GSCAlpha() { return _instance._gsc_alpha; }
};

class GSCLogger {
  private:
    static int _size;
    static void Log(const Eigen::VectorXd& row);
  
  public:
    static void Initialize(int size);
    static void LogFrequencies(const Eigen::VectorXd& frequencies);
    static void LogConverged(bool conv);
};

}
}

#endif	/* _VOTCA_XTP_CUSTOM_TOOLS_H */
