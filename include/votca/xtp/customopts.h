#ifndef _VOTCA_XTP_CUSTOM_OPTS_H
#define	_VOTCA_XTP_CUSTOM_OPTS_H

#include <votca/tools/property.h>

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

}
}

#endif	/* _VOTCA_XTP_CUSTOM_OPTS_H */
