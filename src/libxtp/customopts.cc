#include <votca/xtp/customopts.h>
#include <iostream>

namespace votca {
namespace xtp {

CustomOpts CustomOpts::_instance;

void CustomOpts::Load() {
  tools::Property options;
  try {
    load_property_from_xml(options, "customopts.xml");
  } catch(...) {
    return;
  }
  CustomOpts::_instance.Parse(options);
  std::cout << std::endl << "************************************************";
  std::cout << std::endl << "Loaded custom options";
  CustomOpts::_instance.Report();
  std::cout << std::endl << "************************************************";
}

void CustomOpts::Parse(tools::Property& options) {
  _hedin = options.ifExistsReturnElseReturnDefault<bool>("customopts.hedin", _hedin);
  _gsc_export = options.ifExistsReturnElseReturnDefault<bool>("customopts.gsc_export", _gsc_export);
  _gsc_alpha = options.ifExistsReturnElseReturnDefault<double>("customopts.gsc_alpha", _gsc_alpha);
}

void CustomOpts::Report() {
  std::cout << std::endl << "Hedin:      " << _hedin;
  std::cout << std::endl << "GSC Export: " << _gsc_export;
  std::cout << std::endl << "GSC Alpha:  " << _gsc_alpha;
}

}
}
