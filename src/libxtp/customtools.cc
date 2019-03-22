#include <votca/xtp/customtools.h>
#include <iostream>
#include <fstream>
#include <string.h>

namespace votca {
namespace xtp {

/* CustomOpts */

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
  _sigc_export = options.ifExistsReturnElseReturnDefault<Eigen::Vector3d>("customopts.sigc_export", Eigen::Vector3d::Zero());
}

void CustomOpts::Report() {
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", "");
  std::cout << std::endl << "Hedin:       " << _hedin;
  std::cout << std::endl << "GSC Export:  " << _gsc_export;
  std::cout << std::endl << "GSC Alpha:   " << _gsc_alpha;
  std::cout << std::endl << "SigC export: " << _sigc_export.format(fmt);
}

/* GSCLogger */

int GSCLogger::_size = 0;

void GSCLogger::Initialize(int size) {
  GSCLogger::_size = size; // Number of states (= number of columns)
  std::remove("gsc.log"); // Remove existing log file
}

void GSCLogger::LogFrequencies(const Eigen::VectorXd& frequencies) {
  GSCLogger::Log(frequencies);
}

void GSCLogger::LogConverged(bool conv) {
  if (conv) {
    GSCLogger::Log(Eigen::VectorXd::Ones(GSCLogger::_size));
  } else {
    GSCLogger::Log(Eigen::VectorXd::Zero(GSCLogger::_size));
  }
}

// Keep file in comma-separated matrix style so that it is easy to read
void GSCLogger::Log(const Eigen::VectorXd& row) {
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", "");
  std::ofstream gsc_log;
  gsc_log.open("gsc.log", std::ios_base::app);
  gsc_log << row.format(fmt) << std::endl;
  gsc_log.close();
}

}
}
