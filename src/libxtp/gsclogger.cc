#include <votca/xtp/gsclogger.h>
#include <fstream>

namespace votca {
namespace xtp {

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
  std::ofstream g_sc_log;
  g_sc_log.open("gsc.log", std::ios_base::app);
  g_sc_log << row.format(fmt) << std::endl;
  g_sc_log.close();
}

}
}
