#ifndef _VOTCA_XTP_GSC_LOGGER_H
#define	_VOTCA_XTP_GSC_LOGGER_H

#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {

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

#endif	/* _VOTCA_XTP_GSC_LOGGER_H */
