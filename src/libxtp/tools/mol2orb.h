

#pragma once
#ifndef VOTCA_XTP_MOL2ORB_PRIVATE_H
#define VOTCA_XTP_MOL2ORB_PRIVATE_H

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/qmtool.h"
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {
class Orbitals;

class Mol2Orb : public QMTool {
 public:
  Mol2Orb() = default;

  ~Mol2Orb() override = default;

  std::string Identify() override { return "mol2orb"; }

  void Initialize(const tools::Property& user_options) override;
  bool Evaluate() override;

 private:
   // clang-format off
  std::array<Index,25> _multipliers={
            1, //s
            1,1,1, //p
            1,1,1,1,1, //d
            1,1,1,1,1,-1,-1, //f 
            1,1,1,1,1,-1,-1,-1,-1 //g
            };
  std::vector<std::array<int, 2>> _twoCyclesS{};
  std::vector<std::array<int, 2>> _twoCyclesP{
    std::array<int, 2>{0, 2}
  };
  std::vector<std::array<int, 2>> _twoCyclesD{
    std::array<int, 2>{1, 2},
    std::array<int, 2>{3, 4}
    };
  std::vector<std::array<int, 2>> _twoCyclesF{
    std::array<int, 2>{1, 2},
    std::array<int, 2>{3, 4},
    std::array<int, 2>{5, 6}
  };
  std::vector<std::array<int, 2>> _twoCyclesG{
    std::array<int, 2>{1, 2},
    std::array<int, 2>{3, 4},
    std::array<int, 2>{5, 6},
    std::array<int, 2>{7,8}
  };
  // clang-format on
  std::string _moldenfile;
  std::string _orbfile;
  std::string _xyzfile;
  std::string _basisset_name;
  std::string _aux_basisset_name;
  AOBasis _basis;
  AOBasis _auxbasis;
  Logger _log;

  inline std::string readAtoms(QMMolecule& mol, std::string units,
                               std::ifstream& input_file) const;
  inline std::string readMOs(Orbitals& orbitals, std::ifstream& input_file);
  inline std::vector<std::array<int, 2>> getTwoCycles(Index numFunc);
  void setupBasis(Orbitals& orbitals);
  void reorderOrbitals(Eigen::MatrixXd& v);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MOL2ORB_PRIVATE_H