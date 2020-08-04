#include "orb2mol.h"

#include <boost/format.hpp>

namespace votca {
namespace xtp {

void Orb2Mol::Initialize(const tools::Property& user_options) {
  tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);

  _job_name = options.ifExistsReturnElseReturnDefault<std::string>("job_name",
                                                                   _job_name);
  _moldenfile = _job_name + ".molden.input";
  _orbfile = _job_name + ".orb";
  _xyzfile = _job_name + ".xyz";
}

void Orb2Mol::writeAtoms(Orbitals& orbitals, std::ofstream& outFile) {
  for (auto& atom : orbitals.QMAtoms()) {
    Eigen::Vector3d pos = atom.getPos();
    outFile << boost::format("%-4s %5d %5d %22.12e %22.10e %22.10e\n") %
                   atom.getElement() % (atom.getId() + 1) %
                   atom.getPureNucCharge() % pos[0] % pos[1] % pos[2];
  }
}

void Orb2Mol::writeMOs(Orbitals& orbitals, std::ofstream& outFile) { ; }

void Orb2Mol::writeBasisSet(Orbitals& orbitals, std::ofstream& outFile) {
  if (orbitals.hasDFTbasisName()) {
    BasisSet bs;
    bs.Load(orbitals.getDFTbasisName());

    for (auto& atom : orbitals.QMAtoms()) {
      const Element& element = bs.getElement(atom.getElement());
      // The 0 in the format string of the next line is meaningless it
      // is included for backwards compatibility of molden files
      outFile <<  boost::format("%4d 0 \n") % (atom.getId() + 1);
      for (const Shell& shell : element) {
        for (const char& subtype : shell.getType()) {
          // The 1.0 at the end of the next line is meaningless it
          // is included for backwards compatibility of molden files
          outFile << boost::format("%-3s %4d %3.1f \n") % std::tolower(subtype, std::locale())  % shell.getSize() % 1.0; 
          for (const GaussianPrimitive& gaussian : shell) {
            outFile << boost::format("%22.10e %22.10e\n") %
                           gaussian.decay() %
                           gaussian.Contractions()[FindLmax(
                               std::string(1, subtype))];
          }
        }
      }
      outFile << " \n";
    }

  } else {
    throw std::runtime_error(".orb file does not contain a basisset name");
  }
}

bool Orb2Mol::Evaluate() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  Orbitals orbitals;
  XTP_LOG(Log::error, _log) << "Loading data from " << _orbfile << std::flush;
  orbitals.ReadFromCpt(_orbfile);

  std::ofstream outFile(_moldenfile);

  if (outFile.is_open()) {
    // print Header
    outFile << "[Molden Format]\n";
    outFile << "[Title]\n";
    outFile << "Molden file created by VOTCA-XTP for basename: " << _job_name
            << "\n";
    outFile << " \n";

    outFile << "[Atoms] AU\n";
    writeAtoms(orbitals, outFile);

    outFile << "[GTO] \n";
    writeBasisSet(orbitals, outFile);

    // indicate spherical D F and G functions
    outFile << "[5D] \n[7F] \n[9G] \n";

    outFile << "[MO]\n";
    writeMOs(orbitals, outFile);
    return true;
  }
  return false;
}

}  // namespace xtp
}  // namespace votca