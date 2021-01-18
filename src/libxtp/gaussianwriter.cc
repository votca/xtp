#include <sstream>

#include "votca/xtp/basisset.h"
#include "votca/xtp/gaussianwriter.h"
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

namespace votca {
namespace xtp {

/*
 * This function converts VOTCA's L enum to the gaussian equivalent.
 * Gaussian uses a minus sign to indicate spherical shells.
 */
Index GaussianWriter::toGaussianL(L l) const {
  switch (l) {
    case L::S:
      return 0;
    case L::P:
      return 1;
    default:
      return -1*static_cast<Index>(l);
  }
}

std::string GaussianWriter::reorderedMOCoefficients(
    const Orbitals& orbitals) const {
  // Setup the reordering parameters
  std::array<Index, 49> multipliers{{
            1, //s
            1,-1,1, //p
            1,-1,1,-1,1, //d
            1,-1,1,-1,1,-1,1, //f 
            1,-1,1,-1,1,-1,1,-1,1, //g
            1,1,1,1,1,1,1,1,1,1,1, // h
            1,1,1,1,1,1,1,1,1,1,1,1,1 // i
  }};
  // clang-format off
  // the ordering of the m quantumnumbers for every shell in gaussian
  std::array<Index, 49> gaussianOrder={{
            0, //s
            1,0,-1, //p
            0,1,-1,2,-2, //d
            0,1,-1,2,-2,3,-3, //f 
            0,1,-1,2,-2,3,-3,4,-4, //g
            0,1,-1,2,-2,3,-3,4,-4,5,-5, // h
            0,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6 // i
  }};
  // clang-format on
  OrbReorder reorder(gaussianOrder, multipliers,true);
  Eigen::MatrixXd moCoefficients = orbitals.MOs().eigenvectors();
  reorder.reorderOrbitals(moCoefficients, orbitals.SetupDftBasis());

  // put the reordered mos in a string
  std::stringstream mos_string;

  int temp_int = 1;
  for (Index i = 0; i < moCoefficients.rows(); ++i) {
    for (Index j = 0; j < moCoefficients.cols(); ++j) {
      mos_string << boost::format("%16.8e") % moCoefficients(j,i);
      if (temp_int % 5 == 0) {
        mos_string << "\n";
      }
      temp_int++;
    }
  }
  mos_string << ((temp_int - 1) % 5 == 0 ? "" : "\n");

  return mos_string.str();
}

void GaussianWriter::WriteFile(const std::string& basename,
                               const Orbitals& orbitals) const {
  if (!orbitals.hasDFTbasisName()) {
    throw std::runtime_error(".orb file does not contain a basisset name");
  }

  AOBasis basis = orbitals.SetupDftBasis();

  std::ofstream outFile(basename + ".fchk");
  if (outFile.is_open()) {
    int temp_int;
    // job description
    outFile << basename << ", fchk created by VOTCA-XTP\n";
    outFile << "DUMMY_TYPE    DUMMY_METHOD    " << orbitals.getDFTbasisName()
            << "\n";

    // clang-format off
    outFile << boost::format("%-43s%-2s%15d\n") % "number of atoms" % "I" % orbitals.QMAtoms().size();
    outFile << boost::format("%-43s%-2s%15d\n") % "Charge" % "I" % 0;
    outFile << boost::format("%-43s%-2s%15d\n") % "Multiplicity" % "I" % 1;
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of electrons" % "I" % (2*orbitals.getNumberOfAlphaElectrons());
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of alpha electrons" % "I" % orbitals.getNumberOfAlphaElectrons();
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of beta electrons" % "I" % orbitals.getNumberOfAlphaElectrons();
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of basis functions" % "I" % basis.AOBasisSize();
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of independent functions" % "I" % basis.AOBasisSize();
    // clang-format on
    // ATOMIC NUMBERS
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Atomic numbers" % "I" %
                   orbitals.QMAtoms().size();
    temp_int = 1;
    for (const auto& atom : orbitals.QMAtoms()) {
      outFile << boost::format("%12d") % atom.getElementNumber();
      if (temp_int % 6 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 6 == 0 ? "" : "\n");
    // NUCLEAR CHARGES
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Nuclear charges" % "R" %
                   orbitals.QMAtoms().size();
    temp_int = 1;
    for (const auto& atom : orbitals.QMAtoms()) {
      outFile << boost::format("%16.8e") % (double)atom.getNuccharge();
      if (temp_int % 5 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // CURRENT CARTESIAN COORDINATES
    outFile << boost::format("%-43s%-2s N=  %10d\n") %
                   "Current cartesian coordinates" % "R" %
                   (3 * orbitals.QMAtoms().size());
    temp_int = 1;
    for (const auto& atom : orbitals.QMAtoms()) {
      for (int i = 0; i < 3; ++i) {
        outFile << boost::format("%16.8e") % atom.getPos()(i);
        if (temp_int % 5 == 0) {
          outFile << "\n";
        }
        temp_int++;
      }
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // NUMBER OF PRIMITIVE SHELLS
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of primitive shells" %
                   "I" % basis.getNumberOfPrimitives();
    // NUMBER OF CONTRACTED SHELLS
    outFile << boost::format("%-43s%-2s%15d\n") %
                   "Number of contracted shells" % "I" % basis.getNumofShells();
    // PURE/CARTESIAN D
    outFile << boost::format("%-43s%-2s%15d\n") % "Pure/Cartesian d shells " %
                   "I" % 0;
    // PURE/CARTESIAN F
    outFile << boost::format("%-43s%-2s%15d\n") % "Pure/Cartesian f shells " %
                   "I" % 0;
    // HIGHEST ANGULAR MOMENTUM
    outFile << boost::format("%-43s%-2s%15d\n") % "Highest angular momentum " %
                   "I" % basis.getMaxL();
    // HIGHEST ANGULAR MOMENTUM
    outFile << boost::format("%-43s%-2s%15d\n") %
                   "Largest degree of contraction " % "I" % basis.getMaxNprim();
    // SHELL TYPES
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Shell types" % "I" %
                   (basis.getNumofShells());
    temp_int = 1;
    for (const auto& shell : basis) {
      outFile << boost::format("%12d") % toGaussianL(shell.getL());
      if (temp_int % 6 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 6 == 0 ? "" : "\n");
    // NR OF PRIMITIVES PER SHELL
    outFile << boost::format("%-43s%-2s N=  %10d\n") %
                   "Number of primitives per shell" % "I" %
                   (basis.getNumofShells());
    temp_int = 1;
    for (const AOShell& shell : basis) {
      outFile << boost::format("%12d") % shell.getSize();
      if (temp_int % 6 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 6 == 0 ? "" : "\n");
    // SHELL TO ATOM MAP
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Shell to atom map" %
                   "I" % (basis.getNumofShells());
    temp_int = 1;
    for (const AOShell& shell : basis) {
      // Gaussian indices start at 1, hence the + 1
      outFile << boost::format("%12d") % (shell.getAtomIndex() + 1);
      if (temp_int % 6 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 6 == 0 ? "" : "\n");
    // PRIMITIVE EXPONENTS
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Primitive exponents" %
                   "R" % basis.getNumberOfPrimitives();
    temp_int = 1;
    for (const AOShell& shell : basis) {
      for (const AOGaussianPrimitive& prim : shell) {
        outFile << boost::format("%16.8e") % prim.getDecay();
        if (temp_int % 5 == 0) {
          outFile << "\n";
        }
        temp_int++;
      }
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // CONTRACTION COEFFICIENTS
    outFile << boost::format("%-43s%-2s N=  %10d\n") %
                   "Contraction coefficients" % "R" %
                   basis.getNumberOfPrimitives();
    temp_int = 1;
    for (const AOShell& shell : basis) {
      for (const AOGaussianPrimitive& prim : shell) {
        outFile << boost::format("%16.8e") % prim.getContraction();
        if (temp_int % 5 == 0) {
          outFile << "\n";
        }
        temp_int++;
      }
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // SHELL COORDINATES
    outFile << boost::format("%-43s%-2s N=  %10d\n") %
                   "Coordinates of each shell" % "R" %
                   (3 * basis.getNumofShells());
    temp_int = 1;
    for (const AOShell& shell : basis) {
      for (int i = 0; i < 3; ++i) {
        outFile << boost::format("%16.8e") % shell.getPos()(i);
        if (temp_int % 5 == 0) {
          outFile << "\n";
        }
        temp_int++;
      }
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // TOTAL ENERGY
    outFile << boost::format("%-43s%-2s%22.15e\n") % "Total Energy" % "R" %
                   orbitals.getDFTTotalEnergy();
    // ALPHA ORBITAL ENERGIES
    outFile << boost::format("%-43s%-2s N=  %10d\n") %
                   "Alpha Orbital Energies" % "R" %
                   orbitals.MOs().eigenvalues().size();
    temp_int = 1;
    for (Index i = 0; i < orbitals.MOs().eigenvalues().size(); ++i) {
      outFile << boost::format("%16.8e") % orbitals.MOs().eigenvalues()[i];
      if (temp_int % 5 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // ALPHA ORBITAL ENERGIES
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Alpha MO coefficients" %
                   "R" %
                   (orbitals.MOs().eigenvalues().size() *
                    orbitals.MOs().eigenvalues().size());
    outFile << reorderedMOCoefficients(orbitals);
    // SCF DENSITY
  }
}

}  // namespace xtp
}  // namespace votca