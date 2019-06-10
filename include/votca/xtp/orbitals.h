/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __VOTCA_XTP_ORBITALS_H
#define __VOTCA_XTP_ORBITALS_H

#include <votca/xtp/eigen.h>

#include <votca/xtp/qmatom.h>

#include <votca/tools/globals.h>
#include <votca/tools/property.h>
#include <votca/tools/vec.h>
#include <votca/xtp/checkpoint.h>

#include <boost/format.hpp>
#include <votca/ctp/logger.h>
#include <votca/tools/constants.h>
#include <votca/xtp/qmstate.h>

#include "basisset.h"

namespace votca {
namespace xtp {

/**
 * \brief container for molecular orbitals
 *
 * The Orbitals class stores orbital id, energy, MO coefficients, basis set
 *
 */
class Orbitals {
 public:
  Orbitals();

  Orbitals(const Orbitals& orbital);

  Orbitals& operator=(const Orbitals& orbital);

  ~Orbitals();

  // functions for analyzing fragment charges via Mulliken populations
  static Eigen::VectorXd LoewdinPopulation(const Eigen::MatrixXd& densitymatrix,
                                           const Eigen::MatrixXd& overlapmatrix,
                                           int frag);

  bool hasBasisSetSize() const { return (_basis_set_size > 0) ? true : false; }

  int getBasisSetSize() const { return _basis_set_size; }

  void setBasisSetSize(int basis_set_size) { _basis_set_size = basis_set_size; }

  int getLumo() const { return _occupied_levels; }

  int getHomo() const { return _occupied_levels - 1; }
  // access to DFT number of levels, new, tested

  bool hasNumberOfLevels() const {
    return ((_occupied_levels > 0) ? true : false);
  }

  void setNumberOfOccupiedLevels(int occupied_levels);

  // access to DFT number of electrons, new, tested

  bool hasNumberOfAlphaElectrons() const {
    return (_number_alpha_electrons > 0) ? true : false;
  }

  int getNumberOfAlphaElectrons() const { return _number_alpha_electrons; };

  void setNumberOfAlphaElectrons(int electrons) {
    _number_alpha_electrons = electrons;
  }

  bool hasECPName() const { return (_ECP != "") ? true : false; }

  const std::string& getECPName() const { return _ECP; };

  void setECPName(const std::string& ECP) { _ECP = ECP; };

  // access to QM package name, new, tested

  bool hasQMpackage() const { return (!_qm_package.empty()); }

  const std::string& getQMpackage() const { return _qm_package; }

  void setQMpackage(const std::string& qmpackage) { _qm_package = qmpackage; }

  // access to DFT AO overlap matrix, new, tested

  bool hasAOOverlap() const { return (_overlap.rows() > 0) ? true : false; }

  const Eigen::MatrixXd& AOOverlap() const { return _overlap; }

  Eigen::MatrixXd& AOOverlap() { return _overlap; }

  // access to DFT molecular orbital energies, new, tested

  bool hasMOEnergies() const {
    return (_mo_energies.size() > 0) ? true : false;
  }

  const Eigen::VectorXd& MOEnergies() const { return _mo_energies; }

  Eigen::VectorXd& MOEnergies() { return _mo_energies; }

  // access to DFT molecular orbital energy of a specific level (in eV)
  double getEnergy(int level) const {
    return (hasMOEnergies()) ? votca::tools::conv::hrt2ev * _mo_energies[level]
                             : 0;
  }

  // access to DFT molecular orbital coefficients, new, tested
  bool hasMOCoefficients() const {
    return (_mo_coefficients.cols() > 0) ? true : false;
  }

  const Eigen::MatrixXd& MOCoefficients() const { return _mo_coefficients; }

  Eigen::MatrixXd& MOCoefficients() { return _mo_coefficients; }

  // determine (pseudo-)degeneracy of a DFT molecular orbital
  std::vector<int> CheckDegeneracy(int level, double energy_difference) const;

  // access to QM atoms
  bool hasQMAtoms() { return (_atoms.size() > 0) ? true : false; }

  const std::vector<QMAtom*>& QMAtoms() const { return _atoms; }

  std::vector<QMAtom*>& QMAtoms() { return _atoms; }

  // access to classical self-energy in MM environment, new, tested
  bool hasSelfEnergy() const { return (_self_energy != 0.0) ? true : false; }

  double getSelfEnergy() const { return _self_energy; }

  void setSelfEnergy(double selfenergy) { _self_energy = selfenergy; }

  // access to QM total energy, new, tested
  bool hasQMEnergy() const { return (_qm_energy != 0.0) ? true : false; }

  double getQMEnergy() const { return _qm_energy; }

  void setQMEnergy(double qmenergy) { _qm_energy = qmenergy; }

  // access to DFT basis set name

  bool hasDFTbasisName() const { return (!_dftbasis.empty()) ? true : false; }

  void setDFTbasisName(const std::string basis) { _dftbasis = basis; }

  const std::string& getDFTbasisName() const { return _dftbasis; }

  /*
   *  ======= GW-BSE related functions =======
   */

  // access to exchange-correlation AO matrix, new, tested

  bool hasAOVxc() const { return (_vxc.rows() > 0) ? true : false; }

  Eigen::MatrixXd& AOVxc() { return _vxc; }

  const Eigen::MatrixXd& AOVxc() const { return _vxc; }

  // access to auxiliary basis set name

  bool hasAuxbasisName() const { return (!_auxbasis.empty()) ? true : false; }

  void setAuxbasisName(std::string basis) { _auxbasis = basis; }

  const std::string& getAuxbasisName() const { return _auxbasis; }

  // access to list of indices used in GWA

  bool hasGWAindices() const { return (_qpmax > 0) ? true : false; }

  void setGWindices(int qpmin, int qpmax) {
    _qpmin = qpmin;
    _qpmax = qpmax;
  }

  int getGWAmin() const { return _qpmin; }

  int getGWAmax() const { return _qpmax; }

  // access to list of indices used in RPA

  bool hasRPAindices() const { return (_rpamax > 0) ? true : false; }

  void setRPAindices(int rpamin, int rpamax) {
    _rpamin = rpamin;
    _rpamax = rpamax;
  }

  int getRPAmin() const { return _rpamin; }

  int getRPAmax() const { return _rpamax; }

  // access to list of indices used in BSE

  void setTDAApprox(bool usedTDA) { _useTDA = usedTDA; }
  bool getTDAApprox() const { return _useTDA; }

  bool hasBSEindices() const { return (_bse_cmax > 0) ? true : false; }

  void setBSEindices(int vmin, int cmax) {
    _bse_vmin = vmin;
    _bse_vmax = getHomo();
    _bse_cmin = getLumo();
    _bse_cmax = cmax;
    _bse_vtotal = _bse_vmax - _bse_vmin + 1;
    _bse_ctotal = _bse_cmax - _bse_cmin + 1;
    _bse_size = _bse_vtotal * _bse_ctotal;
    return;
  }

  int getBSEvmin() const { return _bse_vmin; }

  int getBSEvmax() const { return _bse_vmax; }

  int getBSEcmin() const { return _bse_cmin; }

  int getBSEcmax() const { return _bse_cmax; }

  double getScaHFX() const { return _ScaHFX; }

  void setScaHFX(double ScaHFX) { _ScaHFX = ScaHFX; }

  // access to perturbative QP energies

  bool hasQPpert() const {
    return (_QPpert_energies.size() > 0) ? true : false;
  }

  const Eigen::MatrixXd& QPpertEnergies() const { return _QPpert_energies; }

  Eigen::MatrixXd& QPpertEnergies() { return _QPpert_energies; }

  // access to diagonalized QP energies and wavefunctions

  bool hasQPdiag() const {
    return (_QPdiag_energies.size() > 0) ? true : false;
  }

  const Eigen::VectorXd& QPdiagEnergies() const { return _QPdiag_energies; }

  Eigen::VectorXd& QPdiagEnergies() { return _QPdiag_energies; }

  const Eigen::MatrixXd& QPdiagCoefficients() const {
    return _QPdiag_coefficients;
  }

  Eigen::MatrixXd& QPdiagCoefficients() { return _QPdiag_coefficients; }

  bool hasBSETriplets() const {
    return (_BSE_triplet_energies.cols() > 0) ? true : false;
  }

  const Eigen::VectorXd& BSETripletEnergies() const {
    return _BSE_triplet_energies;
  }

  Eigen::VectorXd& BSETripletEnergies() { return _BSE_triplet_energies; }

  const Eigen::MatrixXd& BSETripletCoefficients() const {
    return _BSE_triplet_coefficients;
  }

  Eigen::MatrixXd& BSETripletCoefficients() {
    return _BSE_triplet_coefficients;
  }

  const Eigen::MatrixXd& BSETripletCoefficientsAR() const {
    return _BSE_triplet_coefficients_AR;
  }

  Eigen::MatrixXd& BSETripletCoefficientsAR() {
    return _BSE_triplet_coefficients_AR;
  }

  // access to singlet energies and wave function coefficients

  bool hasBSESinglets() const {
    return (_BSE_singlet_energies.cols() > 0) ? true : false;
  }

  const Eigen::VectorXd& BSESingletEnergies() const {
    return _BSE_singlet_energies;
  }

  Eigen::VectorXd& BSESingletEnergies() { return _BSE_singlet_energies; }

  const Eigen::MatrixXd& BSESingletCoefficients() const {
    return _BSE_singlet_coefficients;
  }

  Eigen::MatrixXd& BSESingletCoefficients() {
    return _BSE_singlet_coefficients;
  }

  // for anti-resonant part in full BSE

  const Eigen::MatrixXd& BSESingletCoefficientsAR() const {
    return _BSE_singlet_coefficients_AR;
  }

  Eigen::MatrixXd& BSESingletCoefficientsAR() {
    return _BSE_singlet_coefficients_AR;
  }

  // access to transition dipole moments

  bool hasTransitionDipoles() const {
    return (_transition_dipoles.size() > 0) ? true : false;
  }

  const std::vector<tools::vec>& TransitionDipoles() const {
    return _transition_dipoles;
  }

  std::vector<tools::vec>& TransitionDipoles() { return _transition_dipoles; }

  std::vector<double> Oscillatorstrengths() const;

  Eigen::Vector3d CalcCoM() const;

  Eigen::Vector3d CalcElDipole(const QMState& state);

  // Calculates full electron density for state or transition density, if you
  // want to calculate only the density contribution of hole or electron use
  // DensityMatrixExcitedState
  Eigen::MatrixXd DensityMatrixFull(const QMState& state) const;

  // functions for calculating density matrices
  Eigen::MatrixXd DensityMatrixGroundState() const;
  std::vector<Eigen::MatrixXd> DensityMatrixExcitedState(
      const QMState& state) const;
  Eigen::MatrixXd DensityMatrixQuasiParticle(const QMState& state) const;
  Eigen::MatrixXd CalculateQParticleAORepresentation() const;
  double getTotalStateEnergy(const QMState& state) const;    // Hartree
  double getExcitedStateEnergy(const QMState& state) const;  // Hartree

  // access to fragment charges of singlet excitations
  bool hasFragmentChargesSingEXC() const {
    return (_DqS_frag.size() > 0) ? true : false;
  }

  const std::vector<Eigen::VectorXd>& getFragmentChargesSingEXC() const {
    return _DqS_frag;
  }

  void setFragmentChargesSingEXC(const std::vector<Eigen::VectorXd>& DqS_frag) {
    _DqS_frag = DqS_frag;
  }

  // access to fragment charges of triplet excitations
  bool hasFragmentChargesTripEXC() const {
    return (_DqT_frag.size() > 0) ? true : false;
  }

  const std::vector<Eigen::VectorXd>& getFragmentChargesTripEXC() const {
    return _DqT_frag;
  }

  void setFragmentChargesTripEXC(const std::vector<Eigen::VectorXd>& DqT_frag) {
    _DqT_frag = DqT_frag;
  }

  // access to fragment charges in ground state

  const Eigen::VectorXd& getFragmentChargesGS() const { return _GSq_frag; }

  void setFragmentChargesGS(const Eigen::VectorXd& GSq_frag) {
    _GSq_frag = GSq_frag;
  }

  void setFragment_E_localisation_singlet(
      const std::vector<Eigen::VectorXd>& popE) {
    _popE_s = popE;
  }

  void setFragment_H_localisation_singlet(
      const std::vector<Eigen::VectorXd>& popH) {
    _popH_s = popH;
  }

  void setFragment_E_localisation_triplet(
      const std::vector<Eigen::VectorXd>& popE) {
    _popE_t = popE;
  }

  void setFragment_H_localisation_triplet(
      const std::vector<Eigen::VectorXd>& popH) {
    _popE_s = popH;
  }

  const std::vector<Eigen::VectorXd>& getFragment_E_localisation_singlet()
      const {
    return _popE_s;
  }
  const std::vector<Eigen::VectorXd>& getFragment_H_localisation_singlet()
      const {
    return _popH_s;
  }
  const std::vector<Eigen::VectorXd>& getFragment_E_localisation_triplet()
      const {
    return _popE_t;
  }
  const std::vector<Eigen::VectorXd>& getFragment_H_localisation_triplet()
      const {
    return _popH_t;
  }

  void PrepareDimerGuess(const Orbitals& orbitalsA, const Orbitals& orbitalsB);

  Eigen::VectorXd FragmentNuclearCharges(int frag) const;

  // returns indeces of a re-sorted vector of energies from lowest to highest
  std::vector<int> SortEnergies();

  QMAtom* AddAtom(int AtomID, std::string type, tools::vec pos) {
    QMAtom* pAtom = new QMAtom(AtomID, type, pos);
    _atoms.push_back(pAtom);
    return pAtom;
  }

  QMAtom* AddAtom(QMAtom atom) {
    QMAtom* pAtom = new QMAtom(atom);
    _atoms.push_back(pAtom);
    return pAtom;
  }

  void OrderMOsbyEnergy();

  void WriteXYZ(const std::string& filename,
                std::string header = "GENERATED BY VOTCA::XTP") const;

  void LoadFromXYZ(const std::string& filename);

  void WriteToCpt(const std::string& filename) const;

  void ReadFromCpt(const std::string& filename);

 private:
  void copy(const Orbitals& orbital);

  void WriteToCpt(CheckpointFile f) const;
  void WriteToCpt(CheckpointWriter w) const;

  void ReadFromCpt(CheckpointFile f);
  void ReadFromCpt(CheckpointReader parent);

  Eigen::MatrixXd TransitionDensityMatrix(const QMState& state) const;
  std::vector<Eigen::MatrixXd> DensityMatrixExcitedState_R(
      const QMState& state) const;
  std::vector<Eigen::MatrixXd> DensityMatrixExcitedState_AR(
      const QMState& state) const;
  Eigen::MatrixXd CalcAuxMat_cc(const Eigen::VectorXd& coeffs) const;
  Eigen::MatrixXd CalcAuxMat_vv(const Eigen::VectorXd& coeffs) const;

  int _basis_set_size;
  int _occupied_levels;
  int _number_alpha_electrons;
  std::string _ECP;
  bool _useTDA;

  Eigen::VectorXd _mo_energies;
  Eigen::MatrixXd _mo_coefficients;

  Eigen::MatrixXd _overlap;
  Eigen::MatrixXd _vxc;

  std::vector<QMAtom*> _atoms;

  double _qm_energy;
  double _self_energy;

  // new variables for GW-BSE storage
  int _rpamin;
  int _rpamax;

  int _qpmin;
  int _qpmax;

  int _bse_vmin;
  int _bse_vmax;
  int _bse_cmin;
  int _bse_cmax;
  int _bse_size;
  int _bse_vtotal;
  int _bse_ctotal;

  double _ScaHFX;

  std::string _dftbasis;
  std::string _auxbasis;

  std::string _qm_package;

  // perturbative quasiparticle energies
  Eigen::MatrixXd _QPpert_energies;

  // quasiparticle energies and coefficients after diagonalization
  Eigen::VectorXd _QPdiag_energies;
  Eigen::MatrixXd _QPdiag_coefficients;
  // excitons
  Eigen::VectorXd _BSE_singlet_energies;
  Eigen::MatrixXd _BSE_singlet_coefficients;
  Eigen::MatrixXd _BSE_singlet_coefficients_AR;

  std::vector<tools::vec> _transition_dipoles;
  Eigen::VectorXd _BSE_triplet_energies;
  Eigen::MatrixXd _BSE_triplet_coefficients;
  Eigen::MatrixXd _BSE_triplet_coefficients_AR;

  std::vector<Eigen::VectorXd> _DqS_frag;  // fragment charge changes in exciton

  std::vector<Eigen::VectorXd> _DqT_frag;

  Eigen::VectorXd _GSq_frag;  // ground state effective fragment charges

  std::vector<Eigen::VectorXd> _popE_s;
  std::vector<Eigen::VectorXd> _popE_t;
  std::vector<Eigen::VectorXd> _popH_s;
  std::vector<Eigen::VectorXd> _popH_t;
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_ORBITALS_H */
