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

#include "votca/xtp/ERIs.h"
#include "votca/xtp/adiis.h"
#include "votca/xtp/diis.h"
#include <fstream>
#include <votca/tools/property.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aomatrix3d.h>
#include <votca/xtp/aopotential.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/multishift.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/padeapprox.h>
#include <votca/xtp/qmmolecule.h>
#include <votca/xtp/sternheimer.h>
#include <votca/xtp/vxc_grid.h>
#include <votca/xtp/vxc_potential.h>
namespace votca {
namespace xtp {

void Sternheimer::setUpMatrices() {

  // saving matrices needed from orbitals
  this->_num_occ_lvls = _orbitals.getNumberOfAlphaElectrons();
  this->_basis_size = _orbitals.getBasisSetSize();
  this->_overlap_Matrix = OverlapMatrix();
  this->_density_Matrix = _orbitals.DensityMatrixGroundState();
  this->_mo_coefficients = _orbitals.MOs().eigenvectors();
  this->_mo_energies = _orbitals.MOs().eigenvalues();
  this->_inverse_overlap = _overlap_Matrix.inverse();
  this->_Hamiltonian_Matrix = Hamiltonian();
}

void Sternheimer::configurate(const options_sternheimer& opt) { _opt = opt; }

void Sternheimer::initializeMultishift(Index size) {
  _multishift.setMatrixSize(size);
}

Eigen::MatrixXcd Sternheimer::OverlapMatrix() {

  AOBasis basis = _orbitals.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  return overlap.Matrix().cast<std::complex<double>>();
}

Eigen::MatrixXcd Sternheimer::DensityMatrix() {
  return _orbitals.DensityMatrixGroundState().cast<std::complex<double>>();
}

Eigen::MatrixXcd Sternheimer::Hamiltonian() {

  // rebuilding Hamitonian from eigenvectors and values
  const Eigen::MatrixXd& mo_coefficients = _orbitals.MOs().eigenvectors();
  return (_overlap_Matrix * mo_coefficients *
          _orbitals.MOs().eigenvalues().asDiagonal() *
          mo_coefficients.transpose() * _overlap_Matrix)
      .cast<std::complex<double>>();
}

std::vector<std::complex<double>> Sternheimer::BuildGrid(
    double omega_start, double omega_end, Index steps,
    double imaginary_shift) const {

  std::vector<std::complex<double>> grid;

  double stepsize = (omega_end - omega_start) / steps;
  const double ev2hrt = votca::tools::conv::ev2hrt;
  std::complex<double> d(1, 0);
  std::complex<double> i(0, 1);

  for (Index n = 0; n <= steps; n++) {
    // Converts from input eV to Hartree
    grid.push_back(omega_start * ev2hrt + n * stepsize * ev2hrt +
                   imaginary_shift * i * ev2hrt);
  }
  // Grid has form: stepsize*n+i*imaginary shift
  return grid;
}

Eigen::MatrixXcd Sternheimer::SternheimerLHS(
    const Eigen::MatrixXcd& hamiltonian,
    const Eigen::MatrixXcd& inverse_overlap, double eps,
    std::complex<double> omega, bool pm) const {
  Eigen::MatrixXcd Identity_cmplx =
      Eigen::MatrixXd::Identity(_basis_size, _basis_size)
          .cast<std::complex<double>>();
  std::complex<double> temp = eps + omega;
  if (pm != true) {
    temp = (eps - omega);
  }
  Eigen::MatrixXcd LHS =
      (inverse_overlap * hamiltonian - temp * Identity_cmplx);
  return LHS;
}

Eigen::VectorXcd Sternheimer::SternheimerRHS(
    const Eigen::MatrixXcd& inverse_overlap, const Eigen::MatrixXcd& density,
    const Eigen::MatrixXcd& pertubation, const Eigen::VectorXd& coeff) const {

  Eigen::VectorXcd RHS = -1 * (inverse_overlap - density) * pertubation * coeff;
  return RHS;
}

Eigen::MatrixXcd Sternheimer::DeltaNSC(
    std::complex<double> w, const Eigen::MatrixXcd& perturbation) const {

  // Setting up vectors to store old results for Anderson mixing and initial
  // perturbation
  std::vector<Eigen::MatrixXcd> perturbationVectorInput;
  std::vector<Eigen::MatrixXcd> perturbationVectoroutput;

  perturbationVectorInput.push_back(perturbation);
  Eigen::MatrixXcd perturbationUsed = perturbationVectorInput.back();

  Eigen::MatrixXcd delta_n_out_new =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  Eigen::MatrixXcd delta_n_out_old =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  Eigen::MatrixXcd delta_n_step_one =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);

  // Setting up ERIS for four center integral
  AOBasis dftbasis = _orbitals.SetupDftBasis();
  AOBasis auxbasis = _orbitals.SetupAuxBasis();
  ERIs eris;
  eris.Initialize(dftbasis, auxbasis);

  // Setting up Grid for Fxc functional
  Vxc_Grid grid;
  grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                 dftbasis);
  Vxc_Potential<Vxc_Grid> Vxcpot(grid);
  Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());

  // double alpha = 4*(_mo_energies(_mo_energies.size()-1)-_mo_energies(0));
  double alpha = 1000;
  // Loop until convergence
  for (Index n = 0; n < _opt.max_iterations_sc_sternheimer; n++) {

    // Matrices to store the solutions of the sternheimer equation
    Eigen::MatrixXcd solution_p =
        Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);
    Eigen::MatrixXcd solution_m =
        Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);
    // Loop over all occupied states
    for (Index v = 0; v < _num_occ_lvls; v++) {

      // Building RHS
      Eigen::MatrixXcd RHS =
          SternheimerRHS(_inverse_overlap, _density_Matrix, perturbationUsed,
                         _mo_coefficients.col(v));

      // std::cout<<"RHS "<< RHS <<std::endl;
      // Building LHS with +/- omega and solving the system
      Eigen::MatrixXcd LHS_P = SternheimerLHS(
          _Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, true);
      Eigen::MatrixXcd LHS_M = SternheimerLHS(
          _Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, false);
      if (true) {
        LHS_P = LHS_P + alpha * _density_Matrix.transpose() * _overlap_Matrix;
        LHS_M = LHS_M + alpha * _density_Matrix.transpose() * _overlap_Matrix;
      }
      solution_p.col(v) = LHS_P.colPivHouseholderQr().solve(RHS);
      solution_m.col(v) = LHS_M.colPivHouseholderQr().solve(RHS);
    }
    // Saving previous delta n
    delta_n_out_old = delta_n_out_new;
    // Calculating new delta n
    delta_n_out_new =
        2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) *
            solution_p.transpose() +
        2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) *
            solution_m.transpose();
    // Perfomring the to four center Integrals to update delta V
    Eigen::MatrixXcd contract =
        eris.ContractRightIndecesWithMatrix(delta_n_out_new);

    Eigen::MatrixXcd FxcInt =
        Vxcpot.IntegrateFXC(_density_Matrix, delta_n_out_new);

    // Check if max mixing history is reached and adding new step to history
    if (perturbationVectoroutput.size() > _opt.max_mixing_history - 1) {
      perturbationVectoroutput.erase(perturbationVectoroutput.begin());
    }

    perturbationVectoroutput.push_back((perturbation) + contract + FxcInt);

    double diff =
        (perturbationVectorInput.back() - perturbationVectoroutput.back())
            .squaredNorm();
    //std::cout << n << " " << diff << std::endl;
    if (diff < _opt.tolerance_sc_sternheimer) {
    //  std::cout << "Frequency: " << w << "Converged after " << n + 1
      //          << " iteration." << std::endl;

      // // throw std::exception();
      // Index occ = _orbitals.getNumberOfAlphaElectrons();
      // Eigen::MatrixXcd HmS = _Hamiltonian_Matrix - _overlap_Matrix *
      // _orbitals.MOs().eigenvalues().head(occ).asDiagonal(); Eigen::MatrixXcd
      // moc = _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls);
      // Eigen::MatrixXcd dmoc = solution_p;
      // std::complex<double> pulay1 = -2.0 * (dmoc.cwiseProduct(HmS *
      // moc)).sum(); std::complex<double> pulay2 = -2.0 *
      // (moc.cwiseProduct(HmS*dmoc)).sum(); std::cout << " \n Pulay " << pulay1
      // + pulay2 << std::endl;
      return delta_n_out_new;
    }
    // Mixing if at least in iteration 2
    if (n == 0) {
      perturbationUsed =
          _opt.mixing_constant * perturbationVectoroutput.back() +
          (1 - _opt.mixing_constant) * perturbationVectorInput.back();
      perturbationVectorInput.push_back(perturbationUsed);
    } else {

      perturbationUsed = (NPAndersonMixing(perturbationVectorInput,
                                           perturbationVectoroutput, 0.5));
      if (perturbationVectorInput.size() > _opt.max_mixing_history - 1) {
        perturbationVectorInput.erase(perturbationVectorInput.begin());
      }
      perturbationVectorInput.push_back(perturbationUsed);
    }
  }

  std::cout << "NOT converged the frequency is w = " << w << std::endl;
  return delta_n_step_one;
}

Eigen::MatrixXcd Sternheimer::AndersonMixing(Eigen::MatrixXcd& inNew,
                                             Eigen::MatrixXcd& inOld,
                                             Eigen::MatrixXcd& outNew,
                                             Eigen::MatrixXcd& outOld,
                                             double alpha) const {

  std::complex<double> beta =
      (outNew - inNew).cwiseProduct(outNew - inNew - outOld + inOld).sum() /
      ((outNew - inNew).cwiseProduct((outOld - inOld))).sum();

  Eigen::MatrixXcd nIn = beta * inOld + (1 - beta) * inNew;
  Eigen::MatrixXcd nOut = beta * outOld + (1 - beta) * outNew;

  return alpha * nOut + (1 - alpha) * nIn;
}

Eigen::MatrixXcd Sternheimer::NPAndersonMixing(
    std::vector<Eigen::MatrixXcd>& Input, std::vector<Eigen::MatrixXcd>& Output,
    double alpha) const {

  Eigen::MatrixXcd DeltaN = Output.back() - Input.back();

  // Building Linear System for Coefficients
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(Input.size() - 1, Input.size() - 1);
  Eigen::VectorXd c = Eigen::VectorXd::Zero(Input.size() - 1);

  for (Index m = 1; m < Input.size(); m++) {

    c(m - 1) = (DeltaN - Output.at(Output.size() - 1 - m) +
                Input.at(Input.size() - 1 - m))
                   .cwiseProduct(DeltaN)
                   .sum()
                   .real();
    for (Index j = 1; j < Input.size(); j++) {

      A(m - 1, j - 1) =
          (DeltaN - Output.at(Output.size() - 1 - m) +
           Input.at(Input.size() - 1 - m))
              .cwiseProduct(DeltaN - Output.at(Output.size() - 1 - j) +
                            Input.at(Input.size() - 1 - j))
              .sum()
              .real();
    }
  }
  // Solving the System to obtain coefficients
  Eigen::VectorXcd coefficients = A.fullPivHouseholderQr().solve(c);

  // Mixing the Potentials

  Eigen::MatrixXcd OutMixed = Output.back();
  Eigen::MatrixXcd InMixed = Input.back();

  for (Index n = 1; n < Input.size(); n++) {

    OutMixed += coefficients(n - 1) * (Output.at(Output.size() - 1 - n) -
                                       Output.at(Output.size() - 1));
    InMixed += coefficients(n - 1) *
               (Input.at(Input.size() - 1 - n) - Input.at(Input.size() - 1));
  }

  // Returning the linear Mix of Input and Output
  return alpha * OutMixed + (1 - alpha) * InMixed;
}

Eigen::MatrixXcd Sternheimer::BroydenMixing(
    std::vector<Eigen::MatrixXcd>& Input, std::vector<Eigen::MatrixXcd>& Output,
    double jacobianScaling) const {

  // Approximation of the first inverse Jacobian
  Eigen::MatrixXcd FirstInverseJacobian =
      jacobianScaling * Eigen::MatrixXcd::Identity(_basis_size, _basis_size);

  Index histSize = Input.size();

  Eigen::MatrixXd gamma = Eigen::MatrixXd::Zero(histSize, histSize);
  Eigen::MatrixXd alpha = Eigen::MatrixXd::Zero(histSize - 1, histSize - 1);
  Eigen::MatrixXd beta = Eigen::MatrixXd::Zero(histSize - 1, histSize - 1);
  Eigen::VectorXd weights = Eigen::VectorXd::Zero(histSize);

  for (Index m = 0; m < histSize; m++) {
    weights(m) = 1 / sqrt(abs((Output.at(m) - Input.at(m))
                                  .cwiseProduct((Output.at(m) - Input.at(m)))
                                  .sum()
                                  .real()));
  }
  for (Index m = 0; m < histSize - 1; m++) {
    for (Index l = 0; l < histSize - 1; l++) {
      alpha(m, l) =
          weights(m) * weights(l) *
          ((Output.at(l + 1) - Input.at(l + 1) - Output.at(l) + Input.at(l)) /
           (Output.at(l + 1) - Input.at(l + 1) - Output.at(l) + Input.at(l))
               .norm())
              .cwiseProduct(((Output.at(m + 1) - Input.at(m + 1) -
                              Output.at(m) + Input.at(m)) /
                             (Output.at(m + 1) - Input.at(m + 1) -
                              Output.at(m) + Input.at(m))
                                 .norm()))
              .sum()
              .real();
    }
  }
  beta = (weights(0) * weights(0) *
              Eigen::MatrixXd::Identity(histSize - 1, histSize - 1) +
          alpha)
             .inverse();
  for (Index m = 0; m < histSize; m++) {
    for (Index l = 0; l < histSize - 1; l++) {
      for (Index k = 0; k < m - 1; k++) {
        gamma(m, l) +=
            weights(k) *
            ((Output.at(k + 1) - Input.at(k + 1) - Output.at(k) + Input.at(k)) /
             (Output.at(k + 1) - Input.at(k + 1) - Output.at(k) + Input.at(k))
                 .norm())
                .cwiseProduct(Input.at(m) - Input.at(m - 1))
                .sum()
                .real() *
            beta(k, l);
      }
    }
  }
  Eigen::MatrixXcd BroydenMix =
      Input.back() + FirstInverseJacobian * (Output.back() - Input.back());
  for (Index n = 0; n < histSize - 1; n++) {

    BroydenMix -=
        weights(n) * gamma(histSize, n) * FirstInverseJacobian *
            ((Output.at(n + 1) - Input.at(n + 1) - Output.at(n) + Input.at(n)) /
             (Output.at(n + 1) - Input.at(n + 1) - Output.at(n) + Input.at(n))
                 .norm()) +
        (Input.at(n + 1) - Input.at(n)) /
            (Output.at(n + 1) - Input.at(n + 1) - Output.at(n) + Input.at(n))
                .norm();
  }

  return BroydenMix;
}

std::vector<Eigen::Matrix3cd> Sternheimer::Polarisability() const {

  std::vector<std::complex<double>> frequency_evaluation_grid = BuildGrid(
      _opt.start_frequency_grid, _opt.end_frequency_grid,
      _opt.number_of_frequency_grid_points, _opt.imaginary_shift_pade_approx);

  std::vector<std::complex<double>> output_grid =
      BuildGrid(_opt.start_frequency_grid, _opt.end_frequency_grid,
                _opt.number_output_grid_points, _opt.lorentzian_broadening);

  std::vector<Eigen::Matrix3cd> Polar;
  std::vector<Eigen::Matrix3cd> Polar_pade;

  for (Index i = 0; i < frequency_evaluation_grid.size(); i++) {
    Polar.push_back(Eigen::Matrix3cd::Zero());
  }

  PadeApprox pade_1;
  PadeApprox pade_4;
  PadeApprox pade_6;
  pade_1.initialize(4 * frequency_evaluation_grid.size());
  pade_4.initialize(4 * frequency_evaluation_grid.size());
  pade_6.initialize(4 * frequency_evaluation_grid.size());

  AOBasis basis = _orbitals.SetupDftBasis();
  AODipole dipole;
  dipole.Fill(basis);
#pragma omp parallel for
  for (Index n = 0; n < frequency_evaluation_grid.size(); n++) {
    for (Index i = 0; i < 3; i++) {
      Eigen::MatrixXcd delta_n =
          DeltaNSC(frequency_evaluation_grid[n],
                   -_opt.perturbation_strength * dipole.Matrix()[i]);
      for (Index j = i; j < 3; j++) {
        Polar[n](i, j) = -(delta_n.cwiseProduct(dipole.Matrix()[j])).sum();
      }
    }
    for (Index i = 2; i < 3; i++) {
      for (Index j = i + 1; j < 3; j++) {
        Polar[n](j, i) = conj(Polar[n](i, j));
      }
    }
  }

  for (Index n = 0; n < Polar.size(); n++) {
    pade_1.addPoint(frequency_evaluation_grid[n], Polar[n](0, 0));
    pade_1.addPoint(conj(frequency_evaluation_grid[n]), conj(Polar[n](0, 0)));
    pade_1.addPoint(-frequency_evaluation_grid[n], conj(Polar[n](0, 0)));
    pade_1.addPoint(-conj(frequency_evaluation_grid[n]), Polar[n](0, 0));

    pade_4.addPoint(frequency_evaluation_grid[n], Polar[n](1, 1));
    pade_4.addPoint(conj(frequency_evaluation_grid[n]), conj(Polar[n](1, 1)));
    pade_4.addPoint(-frequency_evaluation_grid[n], conj(Polar[n](1, 1)));
    pade_4.addPoint(-conj(frequency_evaluation_grid[n]), Polar[n](1, 1));

    pade_6.addPoint(frequency_evaluation_grid[n], Polar[n](2, 2));
    pade_6.addPoint(conj(frequency_evaluation_grid[n]), conj(Polar[n](2, 2)));
    pade_6.addPoint(-frequency_evaluation_grid[n], conj(Polar[n](2, 2)));
    pade_6.addPoint(-conj(frequency_evaluation_grid[n]), Polar[n](2, 2));
  }

  for (std::complex<double> w : output_grid) {
    Polar_pade.push_back(Eigen::Matrix3cd::Zero());
    Polar_pade[Polar_pade.size() - 1](0, 0) = pade_1.evaluatePoint(w);
    Polar_pade[Polar_pade.size() - 1](1, 1) = pade_4.evaluatePoint(w);
    Polar_pade[Polar_pade.size() - 1](2, 2) = pade_6.evaluatePoint(w);
  }
  return Polar_pade;
}
void Sternheimer::printIsotropicAverage(
    std::vector<Eigen::Matrix3cd>& polar) const {
  std::vector<std::complex<double>> grid =
      BuildGrid(_opt.start_frequency_grid, _opt.end_frequency_grid,
                _opt.number_output_grid_points, _opt.lorentzian_broadening);
  std::cout << "\n" << "#Freq (ev) \t polarizability_isotropic_average" << std::endl;
  for (Index i = 0; i < polar.size(); i++) {
    std::cout << real(grid.at(i)) * votca::tools::conv::hrt2ev << "\t"
              << real((polar.at(i)(2, 2))) + real(polar.at(i)(1, 1)) +
                     real(polar.at(i)(0, 0)) / 3
              << std::endl;
  }
}
std::vector<double> Sternheimer::getIsotropicAverage(
    std::vector<Eigen::Matrix3cd>& polar) const {

  std::vector<double> iA;

  for (Index i = 0; i < polar.size(); i++) {

    iA.push_back(real((polar.at(i)(2, 2))) + real(polar.at(i)(1, 1)) +
                 real(polar.at(i)(0, 0)) / 3);
  }
  return iA;
}
std::vector<Eigen::Vector3cd> Sternheimer::EnergyGradient() const {

  QMMolecule mol = _orbitals.QMAtoms();

  // Setting up Grid for Fxc functional

  std::cout << "SetUP ERIS" << std::endl;
  AOBasis dftbasis = _orbitals.SetupDftBasis();
  AOBasis auxbasis = _orbitals.SetupAuxBasis();
  ERIs eris;
  eris.Initialize(dftbasis, auxbasis);

  std::cout << "SetUP Grid" << std::endl;
  Vxc_Grid grid;
  grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                 dftbasis);
  Vxc_Potential<Vxc_Grid> Vxcpot(grid);
  Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());

  // QMMolecule molecule;
  Index number_of_atoms = mol.size();

  std::vector<Eigen::Vector3cd> EnergyGrad;

  AO3ddipole ao3dDipole;
  // Loop over Nuclei

  for (int k = 0; k < number_of_atoms; k++) {

    // ao3dDipole.FillPotential(dftbasis, r);
    std::cout << "Loop over Atom " << k << std::endl;

    ao3dDipole.setCenter(mol.at(k).getPos());
    ao3dDipole.Fill(dftbasis);

    double sign = 1.0;

    EnergyGrad.push_back(Eigen::Vector3d::Zero());

    for (int a = 0; a < 3; a++) {

      Eigen::MatrixXcd DeltaN = DeltaNSC(0.0, sign * ao3dDipole.Matrix()[a]);
      Eigen::MatrixXcd contract = eris.ContractRightIndecesWithMatrix(DeltaN);
     Eigen::MatrixXcd FxcInt = Vxcpot.IntegrateFXC(_density_Matrix, DeltaN);
      Eigen::MatrixXcd DeltaV =
          sign * ao3dDipole.Matrix()[a] + contract + FxcInt;
      EnergyGrad[k][a] = _density_Matrix.transpose()
                             .cwiseProduct(DeltaV * mol.at(k).getNuccharge())
                             .sum();
    }

    for (int l = 0; l < number_of_atoms; l++) {
      if (l != k) {
        Eigen::Vector3d distance = (mol.at(k).getPos() - mol.at(l).getPos());
        EnergyGrad[k] += mol.at(k).getNuccharge() * mol.at(l).getNuccharge() *
                         distance / std::pow(distance.norm(), 3);
      }
    }
  }

  for (int i = 0; i < EnergyGrad.size(); i++) {

    // std::cout<<"Atom Type : "<<mol.at(i).getElement()<<" Atom Index :
    // "<<i<<std::endl; std::cout<<"Gradient =
    // "<<std::endl<<EnergyGrad[i]<<std::endl<<std::endl;
    std::cout << EnergyGrad[i][0].real() << "\t" << EnergyGrad[i][1].real()
              << "\t" << EnergyGrad[i][2].real() << std::endl;
  }

  return EnergyGrad;
}

Eigen::MatrixXcd Sternheimer::GreensFunctionLHS(std::complex<double> w) const {
  return _Hamiltonian_Matrix - w * _overlap_Matrix;
}

Eigen::MatrixXcd Sternheimer::AnalyticGreensfunction(
    std::complex<double> w) const {
  return GreensFunctionLHS(w).colPivHouseholderQr().solve(
      -1 * Eigen::MatrixXcd::Identity(_basis_size, _basis_size));
}
double Sternheimer::Lorentzian(double center, std::complex<double> freq) const {
  double gamma = 0.001;

  // double prefactor = 1/(votca::tools::conv::Pi*gamma);
  double enumerator = gamma * gamma;
  double denominator = (std::pow(abs(freq - center), 2) + gamma * gamma);
  return (enumerator / denominator);
}

Eigen::MatrixXcd Sternheimer::NonAnalyticGreensfunction(
    std::complex<double> freq) const {
  Eigen::MatrixXcd NonAnalyticGreens =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);

  std::complex<double> factor(0, 2 * votca::tools::conv::Pi);

  for (Index v = 0; v < _num_occ_lvls; v++) {
    NonAnalyticGreens += factor * Lorentzian(_mo_energies[v], freq) *
                         _mo_coefficients.col(v) *
                         _mo_coefficients.col(v).transpose();
  }

  return NonAnalyticGreens;
}

Eigen::MatrixXcd Sternheimer::GreensFunction(
    std::complex<double> frequency) const {

  Eigen::MatrixXcd GreensFunctionMatrix =
      NonAnalyticGreensfunction(frequency) + AnalyticGreensfunction(frequency);
  return GreensFunctionMatrix;
}
Eigen::MatrixXcd Sternheimer::CoulombMatrix(Eigen::Vector3d gridpoint) const {

  AOBasis basis = _orbitals.SetupDftBasis();
  AOMultipole aoesp;
  aoesp.FillPotential(basis, gridpoint);
  return aoesp.Matrix();
}

Eigen::MatrixXcd Sternheimer::ScreenedCoulomb(
    Eigen::Vector3d gridpoint1, std::complex<double> frequency) const {

  AOBasis dftbasis = _orbitals.SetupDftBasis();
  AOBasis auxbasis = _orbitals.SetupAuxBasis();
  ERIs eris;
  eris.Initialize(dftbasis, auxbasis);

  Vxc_Grid grid;
  grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                 dftbasis);
  Vxc_Potential<Vxc_Grid> Vxcpot(grid);
  Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());
  Eigen::MatrixXcd coulombmatrix = CoulombMatrix(gridpoint1);
  Eigen::MatrixXcd DeltaN = DeltaNSC(frequency, coulombmatrix);
  Eigen::MatrixXcd contract = eris.ContractRightIndecesWithMatrix(DeltaN);
  Eigen::MatrixXcd FxcInt = Vxcpot.IntegrateFXC(_density_Matrix, DeltaN);
  Eigen::MatrixXcd DeltaV = coulombmatrix + contract + FxcInt;

  return DeltaV;
}

Eigen::VectorXd Sternheimer::EvaluateBasisAtPosition(
    const AOBasis& dftbasis, const Eigen::Vector3d& pos) const {

  // get value of orbitals at each gridpoint
  Eigen::VectorXd tmat = Eigen::VectorXd::Zero(dftbasis.AOBasisSize());
  for (const AOShell& shell : dftbasis) {
    const double decay = shell.getMinDecay();
    const Eigen::Vector3d& shellpos = shell.getPos();
    Eigen::Vector3d dist = shellpos - pos;
    double distsq = dist.squaredNorm();
    // if contribution is smaller than -ln(1e-10), calc density
    if ((decay * distsq) < 20.7) {
      Eigen::VectorBlock<Eigen::VectorXd> tmat_block =
          tmat.segment(shell.getStartIndex(), shell.getNumFunc());
      shell.EvalAOspace(tmat_block, pos);
    }
  }
  return tmat;
}

std::complex<double> Sternheimer::SelfEnergy_at_r(double omega,
                                                  Eigen::Vector3d gridpoint1,
                                                  Index n, Index m) const {

  std::vector<std::complex<double>> w = BuildGrid(
      _opt.start_frequency_grid, _opt.end_frequency_grid,
      _opt.number_of_frequency_grid_points, _opt.imaginary_shift_pade_approx);

  AOBasis basis = _orbitals.SetupDftBasis();
  Eigen::VectorXd chi = EvaluateBasisAtPosition(basis, gridpoint1);
  Eigen::MatrixXd left = chi * chi.transpose();
  Eigen::MatrixXcd right = Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  // Integral over frequencies omega: trapezoidal rule
  double delta = (_opt.end_frequency_grid - _opt.start_frequency_grid) /
                 _opt.number_of_frequency_grid_points;
  Index i = 0;

  for (std::complex<double> omega_p : w) {

    if (i == 0) {
      right += 0.5 * GreensFunction(omega + omega_p) *
               ScreenedCoulomb(gridpoint1, omega_p);
    } else if (i == _opt.number_of_frequency_grid_points - 1) {
      right += 0.5 * GreensFunction(omega + omega_p) *
               ScreenedCoulomb(gridpoint1, omega_p);
    } else {
      right += GreensFunction(omega + omega_p) *
               ScreenedCoulomb(gridpoint1, omega_p);
    }
    i++;
  }

  return _mo_coefficients.col(n).transpose() * left * right *
         _mo_coefficients.col(m);
}

std::complex<double> Sternheimer::SelfEnergy(double omega, Index n,
                                             Index m) const {
  AOBasis basis = _orbitals.SetupDftBasis();
  Vxc_Grid grid;
  grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                 basis);

  std::complex<double> prefactor(0., 1.0 / (2 * std::acos(-1.0)));  // i/2pi

  std::complex<double> corrections(0.0, 0.0);

  //for (Index i = 0; i < grid.getBoxesSize(); ++i) {
  for (Index i = 0; i < 3; ++i) {  
    const GridBox& box = grid[i];
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    for (Index p = 0; p < box.size(); p++) {
      std::cout << p << std::endl;
      corrections += weights[p] * SelfEnergy_at_r(omega, points[p], n, m);
    }
  }
  std::cout << "\n State:\t" << n << " " << m << "\t sigma_c \t" << corrections
            << std::endl;
return corrections;
}

}  // namespace xtp

}  // namespace votca