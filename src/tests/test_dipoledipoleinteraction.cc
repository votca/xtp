
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE dipoledipoleinteraction_test

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <votca/xtp/dipoledipoleinteraction.h>
#include <votca/xtp/eigen.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(dipoledipoleinteraction_test)

BOOST_AUTO_TEST_CASE(dipoledipoleinteraction_test) {

  PolarSite zero(0, "H", Eigen::Vector3d::UnitX());
  PolarSegment seg0("zero", 0);
  seg0.push_back(zero);
  PolarSite one(1, "C", Eigen::Vector3d::UnitZ());
  PolarSegment seg1("one", 1);
  seg1.push_back(one);

  std::vector<PolarSegment> segs;
  segs.push_back(seg0);
  segs.push_back(seg1);
  double tholedamp = 0.39;
  eeInteractor interactor(tholedamp);

  DipoleDipoleInteraction dipdip(interactor, segs);

  // building reference
  Eigen::MatrixXd ref = Eigen::MatrixXd::Zero(6, 6);
  PolarSegment seg_ref("ref", 100);
  seg_ref.push_back(zero);
  seg_ref.push_back(one);
  for (int i = 1; i < seg_ref.size(); i++) {
    for (int j = 0; j < i; j++) {
      ref.block<3, 3>(3 * i, 3 * j) =
          interactor.FillTholeInteraction_diponly(seg_ref[i], seg_ref[j]);
      ref.block<3, 3>(3 * j, 3 * i) =
          interactor.FillTholeInteraction_diponly(seg_ref[j], seg_ref[i]);
    }
  }

  for (int i = 0; i < seg_ref.size(); i++) {
    ref.block<3, 3>(3 * i, 3 * i) = seg_ref[i].getPInv();
  }
  // building matrix via (i,j) operator
  Eigen::MatrixXd elementwise = Eigen::MatrixXd::Zero(6, 6);
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      elementwise(i, j) = dipdip(i, j);
    }
  }
  bool op_check = elementwise.isApprox(ref, 1e-6);
  if (!op_check) {
    std::cout << "ref" << std::endl;
    std::cout << ref << std::endl;
    std::cout << "operator" << std::endl;
    std::cout << elementwise << std::endl;
  }

  // building matrix via matrix vector product
  Eigen::MatrixXd gemv = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd ident = Eigen::MatrixXd::Identity(6, 6);
  for (int i = 0; i < 6; i++) {
    gemv.col(i) = dipdip * ident.col(i);
  }
  bool gemv_check = gemv.isApprox(ref, 1e-6);
  if (!gemv_check) {
    std::cout << "ref" << std::endl;
    std::cout << ref << std::endl;
    std::cout << "gemv" << std::endl;
    std::cout << gemv << std::endl;
  }
  // building matrix via iterator product
  Eigen::MatrixXd iterator = Eigen::MatrixXd::Zero(6, 6);
  for (int k = 0; k < dipdip.outerSize(); ++k) {
    for (DipoleDipoleInteraction::InnerIterator it(dipdip, k); it; ++it) {
      iterator(it.row(), k) = it.value();
    }
  }
  bool iterator_check = iterator.isApprox(ref, 1e-6);
  if (!iterator_check) {
    std::cout << "ref" << std::endl;
    std::cout << ref << std::endl;
    std::cout << "iterator" << std::endl;
    std::cout << gemv << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()