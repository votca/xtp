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

#pragma once
#ifndef VOTCA_XTP_ATOMCONTAINER_H
#define VOTCA_XTP_ATOMCONTAINER_H
#include <limits>
#include <typeinfo>
#include <votca/tools/elements.h>
#include <votca/xtp/checkpoint.h>
#include <votca/xtp/eigen.h>

/**
 * \brief Basic Container for QMAtoms,PolarSites and Atoms
 *
 *
 *
 */

namespace votca {
namespace xtp {

template <class T>
class AtomContainer {
 public:
  AtomContainer(std::string type, int id) : _type(type), _id(id){};

  AtomContainer(CheckpointReader& r) { this->ReadFromCpt(r); }
  virtual ~AtomContainer(){};

  typedef typename std::vector<T>::iterator iterator;

  const std::string& getType() const { return _type; }

  void setType(std::string type) { _type = type; }

  int getId() const { return _id; }

  int size() const { return _atomlist.size(); }

  void push_back(const T& atom) {
    _atomlist.push_back(atom);
    calcPos();
  }
  void push_back(T&& atom) {
    _atomlist.push_back(atom);
    calcPos();
  }

  void AddContainer(const AtomContainer<T>& container) {
    _type += "_" + container._type;
    _atomlist.insert(_atomlist.end(), container._atomlist.begin(),
                     container._atomlist.end());
    calcPos();
  }

  const T& at(int index) const { return _atomlist.at(index); }
  T& at(int index) { return _atomlist.at(index); }

  const T& operator[](int index) const { return _atomlist[index]; }
  T& operator[](int index) { return _atomlist[index]; }

  typename std::vector<T>::iterator begin() { return _atomlist.begin(); }
  typename std::vector<T>::iterator end() { return _atomlist.end(); }

  typename std::vector<T>::const_iterator begin() const {
    return _atomlist.begin();
  }
  typename std::vector<T>::const_iterator end() const {
    return _atomlist.end();
  }

  const Eigen::Vector3d& getPos() const { return _pos; }

  // calculates the lowest and highest point in the cube, sorrounding the
  // molecule
  std::pair<Eigen::Vector3d, Eigen::Vector3d> CalcSpatialMinMax() const {
    std::pair<Eigen::Vector3d, Eigen::Vector3d> result;
    Eigen::Vector3d min =
        std::numeric_limits<double>::max() * Eigen::Vector3d::Ones();
    Eigen::Vector3d max =
        std::numeric_limits<double>::min() * Eigen::Vector3d::Ones();
    for (const T& atom : _atomlist) {
      const Eigen::Vector3d& pos = atom.getPos();
      if (pos.x() < min.x()) min.x() = pos.x();
      if (pos.x() > max.x()) max.x() = pos.x();
      if (pos.y() < min.y()) min.y() = pos.y();
      if (pos.y() > max.y()) max.y() = pos.y();
      if (pos.z() < min.z()) min.z() = pos.z();
      if (pos.z() > max.z()) max.z() = pos.z();
    }
    result.first = min;
    result.second = max;
    return result;
  }

  std::vector<std::string> FindUniqueElements() const {
    std::vector<std::string> result;
    for (const T& atom : _atomlist) {
      if (std::find(result.begin(), result.end(), atom.getElement()) ==
          result.end()) {
        result.push_back(atom.getElement());
      }
    }
    return result;
  }

  void Translate(const Eigen::Vector3d& shift) {
    for (T& atom : _atomlist) {
      atom.Translate(shift);
    }
    _pos += shift;
  }

  void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& ref_pos) {
    for (T& atom : _atomlist) {
      atom.Rotate(R, ref_pos);
    }
    calcPos();
  }

  virtual void WriteToCpt(CheckpointWriter& w) const {
    w(_type, "type");
    w(_id, "id");
    w(int(_atomlist.size()), "size");
    T element(0, "H", Eigen::Vector3d::Zero());
    CptTable table =
        w.openTable(element.identify() + "s", element, _atomlist.size());
    std::vector<typename T::data> dataVec(_atomlist.size());
    for (std::size_t i = 0; i < _atomlist.size(); ++i) {
      _atomlist[i].WriteData(dataVec[i]);
    }

    table.write(dataVec);
  }
  virtual void ReadFromCpt(CheckpointReader& r) {
    r(_type, "type");
    r(_id, "id");
    int size = 0;
    r(size, "size");
    if (size == 0) {
      return;
    }
    T element(0, "H", Eigen::Vector3d::Zero());  // dummy element to get
                                                 // .identify for type
    CptTable table = r.openTable(element.identify() + "s", element);
    _atomlist.clear();
    _atomlist.reserve(table.numRows());
    std::vector<typename T::data> dataVec(table.numRows());
    table.read(dataVec);
    for (std::size_t i = 0; i < table.numRows(); ++i) {
      _atomlist.push_back(T(dataVec[i]));
    }
    calcPos();
  }

  void calcPos() {
    tools::Elements element;
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    double totalmass = 0.0;
    for (const T& atom : _atomlist) {
      double mass = element.getMass(atom.getElement());
      totalmass += mass;
      pos += mass * atom.getPos();
    }
    _pos = pos / totalmass;
  }

 protected:
  std::vector<T> _atomlist;
  std::string _type;
  int _id;

 private:
  Eigen::Vector3d _pos = Eigen::Vector3d::Zero();
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ATOMCONTAINER_H
