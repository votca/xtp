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
/// For earlier commit history see ctp commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef VOTCA_XTP_ATOM_H
#define VOTCA_XTP_ATOM_H

#include <exception>
#include <map>
#include <string>
#include <votca/csg/basebead.h>
#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>

namespace CSG = votca::csg;

namespace votca {
namespace xtp {

class Topology;
class Segment;
class Fragment;

/**
 * @brief Atom class stores info relevant to the structure from Molecual
 * Dynamics
 *
 * The atom class is a child of base bead
 */
class Atom : public CSG::BaseBead {
 public:
  Atom(int molecule_id, std::string residue_name, int resnr,
       std::string md_atom_name, int md_atom_id, bool hasQMPart, int qm_atom_id,
       tools::vec qmPos, std::string element, double mass)
      : _resnr(resnr),
        _resname(residue_name),
        _hasQM(hasQMPart),
        _qmId(qm_atom_id),
        _qmPos(qmPos),
        _element(element) {

    molecule_id_.setId(molecule_id);
    setMass(mass);
    HasPos(false);
    setId(md_atom_id);
    setName(md_atom_name);
  }

  /**
   * @brief Atom constructor basic
   *
   * @param[in] atom_id
   * @param[in] atom_name
   */
  Atom(int atom_id, std::string atom_name) : _hasQM(false), _qmId(-1) {

    setId(atom_id);
    setName(atom_name);
  }

  Atom(Atom *stencil)
      : _top(nullptr),
        _resnr(stencil->getResnr()),
        _resname(stencil->getResname()),
        _hasQM(stencil->HasQMPart()),
        _qmId(stencil->getQMId()),
        _qmPos(stencil->getQMPos()),
        _element(stencil->getElement()) {

    setId(stencil->getId());
    setName(stencil->getName() + "_ghost");
    setMass(stencil->getMass());
    setPos(stencil->getPos());
    HasPos(true);
  }

  Atom(){};
  ~Atom() { _Q.clear(); }

  const int &getResnr() const { return _resnr; }

  inline void setTopology(Topology *container) { _top = container; }
  inline void setSegment(Segment *container) { _seg = container; }
  inline void setFragment(Fragment *container) { _frag = container; }

  Topology *getTopology() { return _top; }
  Segment *getSegment() { return _seg; }
  Fragment *getFragment() { return _frag; }

  inline void setResnr(const int &resnr) { _resnr = resnr; }
  inline void setResname(const std::string &resname) { _resname = resname; }
  inline void setWeight(const double &weight) {
    throw std::runtime_error(
        "ERROR atom setWeight is now depricated use setMass.");
  }
  inline void setQMPart(const int &qmid, tools::vec qmPos);
  inline void setQMPos(const tools::vec &qmPos) { _qmPos = qmPos; }
  inline void setElement(const std::string &element) { _element = element; }
  inline void TranslateBy(const tools::vec &shift) {
    CSG::BaseBead::bead_position_ = CSG::BaseBead::bead_position_ + shift;
  }

  inline const int &getResnr() { return _resnr; }
  inline const std::string &getResname() { return _resname; }
  inline const double &getWeight() {
    throw std::runtime_error("getWeight is now depricated use getMass.");
  }
  inline const int &getQMId() { return _qmId; }
  inline const tools::vec &getQMPos() { return _qmPos; }
  inline const std::string &getElement() { return _element; }

  inline const double &getQ(int state) { return _Q.at(state); }
  inline const double &getQ() { return _q->second; }
  inline void setQ(std::map<int, double> Q) { _Q = Q; }
  void chrg(int state) { _q = _Q.find(state); }

  inline void setPTensor(tools::matrix &ptensor) { _ptensor = ptensor; }
  const tools::matrix &getPTensor() { return _ptensor; }

  bool HasQMPart() { return _hasQM; }

 protected:
  Topology *_top;
  Segment *_seg;
  Fragment *_frag;

  int _resnr;
  std::string _resname;

  bool _hasQM;
  int _qmId;
  tools::vec _qmPos;
  std::string _element;

  // charge state of segment => partial charge
  std::map<int, double> _Q;
  std::map<int, double>::iterator _q;
  tools::matrix _ptensor;
};

inline void Atom::setQMPart(const int &qmid, tools::vec qmPos) {
  if (qmid > -1) {
    _hasQM = true;
    _qmId = qmid;
    _qmPos = qmPos;
  } else {
    _hasQM = false;
    _qmId = -1;
  }
}
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ATOM_H
