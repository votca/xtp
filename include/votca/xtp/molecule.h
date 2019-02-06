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

#ifndef VOTCA_XTP_MOLECULE_H
#define VOTCA_XTP_MOLECULE_H

#include "atom.h"
#include <map>
#include <string>
#include <vector>
#include <votca/csg/basemolecule.h>
#include <votca/tools/vec.h>

namespace CSG = votca::csg;

namespace votca {
namespace xtp {

class Topology;
class Fragment;
class Segment;

/**
 * @brief The xtp molecule class is essentially a wrapper for csg base molecule
 *
 * This class facilitates the mapping atoms, note that the inheritance is
 * protected, this is intentenial because the base molecule specifically deals
 * with atoms and not beads. This is of particular importance when making use
 * of the structural methods available to base molecule.
 *
 * E.g. the csg base molecule considers all beads when running structural
 * analysis methods, this is because the beads are likely to be pseudo-atoms
 * (as in the hydrogen atoms may be combined with other atoms).
 *
 * However, the xtp molecule class explicity treates hydrogen atoms, thus to get
 * meaninful data from xtp molecule class when dealing with organic molecules
 * the hydrogen atoms must first be removed from the structure, so that only the
 * backbone is considered.
 */
class Molecule : protected CSG::BaseMolecule<Atom> {
 public:
  Molecule(int id, std::string name) {
    CSG::BaseMolecule<Atom>::id_.setId(id);
    CSG::BaseMolecule<Atom>::name_.setName(name);
  }
  Molecule() {}
  ~Molecule();

  int getId() { return CSG::BaseMolecule<Atom>::id_.getId(); }
  const std::string &getName() { return CSG::BaseMolecule<Atom>::getName(); }

  void AddSegment(Segment *segment);
  void AddFragment(Fragment *fragment);
  void AddAtom(Atom *atom) { CSG::BaseMolecule<Atom>::AddBead(atom); }

  std::vector<Fragment *> &Fragments() { return _fragments; }
  std::vector<Segment *> &Segments() { return _segments; }

  Atom *getAtom(const int &id) { return CSG::BaseMolecule<Atom>::getBead(id); }
  size_t AtomCount() { return CSG::BaseMolecule<Atom>::BeadCount(); }

  inline void setTopology(Topology *container) { _top = container; }
  Topology *getTopology() { return _top; }

 private:
  Topology *_top;

  std::vector<Segment *> _segments;
  std::vector<Fragment *> _fragments;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MOLECULE_H
