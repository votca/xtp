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

#include <votca/xtp/atom.h>
#include <votca/xtp/fragment.h>
#include <votca/xtp/molecule.h>
#include <votca/xtp/segment.h>

using namespace std;
using namespace votca::tools;

namespace votca {
namespace xtp {

/// Destructor
Molecule::~Molecule() {

  _segments.clear();
  _fragments.clear();
}

void Molecule::AddSegment(Segment *segment) {
  _segments.push_back(segment);
  segment->setMolecule(this);
}

void Molecule::AddFragment(Fragment *fragment) {
  _fragments.push_back(fragment);
  fragment->setMolecule(this);
}
}  // namespace xtp
}  // namespace votca
