/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef VOTCA_XTP_MOLECULE_H
#define VOTCA_XTP_MOLECULE_H

#include <string>
#include <vector>
#include <map>
#include <votca/tools/vec.h>

namespace votca { namespace xtp {

class Topology;
class Atom;
class Fragment;
class Segment;

class Molecule {
public:

    Molecule(int id, std::string name) : _id(id), _name(name) {}
    Molecule() { }
   ~Molecule();

    const int       &getId();
    const std::string    &getName();

    void AddSegment( Segment* segment );
    void AddFragment( Fragment* fragment);
    void AddAtom( Atom* atom);

    std::vector< Atom* >     &Atoms() { return _atoms; }
    std::vector< Fragment* > &Fragments() { return _fragments; }
    std::vector< Segment* >  &Segments() { return _segments; }

    Atom           *getAtom(const int &id);
    const std::string   &getAtomType(const int &id);
    const votca::tools::vec       getAtomPosition(const int &id);
    int             NumberOfAtoms();

    inline void setTopology(Topology *container) { _top = container; }
    Topology   *getTopology() { return _top; }

    /// Load molecule coordinates from a file
    void ReadXYZ ( std::string filename );
    void WritePDB( FILE *out );
    
private:

    Topology *_top;

    std::vector < Segment* >   _segments;
    std::vector < Fragment* >  _fragments;
    std::vector < Atom* >      _atoms ;

    int     _id;
    std::string  _name ;

    std::map<std::string, Atom* > _map_AtomName_Atom;

};

}}

#endif // VOTCA_XTP_MOLECULE_H
