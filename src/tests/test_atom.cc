/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE atom_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/atom.h>

using namespace votca::tools;
using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(atom_test)

BOOST_AUTO_TEST_CASE(constructors_test) {
  Atom atom1(1, "C");

  Molecule* molecule_ptr = nullptr;
  string residue_name = "residue_1";
  int residue_number = 11;
  string molecular_dynamics_atom_name = "Carbon";
  int molecular_dynamics_atom_id = 2;
  bool has_quantum_chemistry_part = false;
  int quantum_chemistry_id = 3;
  vec quantum_chemistry_position(4.0, 5.0, 6.0);
  string element = "C";
  double weight = 12.01;
  Atom atom2(molecule_ptr, residue_name, residue_number,
             molecular_dynamics_atom_name, molecular_dynamics_atom_id,
             has_quantum_chemistry_part, quantum_chemistry_id,
             quantum_chemistry_position, element, weight);
}

BOOST_AUTO_TEST_CASE(getters_test) {
  {
    Atom atom1(3, "H");
    BOOST_CHECK_EQUAL(atom1.getId(), 3);
    BOOST_CHECK_EQUAL(atom1.getName(), "H");
  }
  {
    Molecule* molecule_ptr = nullptr;
    string residue_name = "residue_1";
    int residue_number = 11;
    string molecular_dynamics_atom_name = "Carbon";
    int molecular_dynamics_atom_id = 2;
    bool has_quantum_chemistry_part = false;
    int quantum_chemistry_id = 3;
    vec quantum_chemistry_position(4.0, 5.0, 6.0);
    string element = "C";
    double weight = 12.01;

    Atom atom2(molecule_ptr, residue_name, residue_number,
               molecular_dynamics_atom_name, molecular_dynamics_atom_id,
               has_quantum_chemistry_part, quantum_chemistry_id,
               quantum_chemistry_position, element, weight);

    BOOST_CHECK_EQUAL(atom2.getId(), 2);
    BOOST_CHECK_EQUAL(atom2.getName(), molecular_dynamics_atom_name);
    BOOST_CHECK_EQUAL(atom2.getResnr(), residue_number);
    BOOST_CHECK_EQUAL(atom2.getResname(), residue_name);
    BOOST_CHECK_EQUAL(atom2.getMass(), weight);
    BOOST_CHECK_EQUAL(atom2.getElement(), element);
    BOOST_CHECK_EQUAL(atom2.getQMId(), quantum_chemistry_id);
    BOOST_CHECK_EQUAL(atom2.getQMPos(), quantum_chemistry_position);
    BOOST_CHECK(!atom2.HasQMPart());
  }
}

BOOST_AUTO_TEST_CASE(setters_test) {
  {
    string residue_name = "residue_1";
    int residue_number = 11;
    int quantum_chemistry_id = 3;
    vec quantum_chemistry_position(4.0, 5.0, 6.0);
    string element = "C";
    double weight = 12.01;
    vec position(7.0, 8.0, 9.0);
    bool has_position = true;

    Atom atom;
    atom.setResnr(residue_number);
    atom.setResname(residue_name);
    atom.setMass(weight);
    atom.setQMPart(quantum_chemistry_id, quantum_chemistry_position);
    atom.setElement(element);
    atom.setPos(position);
    atom.HasPos(has_position);

    atom.setTopology(nullptr);
    atom.setMolecule(nullptr);
    atom.setSegment(nullptr);
    atom.setFragment(nullptr);

    BOOST_CHECK_EQUAL(atom.getResnr(), residue_number);
    BOOST_CHECK_EQUAL(atom.getResname(), residue_name);
    BOOST_CHECK_EQUAL(atom.getMass(), weight);
    BOOST_CHECK_EQUAL(atom.getElement(), element);
    BOOST_CHECK_EQUAL(atom.getQMId(), quantum_chemistry_id);
    BOOST_CHECK_EQUAL(atom.getQMPos(), quantum_chemistry_position);
    BOOST_CHECK_EQUAL(atom.getPos(), position);
    BOOST_CHECK(atom.HasPos());
    BOOST_CHECK(atom.HasQMPart());

    BOOST_CHECK_EQUAL(atom.getTopology(), nullptr);
    BOOST_CHECK_EQUAL(atom.getMolecule(), nullptr);
    BOOST_CHECK_EQUAL(atom.getSegment(), nullptr);
    BOOST_CHECK_EQUAL(atom.getFragment(), nullptr);
  }
}

BOOST_AUTO_TEST_CASE(test_translateby) {
  Atom atom(1, "Ne");

  vec position(1.0, 2.0, 3.0);
  atom.setPos(position);

  vec translation(1.0, 1.0, 1.0);
  atom.TranslateBy(translation);
  BOOST_CHECK_EQUAL(atom.getPos(), (position + translation));
}

BOOST_AUTO_TEST_SUITE_END()
