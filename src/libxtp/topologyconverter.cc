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

#include "../../include/votca/xtp/topologyconverter.h"

namespace votca {
namespace xtp {

inline ::votca::csg::Topology TopologyConverter::Convert(
    const Topology& from_xtp_top) const {

  ::votca::csg::Topology to_csg_top;
  /// Convert the xtp topology box boundaries to the csg boundaries
  /// by editing the units
  Eigen::Matrix3d box = from_xtp_top.getBox() *
                        converter_.convert(Topology::units::distance_unit,
                                           csg::Topology::units::distance_unit);
  to_csg_top.setBox(box);

  to_csg_top.setStep(from_xtp_top.getStep());
  to_csg_top.setTime(from_xtp_top.getTime() *
                     converter_.convert(Topology::units::time_unit,
                                        csg::Topology::units::time_unit));

  for (const typename Topology::container_t& xtp_cont : from_xtp_top) {
    tools::StructureParameters params = xtp_cont.getParameters();
    params.set(tools::StructureParameter::MoleculeId,
               params.get<int>(tools::StructureParameter::AtomContainerId));
    params.set(
        tools::StructureParameter::MoleculeType,
        params.get<std::string>(tools::StructureParameter::AtomContainerType));
    /// Set the molecule id and type equal to the atom container id and type
    typename ::votca::csg::Topology::container_t& mol =
        to_csg_top.CreateMolecule(params);
    for (const typename Topology::bead_t& bead : xtp_cont) {

      tools::StructureParameters params_bead = bead.getParameters();
      Eigen::Vector3d pos = params_bead.get<Eigen::Vector3d>(
          tools::StructureParameter::XTP_Position);
      pos *= converter_.convert(Topology::units::distance_unit,
                                csg::Topology::units::distance_unit);
      params_bead.set(tools::StructureParameter::CSG_Position, pos);
      typename ::votca::csg::Topology::bead_t& csg_bead =
          to_csg_top.CreateBead(params);
      mol.AddBead(csg_bead);
    }
  }

  return to_csg_top;
}

inline Topology TopologyConverter::Convert(
    const ::votca::csg::Topology& from_csg_top) const {

  Topology to_xtp_top;
  Eigen::Matrix3d box = from_csg_top.getBox() *
                        converter_.convert(csg::Topology::units::distance_unit,
                                           Topology::units::distance_unit);
  to_xtp_top.setBox(box);

  to_xtp_top.setStep(from_csg_top.getStep());
  to_xtp_top.setTime(from_csg_top.getTime() *
                     converter_.convert(csg::Topology::units::time_unit,
                                        Topology::units::time_unit));

  std::vector<std::pair<int, std::string>> id_and_types;
  std::unordered_map<int, std::vector<tools::StructureParameters>>
      mol_id_and_bead_params;
  /// We must collect the molecules and beads in order because the
  /// order that they are created in xtp is important.
  for (const typename ::votca::csg::Topology::container_t& csg_cont :
       from_csg_top) {
    tools::StructureParameters params = csg_cont.getParameters();
    /// Set the molecule id and type equal to the atom container id and type
    int mol_id = params.get<int>(tools::StructureParameter::MoleculeId);
    std::string molecule_type =
        params.get<std::string>(tools::StructureParameter::MoleculeType);
    id_and_types.push_back(std::pair<int, std::string>(mol_id, molecule_type));

    std::vector<int> bead_ids = csg_cont.getBeadIds();
    for (const int& bead_id : bead_ids) {
      const typename ::votca::csg::Topology::bead_t& bead =
          csg_cont.getBead(bead_id);
      tools::StructureParameters params_bead = bead.getParameters();
      mol_id_and_bead_params[mol_id].push_back(params_bead);
    }
  }

  std::sort(id_and_types.begin(), id_and_types.end());
  assert(id_and_types.at(0).first == 0 &&
         "The smallest id of a molecule in the csg topology object is not 0, "
         "this means when converting to xtp topology segment the ids of the "
         "segments will not match those of the csg molecules.");
  for (std::pair<int, std::string>& id_and_type : id_and_types) {
    to_xtp_top.AddSegment(id_and_type.second);

    typename Topology::container_t seg =
        to_xtp_top.getSegment(id_and_type.first);
    for (tools::StructureParameters bead_params :
         mol_id_and_bead_params[id_and_type.first]) {
      Eigen::Vector3d pos = bead_params.get<Eigen::Vector3d>(
          tools::StructureParameter::CSG_Position);
      pos *= converter_.convert(csg::Topology::units::distance_unit,
                                Topology::units::distance_unit);
      bead_params.set(tools::StructureParameter::XTP_Position, pos);
      typename Topology::bead_t bead(bead_params);
      seg.push_back(bead);
    }
  }
  return to_xtp_top;
}
}  // namespace xtp
}  // namespace votca
