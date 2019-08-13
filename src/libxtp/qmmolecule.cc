/*
 *            Copyright 2016 The MUSCET Development Team
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
#include <votca/csg/io/pdbwriter.h>
#include <votca/csg/io/xyzreader.h>
#include <votca/csg/io/xyzwriter.h>
#include <votca/csg/topology.h>
#include <votca/tools/elements.h>
#include <votca/xtp/checkpointreader.h>
#include <votca/xtp/checkpointwriter.h>
#include <votca/xtp/qmmolecule.h>

using namespace std;
using namespace votca::tools;

namespace votca {
namespace xtp {

/*const tools::DistanceUnit QMMolecule::distance_unit =
    tools::DistanceUnit::bohr;
const tools::MassUnit QMMolecule::mass_unit =
    tools::MassUnit::atomic_mass_units;
const tools::TimeUnit QMMolecule::time_unit = tools::TimeUnit::seconds;
const tools::ChargeUnit QMMolecule::charge_unit = tools::ChargeUnit::e;
const tools::EnergyUnit QMMolecule::energy_unit =
    tools::EnergyUnit::hartrees; */

tools::StructureParameters QMMolecule::getParameters() const {
  tools::StructureParameters params;
  params.set(tools::StructureParameter::AtomContainerType, getName());
  params.set(tools::StructureParameter::AtomContainerId, getId());
  return params;
}

void QMMolecule::WriteXYZ(std::string filename, std::string header) const {
  csg::XYZWriter<csg::Topology> writer;
  writer.Open(filename, false);
  writer.WriteHeader(header, size());
  TopologyContainerConverter converter;
  csg::Topology csg_top = converter.Convert(*this);
  // writer.Write(&csg_top);
  writer.Close();
  return;
}

void QMMolecule::LoadFromFile(std::string filename) {
  csg::XYZReader<csg::Topology> reader;
  csg::Topology csg_top;
  reader.ReadTopology(filename, &csg_top);
  TopologyContainerConverter cont_converter;

  // Preserve name and id
  int temp_id = _id;
  string temp_name = _name;
  *this = cont_converter.Convert<QMMolecule>(csg_top);
  _id = temp_id;
  _name = temp_name;
}
}  // namespace xtp
}  // namespace votca
