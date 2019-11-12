/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

#pragma once
#ifndef VOTCA_XTP_SEGMENTMAPPER_H
#define VOTCA_XTP_SEGMENTMAPPER_H

#include <type_traits>
#include <votca/csg/pdbwriter.h>
#include <votca/tools/property.h>
#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmmolecule.h>
#include <votca/xtp/segid.h>
#include <votca/xtp/topology.h>

namespace votca {
namespace xtp {
template <class AtomContainer>
class SegmentMapper {
 public:
  SegmentMapper(Logger& log);

  void LoadMappingFile(const std::string& mapfile);

  AtomContainer map(const Segment& seg, const SegId& segid) const;

  AtomContainer map(const Segment& seg, QMState state) const;

  AtomContainer map(const Segment& seg, const std::string& coordfilename) const;

 private:
  using mapAtom = typename std::iterator_traits<
      typename AtomContainer::iterator>::value_type;

  typedef std::pair<long, std::string> atom_id;

  struct FragInfo {
    std::vector<double> _weights;
    std::vector<atom_id> _mapatom_ids;
    std::vector<atom_id> _mdatom_ids;
    std::vector<Index> _map_local_frame;
  };

  struct Seginfo {
    std::pair<long, Index> minmax;
    std::vector<Index> mdatoms;
    std::vector<FragInfo> fragments;
    bool map2md;
    std::string segname;
    std::vector<double> weights;
    std::vector<atom_id> mapatoms;
    std::map<std::string, std::string> coordfiles;
  };
  std::map<std::string, std::string> _mapatom_xml;
  std::map<std::string, Seginfo> _segment_info;

  Index FindVectorIndexFromAtomId(
      Index atomid, const std::vector<mapAtom*>& fragment_mapatoms) const;

  void ParseFragment(Seginfo& seginfo, const tools::Property& frag);

  template <typename T>
  Eigen::Vector3d CalcWeightedPos(const std::vector<double>& weights,
                                  const T& atoms) const;

  void PlaceMapAtomonMD(const std::vector<mapAtom*>& fragment_mapatoms,
                        const std::vector<const Atom*>& fragment_mdatoms) const;

  void MapMapAtomonMD(const FragInfo& frag,
                      const std::vector<mapAtom*>& fragment_mapatoms,
                      const std::vector<const Atom*>& fragment_mdatoms) const;

  Logger& _log;
  std::pair<long, Index> CalcAtomIdRange(const Segment& seg) const;
  std::pair<long, Index> CalcAtomIdRange(const std::vector<Index>& seg) const;

  atom_id StringToMapIndex(const std::string& map_string) const;

  atom_id StringToMDIndex(const std::string& md_string) const;

  Index getRank(const mapAtom& atom) const { return atom.getRank(); }

  std::vector<double> getWeights(const tools::Property& frag) const;

  std::string getFrame(const tools::Property& frag) const {
    if (frag.exists(_mapatom_xml.at("frame"))) {
      return frag.get(_mapatom_xml.at("frame")).template as<std::string>();
    }
    return frag.get("localframe").template as<std::string>();
  }

  void FillMap() {
    _mapatom_xml["tag"] = "MP";
    _mapatom_xml["name"] = "MPole";
    _mapatom_xml["atoms"] = "mpoles";
    _mapatom_xml["coords"] = "multipoles";
    _mapatom_xml["weights"] = "mp_weights";
    _mapatom_xml["frame"] = "mp_localframe";
  }
};

template <>
inline void SegmentMapper<QMMolecule>::FillMap() {
  _mapatom_xml["tag"] = "QM";
  _mapatom_xml["name"] = "QMAtom";
  _mapatom_xml["atoms"] = "qmatoms";
  _mapatom_xml["coords"] = "qmcoords";
  _mapatom_xml["weights"] = "qm_weights";
  _mapatom_xml["frame"] = "qm_localframe";
}

template <>
inline Index SegmentMapper<QMMolecule>::getRank(const QMAtom&) const {
  return 0;
}

using QMMapper = SegmentMapper<QMMolecule>;
using StaticMapper = SegmentMapper<StaticSegment>;
using PolarMapper = SegmentMapper<PolarSegment>;
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SEGMENTMAPPER_H
