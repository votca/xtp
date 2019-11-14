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
#ifndef _VOTCA_XTP_COORD_CONTAINER
#define _VOTCA_XTP_COORD_CONTAINER
#include <array>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>
#include <votca/tools/types.h>
namespace votca {
namespace xtp {

typedef std::array<Index, 2> BondIdx;
typedef std::array<Index, 3> AngleIdx;
typedef std::array<Index, 4> DihedralIdx;

template <std::size_t SIZE>
std::ostream& operator<<(std::ostream& s, const std::array<Index, SIZE>& a) {
  s << "(";
  for (const Index& i : a) s << i << ",";
  s << "\b)";
  return s;
}

template <std::size_t SIZE>
std::array<Index, SIZE> IdxReverse(std::array<Index, SIZE> idx) {
  std::reverse(idx.begin(), idx.end());
  return idx;
}

template <std::size_t SIZE>
bool IdxContains(const Index i, const std::array<Index, SIZE>& a) {
  return (std::find(a.begin(), a.end(), i) != a.end());
}

template <typename IdxType, typename ValType>
struct CoordContainer {
  CoordContainer() = default;

  std::vector<IdxType> _indices;
  std::map<IdxType, ValType> _values;
  votca::Index num;

  ValType& operator[](const IdxType idx) {
    IdxType xdi = IdxReverse(idx);

    if (Contains(idx)) return _values[idx];
    if (Contains(xdi)) return _values[xdi];

    _values[idx] = 0;
    _indices.emplace_back(idx);
    num = _indices.size();

    return _values[idx];
  }

  ValType operator[](const IdxType idx) const { return Get(idx); }

  bool Contains(const IdxType idx) const {
    IdxType xdi = IdxReverse(idx);
    return (_values.find(idx) != _values.end() ||
            _values.find(xdi) != _values.end());
  }

  ValType Get(const IdxType idx) const {
    auto search = _values.find(idx);
    if (search != _values.end()) {
      return search->second;
    }
    throw std::runtime_error("Could not find any value ");
  }

  Index size() const { return Index(_indices.size()); }

  friend std::ostream& operator<<(std::ostream& stream,
                                  const CoordContainer& cc) {
    for (const auto& i : cc._indices) {
      stream << i << " = " << cc[i] << std::endl;
    }
    return stream;
  }
};
}  // namespace xtp
}  // namespace votca

#endif  //_VOTCA_XTP_COORD_CAONTAINER
