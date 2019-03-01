#ifndef _VOTCA_XTP_COORD_CONTAINER
#define _VOTCA_XTP_COORD_CONTAINER
#include <array>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

// typedef std::tuple<int,int> Bond;
// typedef std::tuple<int,int,int> Angle;
// typedef std::tuple<int,int,int,int> Dihedral;

#define IND_TYPE int

typedef std::array<IND_TYPE, 2> BondIdx;
typedef std::array<IND_TYPE, 3> AngleIdx;
typedef std::array<IND_TYPE, 4> DihedralIdx;

template <std::size_t SIZE>
std::ostream& operator<<(std::ostream& s, const std::array<IND_TYPE, SIZE>& a) {
  s << "(";
  for (const IND_TYPE& i : a) s << i << ",";
  s << "\b)";
  return s;
}

template <std::size_t SIZE>
std::array<IND_TYPE, SIZE> IdxReverse(std::array<IND_TYPE, SIZE> idx) {
  /* std::array<IND_TYPE, SIZE> xdi; */
  /* for (IND_TYPE i = 0; i < idx.size(); ++i){ */
  /*     xdi[i] = idx[idx.size()-1-i]; */
  /* } */
  /* return xdi; */
  std::cout << "reverse" << std::endl;
  std::cout << idx << std::endl;
  std::reverse(idx.begin(), idx.end());
  std::cout << idx << std::endl;
  return idx;
}

template <std::size_t SIZE>
bool IdxContains(const IND_TYPE i, const std::array<IND_TYPE, SIZE> a) {
  return (std::find(a.begin(), a.end(), i) != a.end());
}

template <typename IdxType, typename ValType>
struct CoordContainer {
  CoordContainer(){};

  std::vector<IdxType> _indices;
  std::map<IdxType, ValType> _values;
  size_t num;

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
    std::stringstream m;
    m << "Could not find any value " << std::endl;
    throw std::runtime_error(m.str());
  }

  size_t size() const { return _indices.size(); }

  friend std::ostream& operator<<(std::ostream& stream,
                                  const CoordContainer& cc) {
    for (const auto& i : cc._indices) {
      stream << i << " = " << cc[i] << std::endl;
    }
  }
};

#endif  //_VOTCA_XTP_COORD_CAONTAINER
