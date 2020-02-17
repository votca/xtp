/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

#ifndef __XTP_SETTINGS__H
#define __XTP_SETTINGS__H

#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <votca/tools/property.h>

using votca::tools::Property;

/*
 * \brief Hierarchical representation of a QM Package input
 */
namespace votca {
namespace xtp {

class Settings {
 public:
  // Decompose a Property object into Settings
  Settings() = default;
  ~Settings() = default;

  /**
   * \brief Transform Properties into settings
   * @param Property object
   */
  void read_property(const Property& prop);

  /**
   * \brief read Settings from xml
   * @param path to the xml
   */
  void load_from_xml(const std::string& path);

  /**
   * \brief Merge two settings objects
   * @param other Settings object
   */
  void merge(const Settings& other);

  /**
   * \brief Get a given key
   * @param key
   */
  template <typename T = std::string>
  T get(const std::string& key) const {
    auto it = this->_nodes.find(key);
    if (it != this->_nodes.end()) {
      return it->second.get(key).as<T>();
    } else {
      std::ostringstream oss;
      oss << "Unknown keyword: " << key << "\n";
      throw std::runtime_error(oss.str());
    }
  }

  /**
   * \brief Check that the input is correct
   */
  void validate() const;

  friend std::ostream& operator<<(std::ostream& os, const Settings& sett);

 protected:
  using Settings_map = std::unordered_map<std::string, Property>;
  Settings_map _nodes;
  std::string _executable_path;

  Settings_map::const_iterator search_for_mandatory_keyword(
      const std::string& key) const;

  std::vector<std::string> _general_properties = {
      "auxbasisset",            // string
      "basisset",               // string
      "charge",                 // int
      "convergence_tightness",  // std::string
      "external_charge",        // Eigen::Vector9d
      "optimize",               // boolean
      "polarisation",           // boolean
      "pseudopotential",        // string
      "spin"                    // int
  };

  std::vector<std::string> _mandatory_keyword = {
      "executable",  // string
      "functional",  // string
      "name",        // string, one of: orca, nwchem, gaussian, xtpdft,
  };

  std::unordered_map<std::string, std::vector<std::string>> _keyword_options{
      {"name", {"orca", "gaussian", "nwchem", "xtpdft"}}};
};

}  // namespace xtp
}  // namespace votca
#endif
