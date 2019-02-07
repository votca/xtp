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


#include <votca/xtp/qmmolecule.h>
#include <votca/tools/elements.h>
#include <votca/xtp/checkpointwriter.h>
#include <votca/xtp/checkpointreader.h>

using namespace std;
using namespace votca::tools;

namespace votca { namespace xtp {

  
      void QMMolecule::WriteXYZ(const std::string& filename, std::string header) const{

          std::ofstream out(filename);
          if (!out.is_open()) {
                throw std::runtime_error("Bad file handle: " + filename);
            }
          out<<this->_atomlist.size()<<endl;
          out<<header<<endl;
          for (const QMAtom& atom:_atomlist) {
                const Eigen::Vector3d pos = atom.getPos() * tools::conv::bohr2ang;
                out<<atom.getElement()<<" "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
          }
          out.close();
          return;
        }


     void QMMolecule::LoadFromXYZ(const std::string& filename) {

            string line;
            std::ifstream in;
            string type;
            in.open(filename.c_str(), std::ios::in);
            if (!in) throw runtime_error(string("Error reading coordinates from: ")
                    + filename);
            int atomCount = 0;
            std::getline(in, line);

            Tokenizer tok1(line," \t");
            std::vector<std::string> line1;
            tok1.ToVector(line1);
            if(line1.size()!=1){
              throw std::runtime_error("First line of xyz file should contain number of atoms, nothing else.");
            }
            std::getline(in, line);//Comment line

            if (in.is_open()) {
                while (in.good()) {
                    std::getline(in, line);

                    vector< string > split;
                    Tokenizer toker(line, " \t");
                    toker.ToVector(split);
                    if(split.size()<4){continue;}
                    // Interesting information written here: e.g. 'C 0.000 0.000 0.000'
                    string element = split[0];
                    double x = stod(split[1]);
                    double y = stod(split[2]);
                    double z = stod(split[3]);
                    Eigen::Vector3d pos = {x, y, z};
                    _atomlist.push_back(QMAtom(atomCount, element, pos * tools::conv::ang2bohr));
                    atomCount++;
                }
            } else {
                throw std::runtime_error("No such file: '" + filename + "'.");
            }
            return;
        }

    

}}


