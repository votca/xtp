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
#include <boost/format.hpp>
#include <vector>

#include "votca/xtp/polarsegment.h"

using namespace std;
using namespace votca::tools;

namespace votca { namespace xtp {


//MPS files have a weird format positions can be in bohr or angstroem,
//multipoles are in q*bohr^k, with k rank of multipole and polarisabilities are in angstroem^3
void PolarSegment::LoadFromMPS(const std::string& filename){

    std::string line;
    std::ifstream intt;
    intt.open(filename.c_str());
    double unit_conversion=tools::conv::ang2bohr;

    int readinmultipoles=0;
    int numberofmultipoles=0;
    Eigen::VectorXd multipoles;

    if (!intt.is_open() ) {
        throw runtime_error("File:"+filename+" could not be opened");
    }
    while ( intt.good() ) {

        std::getline(intt, line);
        vector<string> split;
        Tokenizer toker(line, " \t");
        toker.ToVector(split);

        if ( !split.size()      ||
              split[0] == "!"   ||
              split[0].substr(0,1) == "!" ) { continue; }

// ! Interesting information here, e.g.
// ! DCV2T opt
// ! SP        RB3LYP          6-311+G(d,p)
// Units bohr
//
// C          -4.2414603400   -3.8124751600    0.0017575736    Rank  2
//  -0.3853409355
//  -0.0002321905   0.2401559510   0.6602334308
//  -0.7220625314   0.0004894995  -0.0003833545   0.4526409813  -0.50937399
//  P 1.75


        // Units used
        if ( split[0] == "Units") {
            std::string units = split[1];
            if (units != "bohr" && units != "angstrom") {
                throw std::runtime_error( "Unit " + units + " in file "
                                        + filename + " not supported.");
            }
            if (units == "bohr") {
                unit_conversion=1.0;
            }
        }else if ( split.size() == 6 ) {
            // element,  position,  rank limit convert to bohr
            string name = split[0];
            Eigen::Vector3d pos;
            int id=_atomlist.size();
            pos[0] = stod(split[1]);
            pos[1] = stod(split[2]);
            pos[2] = stod(split[3]);
            int rank = boost::lexical_cast<int>(split[5]);
            numberofmultipoles=(rank+1)*(rank+1);
            multipoles=Eigen::VectorXd::Zero(numberofmultipoles);
            pos*=unit_conversion;
            _atomlist.push_back(PolarSite(id,name,pos));
            }
        // 'P', dipole polarizability
        else if ( split[0] == "P") {
            Eigen::Matrix3d p1;
            // Angstroem to bohr
            double pxx, pxy, pxz;
            double      pyy, pyz;
            double           pzz;
            if (split.size() == 7) {
                pxx = stod(split[1]);
                pxy = stod(split[2]);
                pxz = stod(split[3]);
                pyy = stod(split[4]);
                pyz = stod(split[5]);
                pzz = stod(split[6]);
                p1<<pxx,pxy,pxz,
                    pxy,pyy,pyz,
                    pxz,pyz,pzz;
            }else if (split.size() == 2) {
                pxx = stod(split[1]);
                p1=pxx*Eigen::Matrix3d::Identity();;
            }
            else {
                throw std::runtime_error("Invalid line in " + filename
                                         + ": " + line);
            }
            double unit_conversion_3=std::pow(tools::conv::ang2bohr,3);
            p1=p1*unit_conversion_3;
            _atomlist.back().setPolarisation(p1);
        }
        // Multipole lines
        else {
            // stay in bohr
            for (unsigned i = 0; i < split.size(); i++) {
                double qXYZ = stod(split[i]);
                if(multipoles.size()<readinmultipoles){
                    throw std::runtime_error("ReadMpsFile: File"+filename+"is not properly formatted");
                }
                multipoles(readinmultipoles)=qXYZ;
                readinmultipoles++;
            }
            if(readinmultipoles==numberofmultipoles){
                _atomlist.back().setMultipole(multipoles);
                multipoles.resize(0);
                readinmultipoles=0;
            }
        }
    } 
}


double PolarSegment::CalcTotalQ()const{
    double Q=0;
    for(const PolarSite& site:_atomlist){
        Q+=site.getCharge();
    }
    return Q;
}

Eigen::Vector3d PolarSegment::CalcDipole()const{
    Eigen::Vector3d dipole=Eigen::Vector3d::Zero();

    Eigen::Vector3d CoM=this->getPos();
    for(const PolarSite& site:_atomlist){
            dipole += (site.getPos()- CoM) * site.getCharge();
            dipole+=site.getDipole();
    }
    return dipole;
}

void PolarSegment::WriteMPS(const std::string& filename, std::string header) const{

    std::ofstream ofs;
    ofs.open(filename.c_str(), ofstream::out);
    if (!ofs.is_open()) {
        throw runtime_error("Bad file handle: " + filename);
    }

    ofs << (boost::format("! GENERATED BY VOTCA::XTP::%1$s\n") % header);
    ofs << (boost::format("! N=%2$d Q[e]=%1$+1.7f\n")
        % CalcTotalQ() % size());
    ofs << boost::format("Units angstrom\n");

    for(const PolarSite& site:_atomlist){
        ofs <<site.WriteMpsLine("angstrom");
    }
    ofs.close();

}

    

    

}}


