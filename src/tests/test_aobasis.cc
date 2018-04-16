/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE aobasis_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/orbitals.h>
#include <fstream>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(aobasis_test)

BOOST_AUTO_TEST_CASE(FillNormBasis_test) {   
std::ofstream basisfile("notnormalized.xml");
basisfile << "<basis name=\"def2-TZVP\">"<< std::endl;
basisfile << "  <element name=\"Al\">" << std::endl;
basisfile << "    <shell scale=\"1.0\" type=\"D\">"<< std::endl;
basisfile << "      <constant decay=\"1.570000e+00\">"<< std::endl;
basisfile << "        <contractions factor=\"2.000000e-01\" type=\"D\"/>"<<std::endl;
basisfile << "      </constant>"<<std::endl;
basisfile << "      <constant decay=\"3.330000e-01\">"<<std::endl;
basisfile << "        <contractions factor=\"1.000000e+00\" type=\"D\"/>"<<std::endl;
basisfile << "      </constant>"<<std::endl;
basisfile << "    </shell> "<<std::endl;
basisfile << "  </element> "<<std::endl;
basisfile << "</basis> "<<std::endl;
basisfile.close();
   
std::ofstream xyzfile("Al.xyz");
xyzfile << " 1" << std::endl;
xyzfile << " Al" << std::endl;
xyzfile << " Al            .000000     .000000     .000000" << std::endl;
xyzfile.close();

Orbitals orbitals;
orbitals.LoadFromXYZ("Al.xyz");
BasisSet basis;
basis.LoadBasisSet("notnormalized.xml");
AOBasis aobasis;
aobasis.AOBasisFill(&basis,orbitals.QMAtoms());

const AOShell* shell=aobasis.getShell(0);
AOShell::GaussianIterator it;
std::vector<double> ref_results={0.1831079647,0.9155398233};
int i=0;
bool check_norm=true;
for(it=shell->firstGaussian();it<shell->lastGaussian();++it){
  if(std::abs(ref_results[i]-it->getContraction()[2])>1e-7){
   check_norm=false;
   break;
  }
  i++;
}


i=0;
if(!check_norm){
  for(it=shell->firstGaussian();it<shell->lastGaussian();++it){
  std::cout<<"Ref:"<<ref_results[i]<<" result:"<<it->getContraction()[2]<<std::endl;
   i++;
  }
 
}
BOOST_CHECK_EQUAL(check_norm, 1);

}


BOOST_AUTO_TEST_SUITE_END()
