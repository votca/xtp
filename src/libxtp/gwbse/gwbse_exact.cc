/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#include <votca/ctp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/threecenter.h>
#include <votca/xtp/gwbse_exact.h>

namespace votca {
    namespace xtp {

        bool GWBSE_Exact::Evaluate() {
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Running GWBSE_Exact :D "
                    << std::flush;
            
            return false;
            
            // ****** Check input ******
            
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Molecule Coordinates [A] "
                    << std::flush;

            for (QMAtom* atom : _orbitals.QMAtoms()) {

                std::string output = (boost::format("  %1$s"
                        "   %2$+1.4f %3$+1.4f %4$+1.4f")
                        % (atom->getType())
                        % (atom->getPos().getX() * tools::conv::bohr2ang)
                        % (atom->getPos().getY() * tools::conv::bohr2ang)
                        % (atom->getPos().getZ() * tools::conv::bohr2ang)).str();
                
                CTP_LOG(ctp::logDEBUG, *_pLog)
                        << ctp::TimeStamp()
                        << output
                        << std::flush;
                
            }

            // Check which QC program was used for the DFT run
            // -> Implicit info about MO coefficient storage order
            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " DFT data was created by "
                    << _orbitals.getQMpackage()
                    << std::flush;

            // ****** Configure orbitals object ******
            
            _orbitals.setRPAindices(_rpamin, _rpamax);
            _orbitals.setGWAindices(_qpmin, _qpmax);
            _orbitals.setBSEindices(_bse_vmin, _bse_cmax, _bse_maxeigenvectors);
            
            // ****** Set TDA approximation ******

            // ****** Load DFT basis ******
            
            BasisSet dftbs;
            dftbs.LoadBasisSet(_dftbasis_name);
            _orbitals.setDFTbasis(_dftbasis_name);

            CTP_LOG(ctp::logDEBUG, *_pLog)
                    << ctp::TimeStamp()
                    << " Loaded DFT Basis Set "
                    << _dftbasis_name
                    << std::flush;
            
            // ****** Fill DFT AO basis ******
            
            // ****** Load aux. basis ******
            
            // ****** Fill Exchange-Correlation Potential (VXC) matrix ******
            
            // ****** Fill overlap matrix ******
            
            // ****** Fill Coulomb matrix ******
            
            // ****** Prepare 3-center integral object ******
            
            // ****** Initialize QP energies ******
            
            // ****** Start iterations ******

            return false;
        }

    }
}