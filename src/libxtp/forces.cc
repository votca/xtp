/*
 *            Copyright 2009-2017 The VOTCA Development Team
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

#include <votca/tools/linalg.h>
#include <votca/xtp/forces.h>
#include <boost/format.hpp>

namespace votca {
    namespace xtp {

        void Forces::Initialize(Property *options) {

            // checking if there is only one segment
            _nsegments = _segments.size();
            if (_nsegments > 1) throw runtime_error(string("\n Force calculation for more than 1 conjugated segment not supported. Stopping!"));

            // pre-check forces method
            std::vector<string> choices = {"forward", "central"};
            _force_method = options->ifExistsAndinListReturnElseThrowRuntimeError<string>(".method", choices);

            // output level
            _noisy_output = options->ifExistsReturnElseReturnDefault<bool>(".noisy", false); 
            
            
            // precaution in case we implement approx. analytic forces in the future
            if ((_force_method == "forward") || (_force_method == "central")) {
                _displacement = options->ifExistsReturnElseReturnDefault<double>(".displacement", 0.001); // Angstrom
            }

            // check for force removal options
            choices = {"total", "CoM", "none"};
            string _force_removal = options->ifExistsAndinListReturnElseThrowRuntimeError<string>(".removal", choices);
            if (_force_removal == "total") _remove_total_force = true;
            if (_force_removal == "CoM") _remove_CoM_force = true;

            _natoms = _segments[0]->Atoms().size();
            _forces = ub::zero_matrix<double>(_natoms, 3);


            return;

        }

        void Forces::Calculate(const double& energy) {

            ctp::TLogLevel _ReportLevel = _pLog->getReportLevel(); // backup report level
            if ( ! _noisy_output ){
                _pLog->setReportLevel(ctp::logERROR); // go silent for force calculations
            }
            
            //backup current coordinates (WHY?)
            std::vector <ctp::Segment* > _molecule;
            ctp::Segment _current_coordinates(0, "mol");
            _qminterface.Orbitals2Segment(&_current_coordinates, _orbitals);
            _molecule.push_back(&_current_coordinates);

            // displace all atoms in each Cartesian coordinate and get new energy
            std::vector< ctp::Atom* > _atoms;
            std::vector< ctp::Atom* > ::iterator ait;
            _atoms = _current_coordinates.Atoms();

            int _i_atom = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

                if ( _noisy_output ){
                    CTP_LOG(ctp::logINFO, *_pLog) << "FORCES--DEBUG working on atom " << _i_atom << flush;
                }
                
                // Calculate Force on this atom
                ub::matrix_range< ub::matrix<double> > _force = ub::subrange(_forces, _i_atom, _i_atom + 1, 0, 3);
                if (_force_method == "forward") NumForceForward(energy, ait, _force, _molecule);
                if (_force_method == "central") NumForceCentral(energy, ait, _force, _molecule);
                _i_atom++;
            }

            _pLog->setReportLevel(_ReportLevel); // 

            // Remove Total Force, if requested
            if (_remove_total_force) RemoveTotalForce();
            //Report();

            return;
        }

        void Forces::Report() {

            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   ---- FORCES (Hartree/Bohr)   ")).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("        %1$s differences   ") % _force_method).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("        displacement %1$1.4f Angstrom   ") % _displacement).str() << flush;
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("   Atom\t x\t  y\t  z ")).str() << flush;

            for (unsigned _i = 0; _i < _forces.size1(); _i++) {
                CTP_LOG(ctp::logINFO, *_pLog) << (boost::format(" %1$4d    %2$+1.4f  %3$+1.4f  %4$+1.4f")
                        % _i % _forces(_i, 0) % _forces(_i, 1) % _forces(_i, 2)).str() << flush;
            }

            return;

        }

        /* Calculate forces on an atom numerically by forward differences */
        void Forces::NumForceForward(double energy, std::vector< ctp::Atom* > ::iterator ait, ub::matrix_range< ub::matrix<double> >& _force,
                std::vector<ctp::Segment*> _molecule) {

            // get this atoms's current coordinates
            vec _current_pos = (*ait)->getQMPos(); // in nm

            for (unsigned _i_cart = 0; _i_cart < 3; _i_cart++) {

                // get displacement std::vector
                vec _displaced(0, 0, 0);
                if (_i_cart == 0) {
                    _displaced.setX(_displacement * tools::conv::ang2nm); // x, _displacement in Angstrom, now in nm
                }
                if (_i_cart == 1) {
                    _displaced.setY(_displacement * tools::conv::ang2nm); // y, _displacement in in Angstrom, now in nm
                }
                if (_i_cart == 2) {
                    _displaced.setZ(_displacement * tools::conv::ang2nm); // z, _displacement in in Angstrom, now in nm
                }

                // update the coordinate
                vec _pos_displaced = _current_pos + _displaced;

                (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment

                // run DFT and GW-BSE for this geometry
                _gwbse_engine.ExcitationEnergies(_qmpackage, _molecule, _orbitals);

                // get total energy for this excited state
                double energy_displaced = _orbitals->GetTotalEnergy(_spin_type, _opt_state);

                // calculate force and put into matrix
                _force(0, _i_cart) = (energy - energy_displaced) / (_displacement * votca::tools::conv::ang2bohr); // force a.u./a.u.
                (*ait)->setQMPos(_current_pos); // restore original coordinate into segment
            } // Cartesian directions
            return;
        }

        /* Calculate forces on atoms numerically by central differences */
        void Forces::NumForceCentral(double energy, std::vector< ctp::Atom* > ::iterator ait, ub::matrix_range< ub::matrix<double> >& _force,
                std::vector<ctp::Segment*> _molecule) {


            vec _current_pos = (*ait)->getQMPos(); // in nm

            // go through all cartesian components
            for (unsigned _i_cart = 0; _i_cart < 3; _i_cart++) {

                if ( _noisy_output ){
                    CTP_LOG(ctp::logINFO, *_pLog) << "FORCES--DEBUG           Cartesian component " << _i_cart << flush;
                }
                
                // get displacement vector in positive direction
                vec _displaced(0, 0, 0);
                if (_i_cart == 0) {
                    _displaced.setX(_displacement * tools::conv::ang2nm); // x, _displacement in Bohr
                }
                if (_i_cart == 1) {
                    _displaced.setY(_displacement * tools::conv::ang2nm); // y, _displacement in in Angstrom
                }
                if (_i_cart == 2) {
                    _displaced.setZ(_displacement * tools::conv::ang2nm); // z, _displacement in in Angstrom
                }

                // update the coordinate
                vec _pos_displaced = _current_pos + _displaced;
                (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment

                // run DFT and GW-BSE for this geometry
                _gwbse_engine.ExcitationEnergies(_qmpackage, _molecule, _orbitals);

                // get total energy for this excited state
                double energy_displaced_plus = _orbitals->GetTotalEnergy(_spin_type, _opt_state);

                // get displacement vector in negative direction

                // update the coordinate
                _pos_displaced = _current_pos - 2.0 * _displaced;
                (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment

                // run DFT and GW-BSE for this geometry
                _gwbse_engine.ExcitationEnergies(_qmpackage, _molecule, _orbitals);

                // get total energy for this excited state
                double energy_displaced_minus = _orbitals->GetTotalEnergy(_spin_type, _opt_state);

                // calculate force and put into matrix
                _force(0, _i_cart) = 0.5 * (energy_displaced_minus - energy_displaced_plus) / (_displacement * votca::tools::conv::ang2bohr); // force a.u./a.u.

                (*ait)->setQMPos(_current_pos); // restore original coordinate into segment
            }

            return;
        }

        /* Adjust forces so that sum of forces is zero */
        void Forces::RemoveTotalForce() {

            // total force on all atoms
            ub::vector<double> _total_force = TotalForce();

            // zero total force
            for (unsigned _i_atom = 0; _i_atom < _natoms; _i_atom++) {
                for (unsigned _i_cart = 0; _i_cart < 3; _i_cart++) {
                    _forces(_i_atom, _i_cart) -= _total_force(_i_cart) / _natoms;
                }
            }
            return;
        }

        /* Determine Total Force on all atoms */
        ub::vector<double> Forces::TotalForce() {

            ub::vector<double> _total_force(3, 0.0);
            for (unsigned _i_atom = 0; _i_atom < _natoms; _i_atom++) {
                for (unsigned _i_cart = 0; _i_cart < 3; _i_cart++) {
                    _total_force(_i_cart) += _forces(_i_atom, _i_cart);
                }
            }
            return _total_force;
        }
    }
}
