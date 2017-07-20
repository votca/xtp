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

// Overload of uBLAS prod function with MKL/GSL implementations


#include <votca/xtp/gwbseengine.h>
#include <votca/xtp/gwbse.h>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/linalg.h>

#include <votca/ctp/logger.h>



using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;
        // +++++++++++++++++++++++++++++ //
        // GWBSEENGINE MEMBER FUNCTIONS  //
        // +++++++++++++++++++++++++++++ //

        void GWBSEENGINE::Initialize(Property* options, string _archive_filename) {


            _archive_file = _archive_filename;
            string key = Identify();

            // get the tasks
            string _tasks_string = options->get(".tasks").as<string> ();
            _do_guess = false;
            _do_dft_input = false;
            _do_dft_run = false;
            _do_dft_parse = false;
            _do_gwbse = false;

            if (_tasks_string.find("guess") != std::string::npos) _do_guess = true;
            if (_tasks_string.find("input") != std::string::npos) _do_dft_input = true;
            if (_tasks_string.find("dft") != std::string::npos)   _do_dft_run = true;
            if (_tasks_string.find("parse") != std::string::npos) _do_dft_parse = true;
            if (_tasks_string.find("gwbse") != std::string::npos) _do_gwbse = true;

            // XML option file for GWBSE
            string _gwbse_xml = options->get(".gwbse_options").as<string> ();
            load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());

            // DFT log and MO file names
            _MO_file = options->get(".mofile").as<string> ();
            _dftlog_file = options->get(".dftlog").as<string> ();

            // Logger redirection
            _redirect_logger = options->ifExistsReturnElseReturnDefault<bool>(".redirect_logger", false);
            _logger_file = "gwbse.log";
            
            // for requested merged guess, two archived orbitals objects are needed
            if ( _do_guess ){
                _guess_archiveA = options->ifExistsReturnElseThrowRuntimeError<string>(".archiveA");
                _guess_archiveB = options->ifExistsReturnElseThrowRuntimeError<string>(".archiveB");
            }

            return;
        }

        /* 
         *    CALL DFT and GWBSE modules to get excitation energies
         * 
         */


        void GWBSEENGINE::ExcitationEnergies(QMPackage* _qmpackage, vector<ctp::Segment*> _segments, Orbitals* _orbitals) {


            //redirect log, if required
            // define own logger for GW-BSE that is written into a runFolder logfile
            ctp::Logger _gwbse_engine_logger(_pLog->getReportLevel());
            if (_redirect_logger) {
                _gwbse_engine_logger.setMultithreading(false);
                _gwbse_engine_logger.setPreface(ctp::logINFO, (format("\n ...")).str());
                _gwbse_engine_logger.setPreface(ctp::logERROR, (format("\n ...")).str());
                _gwbse_engine_logger.setPreface(ctp::logWARNING, (format("\n ...")).str());
                _gwbse_engine_logger.setPreface(ctp::logDEBUG, (format("\n ...")).str());
                _qmpackage->setLog(&_gwbse_engine_logger);
            }

            if (_do_dft_input) {

                // required for merged guess
                Orbitals *_orbitalsAB = NULL;
                if (_qmpackage->GuessRequested() && _do_guess) { // do not want to do an SCF loop for a dimer
                    if (_redirect_logger) {
                       CTP_LOG(ctp::logINFO, _gwbse_engine_logger) << "Guess requested, reading molecular orbitals" << flush;
                    } else {
                       CTP_LOG(ctp::logINFO, *_pLog) << "Guess requested, reading molecular orbitals" << flush;
                    }
                    Orbitals _orbitalsA, _orbitalsB;
                    _orbitalsAB = new Orbitals();
                    // load the corresponding monomer orbitals and prepare the dimer guess

                    // failed to load; wrap-up and finish current job
                    if (!_orbitalsA.Load(_guess_archiveA)) {
                        throw runtime_error(string("Do input: failed loading orbitals from ") + _guess_archiveA);
                    }

                    if (!_orbitalsB.Load(_guess_archiveB)) {
                        throw runtime_error(string("Do input: failed loading orbitals from ") + _guess_archiveB);

                    }

                    _orbitals->PrepareGuess(&_orbitalsA, &_orbitalsB, _orbitalsAB);

                }
                
                _qmpackage->WriteInputFile(_segments, _orbitalsAB);
            }

            if (_do_dft_run) {

                bool run_success = _qmpackage->Run();
                if (!run_success) {
                    throw runtime_error(string("\n GW-BSE without DFT is difficult. Stopping!"));
                }
            }

            // parse DFT data, if required
            if (_do_dft_parse) {
                if (_redirect_logger) {
                    CTP_LOG(ctp::logINFO, _gwbse_engine_logger) << "Parsing DFT data from " << _dftlog_file << " and " << _MO_file << flush;
                } else {
                    CTP_LOG(ctp::logINFO, *_pLog) << "Parsing DFT data from " << _dftlog_file << " and " << _MO_file << flush;
                }
                _qmpackage->setLogFileName(_dftlog_file);
                _qmpackage->ParseLogFile(_orbitals);
                _qmpackage->setOrbitalsFileName(_MO_file);
                _qmpackage->ParseOrbitalsFile(_orbitals);
                _orbitals->setDFTbasis(_qmpackage->getBasisSetName());
            }

            // if no parsing of DFT data is requested, reload serialized orbitals object
            if (!_do_dft_parse && _do_gwbse) {
                if (_redirect_logger) {
                    CTP_LOG(ctp::logINFO, _gwbse_engine_logger) << "Loading serialized data from " << _archive_file << flush;
                } else {
                    CTP_LOG(ctp::logINFO, *_pLog) << "Loading serialized data from " << _archive_file << flush;
                }
                _orbitals->Load(_archive_file);
            }

            if (_do_gwbse) {
                GWBSE _gwbse = GWBSE(_orbitals);
                _gwbse.setLogger(_pLog);
                if (_redirect_logger) _gwbse.setLogger(&_gwbse_engine_logger);
                _gwbse.Initialize(&_gwbse_options);
                _gwbse.Evaluate();
                if (_redirect_logger) SaveRedirectedLogger(&_gwbse_engine_logger);
            }
            return;
        }

        /* Saves a redirected logger (run output) to file */
        void GWBSEENGINE::SaveRedirectedLogger(ctp::Logger* pLog) {

            // write logger to log file
            ofstream ofs;
            //string gwbse_logfile = "gwbse.log";
            ofs.open(_logger_file.c_str(), ofstream::out);
            if (!ofs.is_open()) {
                throw runtime_error("Bad file handle: " + _logger_file);
            }
            ofs << (*pLog) << endl;
            ofs.close();
            return;
        }



    }
}
