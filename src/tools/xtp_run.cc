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

#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/xtp/sqlapplication.h>
#include <votca/xtp/calculatorfactory.h>
#include <votca/ctp/calculatorfactory.h>


using namespace std;
using namespace votca;

class XtpRun : public xtp::SqlApplication {
public:

  string ProgramName() {
    return "xtp_run";
  }

  void HelpText(ostream &out) {
    out << "Runs excitation/charge transport calculators" << endl;
  }

  void HelpText() {
  };

  void Initialize();
  bool EvaluateOptions();

private:

  //void    PrintDescription(string name, HelpOutputType _help_output_type);

};

namespace propt = boost::program_options;

void XtpRun::Initialize() {
  xtp::Calculatorfactory::RegisterAll();
  ctp::Calculatorfactory::RegisterAll();
  xtp::SqlApplication::Initialize();

  AddProgramOptions("Calculators") ("execute,e", propt::value<string>(),
          "List of calculators separated by ',' or ' '");
  AddProgramOptions("Calculators") ("list,l",
          "Lists all available calculators");
  AddProgramOptions("Calculators") ("description,d", propt::value<string>(),
          "Short description of a calculator");
}

bool XtpRun::EvaluateOptions() {

  string helpdir = "xtp/xml";
  string ctphelpdir = "ctp/xml";
  if (OptionsMap().count("list")) {
    cout << "Available XTP calculators: \n";
    for (xtp::Calculatorfactory::assoc_map::const_iterator iter =
            xtp::Calculators().getObjects().begin();
            iter != xtp::Calculators().getObjects().end(); ++iter) {
      PrintDescription(std::cout, (iter->first), helpdir, Application::HelpShort);
    }
    cout << "Available (wrapped) CTP calculators: \n";
    for (ctp::Calculatorfactory::assoc_map::const_iterator iter =
            ctp::Calculators().getObjects().begin();
            iter != ctp::Calculators().getObjects().end(); ++iter) {
      bool printctp = true;
      std::string ctpcalc = (iter->first).c_str();
      for (xtp::Calculatorfactory::assoc_map::const_iterator xter =
              xtp::Calculators().getObjects().begin();
              xter != xtp::Calculators().getObjects().end(); ++xter) {
        if (ctpcalc.compare((xter->first).c_str()) == 0) {
          printctp = false;
          break;
        }
      }


      if (printctp) {
        PrintDescription(std::cout, (iter->first), ctphelpdir, Application::HelpShort);
      }

    }
    StopExecution();
    return true;
  }


  if (OptionsMap().count("description")) {
    CheckRequired("description", "no calculator is given");
    Tokenizer tok(OptionsMap()["description"].as<string>(), " ,\n\t");
    // loop over the names in the description string
    for (Tokenizer::iterator n = tok.begin(); n != tok.end(); ++n) {
      // loop over calculators
      bool printerror = true;
      for (xtp::Calculatorfactory::assoc_map::const_iterator iter = xtp::Calculators().getObjects().begin();
              iter != xtp::Calculators().getObjects().end(); ++iter) {

        if ((*n).compare((iter->first).c_str()) == 0) {
          PrintDescription(std::cout, (iter->first), helpdir, Application::HelpLong);
          printerror = false;
          break;
        }
      }
      for (ctp::Calculatorfactory::assoc_map::const_iterator iter = ctp::Calculators().getObjects().begin();
              iter != ctp::Calculators().getObjects().end(); ++iter) {

        if ((*n).compare((iter->first).c_str()) == 0) {
          bool printctp = true;
          std::string ctpcalc = (iter->first).c_str();
          for (xtp::Calculatorfactory::assoc_map::const_iterator xter =
                  xtp::Calculators().getObjects().begin();
                  xter != xtp::Calculators().getObjects().end(); ++xter) {
            if (ctpcalc.compare((xter->first).c_str()) == 0) {
              printctp = false;
              break;
            }
          }
          if (printctp) {
            PrintDescription(std::cout, iter->first, "ctp/xml", Application::HelpLong);
            printerror = false;
            break;
          }
        }
      }
      if (printerror) cout << "Calculator " << *n << " does not exist\n";
    }
    StopExecution();
    return true;
  }

  xtp::SqlApplication::EvaluateOptions();
  CheckRequired("options", "Please provide an xml file with calculator options");
  CheckRequired("execute", "Nothing to do here: Abort.");

  Tokenizer calcs(OptionsMap()["execute"].as<string>(), " ,\n\t");
  Tokenizer::iterator it;
  for (it = calcs.begin(); it != calcs.end(); it++) {
    bool _found_calc = false;
    for (xtp::Calculatorfactory::assoc_map::const_iterator iter = xtp::Calculators().getObjects().begin();
            iter != xtp::Calculators().getObjects().end(); ++iter) {

      if ((*it).compare((iter->first).c_str()) == 0) {
        cout << " This is a XTP app" << endl;
        xtp::SqlApplication::AddCalculator(xtp::Calculators().Create((*it).c_str()));
        _found_calc = true;
      }
    }

    if (!_found_calc) {
      for (ctp::Calculatorfactory::assoc_map::const_iterator iter = ctp::Calculators().getObjects().begin();
              iter != ctp::Calculators().getObjects().end(); ++iter) {

        if ((*it).compare((iter->first).c_str()) == 0) {
          cout << " This is a CTP app" << endl;
          xtp::SqlApplication::AddCalculator(ctp::Calculators().Create((*it).c_str()));
        }
      }
    }

    load_property_from_xml(_options, _op_vm["options"].as<string>());
  }
    return true;
  }

  int main(int argc, char** argv) {

    XtpRun xtprun;
    return xtprun.Exec(argc, argv);

  }
