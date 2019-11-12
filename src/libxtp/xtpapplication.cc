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

#include <votca/tools/globals.h>
#include <votca/tools/propertyiomanipulator.h>

#include <boost/format.hpp>
#include <votca/xtp/version.h>
#include <votca/xtp/xtpapplication.h>

namespace votca {
namespace xtp {

XtpApplication::XtpApplication() { ; }

/**
 * \brief Adds program options to the executable
 *
 * Every executable requires option file for calculators it is running
 * It is thus a part of the base XtpApplication class
 *
 */
void XtpApplication::Initialize(void) {
  AddProgramOptions()("options,o", boost::program_options::value<std::string>(),
                      "  calculator options");
}

void XtpApplication::ShowHelpText(std::ostream &out) {
  std::string name = ProgramName();
  if (VersionString() != "") {
    name = name + ", version " + VersionString();
  }
  xtp::HelpTextHeader(name);
  HelpText(out);
  out << "\n\n" << VisibleOptions() << std::endl;
}

}  // namespace xtp
}  // namespace votca
