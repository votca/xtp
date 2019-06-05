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
#include "../../include/votca/xtp/units.h"

namespace votca {                                                               
  namespace xtp {                                                                 

    const tools::DistanceUnit Units::distance_unit =                             
      tools::DistanceUnit::bohr;                                                  
    const tools::MassUnit Units::mass_unit =                                     
      tools::MassUnit::atomic_mass_units;                                         
    const tools::TimeUnit Units::time_unit = tools::TimeUnit::seconds;           
    const tools::ChargeUnit Units::charge_unit = tools::ChargeUnit::e;           
    const tools::EnergyUnit Units::energy_unit =                                 
      tools::EnergyUnit::hartrees; 

  }
}
