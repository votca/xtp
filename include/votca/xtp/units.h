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
                                                                                
#pragma once                                                                    
#ifndef VOTCA_XTP_UNITS_H                                                    
#define VOTCA_XTP_UNITS_H   

#include <votca/tools/unitconverter.h>
namespace votca {
  namespace xtp {

    class Units {
      public:
        static const tools::DistanceUnit distance_unit;                               
        static const tools::MassUnit mass_unit;                                       
        static const tools::TimeUnit time_unit;                                       
        static const tools::ChargeUnit charge_unit;                                   
        static const tools::EnergyUnit energy_unit;                                   
    };
  }
}
#endif // VOTCA_XTP_UNITS_H   
