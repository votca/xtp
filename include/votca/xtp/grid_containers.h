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

#ifndef __XTP_GRID_CONTAINERS__H
#define	__XTP_GRID_CONTAINERS__H

#include <votca/tools/vec.h>
#include <votca/xtp/aobasis.h>



namespace votca { namespace xtp {

    
    
    

        class GridContainers {
        public: 
            
            GridContainers(){};
            // containers for radial grids per element
            struct radial_grid {
                std::vector<double> radius;
                std::vector<double> weight;
            };
       
            std::map<std::string,radial_grid> _radial_grids;

            // containers for spherical grids on a unit sphere per element
            struct spherical_grid{
                std::vector<double> theta;
                std::vector<double> phi;
                std::vector<double> weight;
                
            };
            
            std::map<std::string,spherical_grid> _spherical_grids;
            
            // container for cartesian grid points and weights
            struct integration_grid {
                tools::vec grid_pos;//bohr
                double grid_weight;
            };
            
            
            

        };

    }}
#endif	/* NUMERICAL_INTEGRATION_H */