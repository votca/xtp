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
#ifndef VOTCA_XTP_FILEADAPTER_H
#define VOTCA_XTP_FILEADAPTER_H

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>
#include <unordered_map>
#include <votca/tools/unitconverter.h>
#include <votca/tools/structureparameters.h>

#include <Eigen/Dense>

namespace votca {
  namespace xtp {


    /**
     * @brief File Topology Adapter
     *
     * This class is designed to interface with the topology structures in xtp
     * and in csg such that the xtp classes can use the csg topolology and 
     * trajectory input and output readers and writers
     */
    class FileAdapter {

      public: 

        template<class XTP_Topology_T, class CSG_Topology_T>
          void convertXTPTopologyToCSGTopology(XTP_Topology_T & from_xtp_top, CSG_Topology_T & to_csg_top);

        template<class XTP_Container_T, class CSG_Topology_T>
          void convertXTPContainerToCSGTopology(XTP_Container_T & from_xtp_cont, CSG_Topology_T & to_csg_top);

        template<class CSG_Topology_T, class Container_T>  
          void convertCSGTopologyToXTPContainer(CSG_Topology_T & from_csg_top, Container_T & to_cont);

        template<class CSG_Topology_T, class Container_T>  
          void convertCSGTopologyToXTPContainer(CSG_Topology_T & from_csg_top, Container_T & to_cont,const int & csg_molecule_id);

        template<class CSG_Topology_T, class XTP_Topology_T>  
          void convertCSGTopologyToXTPTopology(CSG_Topology_T & from_csg_top, XTP_Topology_T & to_xtp_top);

      private:

        tools::UnitConverter converter_;

    };

    /*****************************************************************************
     * Converting from XTP to CSG 
     *****************************************************************************/

    template<class XTP_Topology_T, class CSG_Topology_T>
      inline void FileAdapter::convertXTPTopologyToCSGTopology(XTP_Topology_T & from_xtp_top, CSG_Topology_T & to_csg_top){

        /// Convert the xtp topology box boundaries to the csg boundaries
        /// by editing the units
        Eigen::Matrix3d box = from_xtp_top.getBox() * 
          converter_.convert(XTP_Topology_T::distance_unit,CSG_Topology_T::distance_unit); 
        to_csg_top.setBox(box);

        to_csg_top.setStep(from_xtp_top.getStep());
        to_csg_top.setTime(from_xtp_top.getTime() * 
            converter_.convert(XTP_Topology_T::time_unit,CSG_Topology_T::time_unit));

        for( typename XTP_Topology_T::container_t & xtp_cont : from_xtp_top){
          tools::StructureParameters params = xtp_cont.getParameters();
          params.set(tools::StructureParameter::MoleculeId,params.get<int>(tools::StructureParameter::AtomContainerId));
          params.set(tools::StructureParameter::MoleculeType,params.get<std::string>(tools::StructureParameter::AtomContainerType));
          /// Set the molecule id and type equal to the atom container id and type
          typename CSG_Topology_T::container_t  mol = to_csg_top.CreateMolecule(params);  
          for( typename XTP_Topology_T::bead_t & bead : xtp_cont ){

            tools::StructureParameters params_bead = bead.getParameters();
            Eigen::VectorXd pos = params_bead.get<Eigen::VectorXd>(tools::StructureParameter::XTP_Position);
            pos*= converter_.convert(XTP_Topology_T::distance_unit,CSG_Topology_T::distance_unit);
            params_bead.set(tools::StructureParameter::CSG_Position,pos);
            typename CSG_Topology_T::bead_t csg_bead = to_csg_top.CreateBead(params);
            mol.AddBead(&csg_bead);
          }
        }
      }

    template<class XTP_Container_T, class CSG_Topology_T>
      inline void FileAdapter::convertXTPContainerToCSGTopology(XTP_Container_T & from_xtp_cont, CSG_Topology_T & to_csg_top){
        tools::StructureParameters cont_params = from_xtp_cont.getParameters();

        cont_params.set(tools::StructureParameter::MoleculeId,cont_params.get<int>(tools::StructureParameter::AtomContainerId));
        cont_params.set(tools::StructureParameter::MoleculeType,cont_params.get<std::string>(tools::StructureParameter::AtomContainerType));
        typename CSG_Topology_T::container_t mol = to_csg_top.CreateMolecule(cont_params);
        for( const typename XTP_Container_T::bead_t & bead : from_xtp_cont){
          tools::StructureParameters params = bead.getParameters();
          Eigen::VectorXd pos = params.get<Eigen::VectorXd>(tools::StructureParameter::XTP_Position);
          pos*= converter_.convert(XTP_Container_T::distance_unit,CSG_Topology_T::distance_unit);
          params.set(tools::StructureParameter::CSG_Position,pos);
          typename CSG_Topology_T::bead_t csg_bead = to_csg_top.CreateBead(params);
          mol.AddBead(&csg_bead);
        }
      }

    /*****************************************************************************
     * Converting from CSG to XTP 
     *****************************************************************************/
    template<class CSG_Topology_T, class Container_T>  
      inline void FileAdapter::convertCSGTopologyToXTPContainer(CSG_Topology_T & from_csg_top, Container_T & to_cont){
        assert(from_csg_top.MoleculeCount()==1 && "It is only possible to convert a csg topology object to a atom container when there is a single molecule in the topology object, otherwise it is not clear which molecule should be converted"); 

        convertCSGTopologyToXTPContainer(from_csg_top,to_cont,0);
      }

    template<class CSG_Topology_T, class Container_T>  
      inline void FileAdapter::convertCSGTopologyToXTPContainer(CSG_Topology_T & from_csg_top, Container_T & to_cont,const int & csg_molecule_id){
        assert(from_csg_top.MoleculeExist(csg_molecule_id) && "Cannot convert csg molecule object with provided id to a atom container, no csg molecule with the specified id exists"); 

        typename CSG_Topology_T::container_t csg_cont = from_csg_top.getMolecule(csg_molecule_id);
        tools::StructureParameters params = csg_cont.getParameters();
        /// Set the molecule id and type equal to the atom container id and type
        int mol_id = params.get<int>(tools::StructureParameter::MoleculeId);
        std::string molecule_type = params.get<std::string>(tools::StructureParameter::MoleculeType);
        std::unordered_map<int,std::vector<tools::StructureParameters>> mol_id_and_bead_params;

        std::vector<int> bead_ids = csg_cont.getBeadIds();
        for( const int & bead_id : bead_ids ){
          typename CSG_Topology_T::bead_t * bead = csg_cont.getBead(bead_id); 
          tools::StructureParameters params_bead = bead->getParameters();
          mol_id_and_bead_params[mol_id].push_back(params_bead);
        }

        to_cont = Container_T(molecule_type,mol_id);

        for( tools::StructureParameters bead_params : mol_id_and_bead_params[mol_id]){ 
          Eigen::VectorXd pos = bead_params.get<Eigen::VectorXd>(tools::StructureParameter::CSG_Position);
          pos*= converter_.convert(CSG_Topology_T::distance_unit,Container_T::distance_unit);
          bead_params.set(tools::StructureParameter::XTP_Position,pos);
          typename Container_T::bead_t bead(bead_params);
          to_cont.push_back(bead);
        }

      }

    template<class CSG_Topology_T, class XTP_Topology_T>  
      inline void FileAdapter::convertCSGTopologyToXTPTopology(CSG_Topology_T & from_csg_top, XTP_Topology_T & to_xtp_top){

        Eigen::Matrix3d box = from_csg_top.getBox() * 
          converter_.convert(CSG_Topology_T::distance_unit,XTP_Topology_T::distance_unit); 
        to_xtp_top.setBox(box);

        to_xtp_top.setStep(from_csg_top.getStep());
        to_xtp_top.setTime(from_csg_top.getTime() * 
            converter_.convert(CSG_Topology_T::time_unit,XTP_Topology_T::time_unit));

        std::vector<std::pair<int,std::string>> id_and_types;
        std::unordered_map<int,std::vector<tools::StructureParameters>> mol_id_and_bead_params;
        /// We must collect the molecules and beads in order because the 
        /// order that they are created in xtp is important. 
        for( typename CSG_Topology_T::container_t & csg_cont : from_csg_top){
          tools::StructureParameters params = csg_cont.getParameters();
          /// Set the molecule id and type equal to the atom container id and type
          int mol_id = params.get<int>(tools::StructureParameter::MoleculeId);
          std::string molecule_type = params.get<std::string>(tools::StructureParameter::MoleculeType);
          id_and_types.push_back(std::pair<int,std::string>(mol_id,molecule_type));
 
          std::vector<int> bead_ids = csg_cont.getBeadIds();
          for( const int & bead_id : bead_ids ){
            typename CSG_Topology_T::bead_t * bead = csg_cont.getBead(bead_id);
            tools::StructureParameters params_bead = bead->getParameters();
            mol_id_and_bead_params[mol_id].push_back(params_bead);
          }
        }

        std::sort(id_and_types.begin(),id_and_types.end()); 
        assert(id_and_types.at(0).first == 0 && "The smallest id of a molecule in the csg topology object is not 0, this means when converting to xtp topology segment the ids of the segments will not match those of the csg molecules.");
        for( std::pair<int,std::string> & id_and_type : id_and_types){
          to_xtp_top.AddSegment(id_and_type.second);

          typename XTP_Topology_T::container_t seg = to_xtp_top.getSegment(id_and_type.first);
          for( tools::StructureParameters bead_params : mol_id_and_bead_params[id_and_type.first]){ 
          Eigen::VectorXd pos = bead_params.get<Eigen::VectorXd>(tools::StructureParameter::CSG_Position);
          pos*= converter_.convert(CSG_Topology_T::distance_unit,XTP_Topology_T::distance_unit);
          bead_params.set(tools::StructureParameter::XTP_Position,pos);
            typename XTP_Topology_T::bead_t bead(bead_params);
            seg.push_back(bead);
          }
        } 
      }

  }
}

#endif // VOTCA_XTP_FILEADAPTER_H
