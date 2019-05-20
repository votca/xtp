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
#include <votca/csg/topology.h>
#include <Eigen/Dense>
#include "topology.h"

namespace votca {
  namespace xtp {


    /**
     * @brief Topology Converter
     *
     * This class is designed to interface with the topology structures in xtp
     * and in csg such that the xtp classes can use the csg topolology and 
     * trajectory input and output readers and writers
     */
    class TopologyConverter {

      public: 

        ::votca::csg::Topology Convert(const Topology & from_xtp_top) const;

        Topology Convert(const ::votca::csg::Topology & from_csg_top) const;

      private:

        tools::UnitConverter converter_;

    };

    /**
     * @brief Designed to convert between a xtp atom container or its children
     * and a csg::Topology object
     */
    class TopologyContainerConverter {

      public: 
        template<class XTP_Container_T>  
          XTP_Container_T Convert(const ::votca::csg::Topology & from_csg_top, const int & csg_molecule_id) const;

        template<class XTP_Container_T>
          ::votca::csg::Topology Convert(const XTP_Container_T & from_xtp_cont) const;

        template<class XTP_Container_T>  
          XTP_Container_T Convert(const ::votca::csg::Topology & from_csg_top) const;


      private:

        tools::UnitConverter converter_;

    };
    /*****************************************************************************
     * Converting from XTP to CSG 
     *****************************************************************************/

      inline ::votca::csg::Topology TopologyConverter::Convert(const Topology & from_xtp_top) const{

        ::votca::csg::Topology to_csg_top;
        /// Convert the xtp topology box boundaries to the csg boundaries
        /// by editing the units
        Eigen::Matrix3d box = from_xtp_top.getBox() * 
          converter_.convert(Topology::distance_unit,csg::Topology::distance_unit); 
        to_csg_top.setBox(box);

        to_csg_top.setStep(from_xtp_top.getStep());
        to_csg_top.setTime(from_xtp_top.getTime() * 
            converter_.convert(Topology::time_unit,csg::Topology::time_unit));

        for( const typename Topology::container_t & xtp_cont : from_xtp_top){
          tools::StructureParameters params = xtp_cont.getParameters();
          params.set(tools::StructureParameter::MoleculeId,params.get<int>(tools::StructureParameter::AtomContainerId));
          params.set(tools::StructureParameter::MoleculeType,params.get<std::string>(tools::StructureParameter::AtomContainerType));
          /// Set the molecule id and type equal to the atom container id and type
          typename ::votca::csg::Topology::container_t & mol = to_csg_top.CreateMolecule(params);  
          for( const typename Topology::bead_t & bead : xtp_cont ){

            tools::StructureParameters params_bead = bead.getParameters();
            Eigen::Vector3d pos = params_bead.get<Eigen::Vector3d>(tools::StructureParameter::XTP_Position);
            pos*= converter_.convert(Topology::distance_unit,csg::Topology::distance_unit);
            params_bead.set(tools::StructureParameter::CSG_Position,pos);
            typename ::votca::csg::Topology::bead_t & csg_bead = to_csg_top.CreateBead(params);
            mol.AddBead(csg_bead);
          }
        }

        return to_csg_top;
      }

      inline Topology TopologyConverter::Convert(const ::votca::csg::Topology & from_csg_top) const {

        Topology to_xtp_top;
        Eigen::Matrix3d box = from_csg_top.getBox() * 
          converter_.convert(csg::Topology::distance_unit,Topology::distance_unit); 
        to_xtp_top.setBox(box);

        to_xtp_top.setStep(from_csg_top.getStep());
        to_xtp_top.setTime(from_csg_top.getTime() * 
            converter_.convert(csg::Topology::time_unit,Topology::time_unit));

        std::vector<std::pair<int,std::string>> id_and_types;
        std::unordered_map<int,std::vector<tools::StructureParameters>> mol_id_and_bead_params;
        /// We must collect the molecules and beads in order because the 
        /// order that they are created in xtp is important. 
        for( const typename ::votca::csg::Topology::container_t & csg_cont : from_csg_top){
          tools::StructureParameters params = csg_cont.getParameters();
          /// Set the molecule id and type equal to the atom container id and type
          int mol_id = params.get<int>(tools::StructureParameter::MoleculeId);
          std::string molecule_type = params.get<std::string>(tools::StructureParameter::MoleculeType);
          id_and_types.push_back(std::pair<int,std::string>(mol_id,molecule_type));
 
          std::vector<int> bead_ids = csg_cont.getBeadIds();
          for( const int & bead_id : bead_ids ){
            const typename ::votca::csg::Topology::bead_t & bead = csg_cont.getBead(bead_id);
            tools::StructureParameters params_bead = bead.getParameters();
            mol_id_and_bead_params[mol_id].push_back(params_bead);
          }
        }

        std::sort(id_and_types.begin(),id_and_types.end()); 
        assert(id_and_types.at(0).first == 0 && "The smallest id of a molecule in the csg topology object is not 0, this means when converting to xtp topology segment the ids of the segments will not match those of the csg molecules.");
        for( std::pair<int,std::string> & id_and_type : id_and_types){
          to_xtp_top.AddSegment(id_and_type.second);

          typename Topology::container_t seg = to_xtp_top.getSegment(id_and_type.first);
          for( tools::StructureParameters bead_params : mol_id_and_bead_params[id_and_type.first]){ 
          Eigen::Vector3d pos = bead_params.get<Eigen::Vector3d>(tools::StructureParameter::CSG_Position);
          pos*= converter_.convert(csg::Topology::distance_unit,Topology::distance_unit);
          bead_params.set(tools::StructureParameter::XTP_Position,pos);
            typename Topology::bead_t bead(bead_params);
            seg.push_back(bead);
          }
        }
        return to_xtp_top;
      }

      /*****************************************************************************
     * Converting from CSG to XTP 
     *****************************************************************************/
  
      template<class XTP_Container_T>
      inline ::votca::csg::Topology TopologyContainerConverter::Convert(const XTP_Container_T & from_xtp_cont) const{

        ::votca::csg::Topology to_csg_top;
        tools::StructureParameters cont_params = from_xtp_cont.getParameters();

        cont_params.set(tools::StructureParameter::MoleculeId,cont_params.get<int>(tools::StructureParameter::AtomContainerId));
        cont_params.set(tools::StructureParameter::MoleculeType,cont_params.get<std::string>(tools::StructureParameter::AtomContainerType));
        typename ::votca::csg::Topology::container_t & mol = to_csg_top.CreateMolecule(cont_params);
        for( const typename XTP_Container_T::bead_t & bead : from_xtp_cont){
          tools::StructureParameters params = bead.getParameters();
          tools::byte_t symmetry = 1;                                               
          double mass = 0.0;                                                        
          double charge = 0.0;                                                      
          // Because xyz files do not contain any bond info each atom is assumed       
          // to be in its own molecule                                              
          int residue_id = tools::topology_constants::unassigned_residue_id;                         
          std::string residue_type = tools::topology_constants::unassigned_residue_type;             

          params.set(tools::StructureParameter::Symmetry, symmetry);                
          params.set(tools::StructureParameter::CSG_Mass, mass);                    
          params.set(tools::StructureParameter::CSG_Charge, charge);                
          params.set(tools::StructureParameter::MoleculeId, mol.getId());           
          params.set(tools::StructureParameter::MoleculeType, mol.getType());          
          params.set(tools::StructureParameter::ResidueId, residue_id);             
          params.set(tools::StructureParameter::ResidueType, residue_type);         
          Eigen::Vector3d pos = params.get<Eigen::Vector3d>(tools::StructureParameter::XTP_Position);
          pos*= converter_.convert(XTP_Container_T::distance_unit,csg::Topology::distance_unit);
          params.set(tools::StructureParameter::CSG_Position,pos);
          csg::Topology::bead_t & csg_bead = to_csg_top.CreateBead(params);
          mol.AddBead(csg_bead);

        }
        return to_csg_top;
      }

      template<class Container_T>  
      inline Container_T TopologyContainerConverter::Convert(const ::votca::csg::Topology & from_csg_top,const int & csg_molecule_id) const {
        assert(from_csg_top.MoleculeExist(csg_molecule_id) && "Cannot convert csg molecule object with provided id to a atom container, no csg molecule with the specified id exists"); 

        typename ::votca::csg::Topology::container_t csg_cont = from_csg_top.getMolecule(csg_molecule_id);
        tools::StructureParameters params = csg_cont.getParameters();
        /// Set the molecule id and type equal to the atom container id and type
        int mol_id = params.get<int>(tools::StructureParameter::MoleculeId);
        std::unordered_map<int,std::vector<tools::StructureParameters>> mol_id_and_bead_params;

        std::vector<int> bead_ids = csg_cont.getBeadIds();
        for( const int & bead_id : bead_ids ){
          typename ::votca::csg::Topology::bead_t & bead = csg_cont.getBead(bead_id); 
          tools::StructureParameters params_bead = bead.getParameters();
          mol_id_and_bead_params[mol_id].push_back(params_bead);
        }

        Container_T to_cont(params);

        for( tools::StructureParameters bead_params : mol_id_and_bead_params[mol_id]){ 
          Eigen::Vector3d pos = bead_params.get<Eigen::Vector3d>(tools::StructureParameter::CSG_Position);
          pos*= converter_.convert(csg::Topology::distance_unit,Container_T::distance_unit);
          bead_params.set(tools::StructureParameter::XTP_Position,pos);
          typename Container_T::bead_t bead(bead_params);
          to_cont.push_back(bead);
        }

        return to_cont;
      }

      template<class Container_T>  
      inline Container_T TopologyContainerConverter::Convert(const ::votca::csg::Topology & from_csg_top) const {

        /// Here all the atoms are places in a single container, no assumptions
        /// are made as to whether the atoms/beads belong to a single molecule
        /// or not.

        std::vector<std::pair<int,tools::StructureParameters>> bead_id_and_params;
        for( const typename ::votca::csg::Topology::container_t & csg_cont : from_csg_top){

          std::vector<int> bead_ids = csg_cont.getBeadIds();
          for( const int & bead_id : bead_ids ){
            const typename ::votca::csg::Topology::bead_t & bead = csg_cont.getBead(bead_id);
            tools::StructureParameters params_bead = bead.getParameters();
            bead_id_and_params.push_back(
                std::pair<int,tools::StructureParameters>(bead_id,params_bead));
          }
        }

        tools::StructureParameters cont_params;
        cont_params.set(tools::StructureParameter::AtomContainerType,
            tools::topology_constants::unassigned_atom_container_type);
        cont_params.set(tools::StructureParameter::AtomContainerId,
            tools::topology_constants::unassigned_atom_container_id);

        Container_T to_cont(cont_params);
        std::sort(bead_id_and_params.begin(),bead_id_and_params.end(),
            [](const std::pair<int,tools::StructureParameters> & left,
               const std::pair<int,tools::StructureParameters> & right)
            { return left.first < right.first; }); 

        assert(bead_id_and_params.at(0).first == 0 && "The smallest id of a bead in the csg topology object is not 0, this means when converting to xtp topology atom container the ids of the segments will not match those of the csg molecules.");
        for( std::pair<int,tools::StructureParameters> & id_and_params : bead_id_and_params){
          Eigen::Vector3d pos = id_and_params.second.get<Eigen::Vector3d>(tools::StructureParameter::CSG_Position);
          pos*= converter_.convert(csg::Topology::distance_unit,Topology::distance_unit);
          id_and_params.second.set(tools::StructureParameter::XTP_Position,pos);
          std::cout << "Dist " << id_and_params.second.get<Eigen::Vector3d>(tools::StructureParameter::XTP_Position) << std::endl;
          std::cout << "Elem " << id_and_params.second.get<std::string>(tools::StructureParameter::Element) << std::endl;
          std::cout << "type " << id_and_params.second.get<std::string>(tools::StructureParameter::BeadType) << std::endl;
          to_cont.push_back(typename Container_T::bead_t(id_and_params.second));
        }
        return to_cont;
      }

  }
}

#endif // VOTCA_XTP_FILEADAPTER_H
