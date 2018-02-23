/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef __mast_beam_thermally_stressed_modal_analysis_h__
#define __mast_beam_thermally_stressed_modal_analysis_h__


// C++ includes
#include <memory>
#include <vector>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/fe_type.h"
#include "libmesh/dof_map.h"



namespace MAST {
    
    // Forward declerations
    class StructuralSystemInitialization;
    class StructuralDiscipline;
    class Parameter;
    class ConstantFieldFunction;
    class IsotropicMaterialPropertyCard;
    class Solid1DSectionElementPropertyCard;
    class DirichletBoundaryCondition;
    class BoundaryConditionBase;
    class StressStrainOutputBase;
    class NonlinearSystem;
    
    
    struct BeamThermallyStressedModalAnalysis {
        
        
        BeamThermallyStressedModalAnalysis();
        
        
        ~BeamThermallyStressedModalAnalysis();
        
        
        /*!
         *   initializes the object for specified characteristics
         */
        void init(libMesh::ElemType etype, bool if_nonlin);

        
        /*!
         *   @returns a pointer to the parameter of the specified name.
         *   If no parameter exists by the specified name, then a \p nullptr
         *   pointer is returned and a message is printed with a valid list
         *   of parameters.
         */
        MAST::Parameter* get_parameter(const std::string& nm);
        
        /*!
         *  solves the system and returns the final solution
         */
        const libMesh::NumericVector<Real>&
        solve(bool if_write_output = false,
              std::vector<Real>* eig = nullptr);
        
        
        /*!
         *  solves the sensitivity of system and returns the final solution
         */
        const libMesh::NumericVector<Real>&
        sensitivity_solve(MAST::Parameter& p,
                          std::vector<Real>& eig);
        
        
        /*!
         *   clears the stress data structures for a followup analysis
         */
        void clear_stresss();
        
        bool _initialized;
        
        // length of domain
        Real _length;

        // create the mesh
        libMesh::SerialMesh*           _mesh;
        
        // create the equation system
        libMesh::EquationSystems*      _eq_sys;
        
        // create the libmesh system
        MAST::NonlinearSystem*  _sys;
        
        // initialize the system to the right set of variables
        MAST::StructuralSystemInitialization* _structural_sys;
        MAST::StructuralDiscipline*           _discipline;
        
        // create the property functions and add them to the
        MAST::Parameter
        *_thy,
        *_thz,
        *_alpha,
        *_rho,
        *_E,
        *_nu,
        *_temp,
        *_press,
        *_zero;
        
        MAST::ConstantFieldFunction
        *_thy_f,
        *_thz_f,
        *_alpha_f,
        *_rho_f,
        *_E_f,
        *_nu_f,
        *_hyoff_f,
        *_hzoff_f,
        *_temp_f,
        *_ref_temp_f,
        *_press_f;
        
        // create the material property card
        MAST::IsotropicMaterialPropertyCard*     _m_card;
        
        // create the element property card
        MAST::Solid1DSectionElementPropertyCard* _p_card;
        
        // create the Dirichlet boundary condition on left edge
        MAST::DirichletBoundaryCondition*     _dirichlet_left;
        
        // create the Dirichlet boundary condition on right edge
        MAST::DirichletBoundaryCondition*     _dirichlet_right;
        
        // create the pressure boundary condition
        MAST::BoundaryConditionBase*             _T_load;
        MAST::BoundaryConditionBase*             _p_load;
        
        // vector of parameters to evaluate sensitivity wrt
        std::vector<MAST::Parameter*> _params_for_sensitivity;
        
        // output quantity objects to evaluate stress
        std::vector<MAST::StressStrainOutputBase*>  _outputs;
    };
}



#endif //  __mast_beam_thermally_stressed_modal_analysis_h__