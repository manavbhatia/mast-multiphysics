/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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


#ifndef __mast_test_build_structural_element_2d_h__
#define __mast_test_build_structural_element_2d_h__

// C++ includes
#include <memory>
#include <vector>

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
    class Solid2DSectionElementPropertyCard;
    class PistonTheoryBoundaryCondition;
    
    
    struct BuildStructural2DElem {
        
        
        BuildStructural2DElem();
        
        
        ~BuildStructural2DElem();
        
        
        // create the mesh
        libMesh::SerialMesh*           _mesh;
        
        // create the equation system
        libMesh::EquationSystems*      _eq_sys;
        
        // create the libmesh system
        libMesh::NonlinearImplicitSystem*  _sys;
        
        // initialize the system to the right set of variables
        MAST::StructuralSystemInitialization* _structural_sys;
        MAST::StructuralDiscipline*           _discipline;
        
        // create the property functions and add them to the
        MAST::Parameter
        *_thz,
        *_E,
        *_nu,
        *_kappa,
        *_zero;
        
        MAST::ConstantFieldFunction
        *_thz_f,
        *_E_f,
        *_nu_f,
        *_kappa_f,
        *_hzoff_f;
        
        // create the material property card
        MAST::IsotropicMaterialPropertyCard*     _m_card;
        
        // create the element property card
        MAST::Solid2DSectionElementPropertyCard* _p_card;
        
        // create the boundary condition
        MAST::PistonTheoryBoundaryCondition*    _p_theory;
        
        // vector of parameters to evaluate sensitivity wrt
        std::vector<MAST::Parameter*> _params_for_sensitivity;
    };
}


#endif // __mast_test_build_structural_element_2d_h__