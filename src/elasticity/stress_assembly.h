/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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

#ifndef __mast__stress_assembly__
#define __mast__stress_assembly__

// MAST includes
#include "base/assembly_base.h"

// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"


namespace MAST {
    
    // Forward declerations
    class StressStrainOutputBase;
    
    
    class StressAssembly:
    public MAST::AssemblyBase {
    public:
        
        
        /*!
         *   constructor associates this assembly object with the system
         */
        StressAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~StressAssembly();
        
        
        /*!
         *   attaches a system to this discipline, and vice-a-versa
         */
        virtual void
        attach_discipline_and_system(MAST::StressStrainOutputBase& elem_ops,
                                     MAST::PhysicsDisciplineBase& discipline,
                                     MAST::SystemInitialization& system);
        
        
        /*!
         *   Reattaches to the same system that was attached earlier.
         *
         *   This cannot be called if the clear_discipline_and_system() method
         *   has been called.
         */
        virtual void
        reattach_to_system();
        
        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system();
        
        
        /*!
         *   updates the stresses and strains for the specified solution
         *   vector \p X. Only the maximum values out of each element are
         *   updated. This will put the stress data in the System::solution
         *   vector related to stress/strain values.
         */
        void update_stress_strain_data(const libMesh::NumericVector<Real>& X) const;

        /*!
         *   updates the sensitivity of stresses and strains for the specified
         *   solution vector \p X and its sensitivity, \p dXdp, with respect
         *   to parameter \p p. Only the maximum values out of each element are
         *   updated.
         */
        void update_stress_strain_sensitivity_data(const libMesh::NumericVector<Real>& X,
                                                   const libMesh::NumericVector<Real>& dXdp,
                                                   const MAST::FunctionBase& p,
                                                   libMesh::NumericVector<Real>& dsigmadp) const;


    protected:
        

        MAST::StressStrainOutputBase* _stress_ops;
        
    };
}


#endif //__mast__stress_assembly__
