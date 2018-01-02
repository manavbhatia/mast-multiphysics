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

#ifndef __mast_nonlinear_implicit_assembly_elem_operation_h__
#define __mast_nonlinear_implicit_assembly_elem_operation_h__

// MAST includes
#include "base/assembly_elem_operation.h"
#include "base/mast_data_types.h"


namespace MAST {
    
    class NonlinearImplicitAssemblyElemOperations:
    public MAST::AssemblyElemOperations {
        
    public:
        NonlinearImplicitAssemblyElemOperations();
        
        virtual ~NonlinearImplicitAssemblyElemOperations();
        
        /*!
         *   sets the element solution(s) before calculations
         */
        virtual void set_elem_sol(MAST::ElementBase& elem,
                                  const RealVectorX& sol);
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void elem_calculations(MAST::ElementBase& elem,
                                       bool if_jac,
                                       RealVectorX& vec,
                                       RealMatrixX& mat) = 0;
        
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector quantity in \par vec. The vector quantity only
         *   include the \f$ [J] \{dX\} f$ components, so the inherited classes
         *   must ensure that no component of constant forces (traction/body
         *   forces/etc.) are added to this vector.
         */
        virtual void
        elem_linearized_jacobian_solution_product(MAST::ElementBase& elem,
                                                  RealVectorX& vec) = 0;
        
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                   RealVectorX& vec) = 0;
        
        
        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \par elem,
         *   and returns the matrix in \par vec .
         */
        virtual void
        elem_second_derivative_dot_solution_assembly(MAST::ElementBase& elem,
                                                     RealMatrixX& mat) = 0;
        
        
        /*!
         *    a helper function to evaluate the numerical Jacobian
         *    and compare it with the analytical Jacobian.
         */
        void check_element_numerical_jacobian(MAST::ElementBase& e,
                                              RealVectorX& sol);

    protected:
        
    };
}


#endif // __mast_nonlinear_implicit_assembly_elem_operation_h__

