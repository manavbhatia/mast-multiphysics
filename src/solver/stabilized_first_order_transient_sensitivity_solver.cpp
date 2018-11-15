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


// MAST includes
#include "solver/stabilized_first_order_transient_sensitivity_solver.h"
#include "base/transient_assembly_elem_operations.h"
#include "base/elem_base.h"
#include "base/nonlinear_system.h"
#include "base/assembly_base.h"
#include "base/system_initialization.h"
#include "base/output_assembly_elem_operations.h"

// libMesh includes
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/linear_solver.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/enum_norm_type.h"


MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
StabilizedFirstOrderNewmarkTransientSensitivitySolver():
MAST::TransientSolverBase(1, 2),
max_amp         (1.),
beta            (1.),
max_index       (0),
_assemble_mass  (false),
_t0             (0.),
_index0         (0),
_index1         (0) {
    
}


MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
~StabilizedFirstOrderNewmarkTransientSensitivitySolver() {
    
}


void
MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
sensitivity_solve(MAST::AssemblyBase& assembly,
                  const MAST::FunctionBase& f) {
    
    libmesh_assert_greater(  max_amp, 0.);
    libmesh_assert_greater(max_index,  0);
    
    
    // Log how long the linear solve takes.
    LOG_SCOPE("sensitivity_solve()", "StabilizedSensitivity");
    
    MAST::NonlinearSystem& sys = assembly.system();
    
    assembly.set_elem_operation_object(*this);
    
    libMesh::SparseMatrix<Real>
    &J_int = *sys.matrix_A,
    &M_avg = *sys.matrix_B;
    
    libMesh::NumericVector<Real>
    &sol    = this->solution(),
    &dsol0  = this->solution_sensitivity(1),    // previous solution
    &dsol1  = this->solution_sensitivity(),     // current solution
    &rhs    = sys.add_sensitivity_rhs();
    
    std::unique_ptr<libMesh::NumericVector<Real>>
    f_int(rhs.zero_clone().release()),
    vec1(rhs.zero_clone().release());
    
    J_int.zero();
    M_avg.zero();
    rhs.zero();
    J_int.close();
    M_avg.close();
    rhs.close();
    
    _index0     = _index1;
    
    Real
    amp         = 0.;
    
    _t0         = sys.time;
    sys.time   += this->dt;
    _index1++;
    
    bool
    continue_it = true;
    
    // now iterate till a suitable time-step has been found
    while (continue_it) {
        
        // update the solution vector
        std::ostringstream oss;
        oss << "output_sol_t_" << _index1;
        sys.read_in_vector(sol, "data_M0p25_forced", oss.str(), true);

        // assemble the Jacobian matrix
        _assemble_mass = false;
        assembly.residual_and_jacobian(sol, nullptr, sys.matrix, sys);
        J_int.add(this->dt, *sys.matrix);
        
        // assemble the mass matrix
        _assemble_mass = true;
        assembly.residual_and_jacobian(sol, nullptr, sys.matrix, sys);
        Mat M_avg_mat = dynamic_cast<libMesh::PetscMatrix<Real>&>(M_avg).mat();
        PetscErrorCode ierr = MatScale(M_avg_mat, (sys.time - _t0 - this->dt)/(sys.time - _t0));
        CHKERRABORT(sys.comm().get(), ierr);
        M_avg.add(this->dt/(sys.time - _t0), *sys.matrix);
        
        // assemble the forcing vector
        assembly.sensitivity_assemble(f, rhs);
        rhs.scale(-1.);
        f_int->add(this->dt, rhs);  // int_0^t F dt
        f_int->close();
        rhs.zero();
        M_avg.vector_mult(rhs, dsol0); // M_avg x0
        rhs.add(1., *f_int);          // M_avg x0 + int_0^t F dt
        rhs.close();
        
        // close the quantities
        J_int.close();
        M_avg.close();
        
        // now compute the matrix for the linear solution
        // Note that the stabilized solver for system definde as
        // M x_dot = f(x) will define the Jacobian as  (M_avg - Jint).
        // However, since MAST defines the system as  M x_dot - f(x) = 0
        // and returns the Jacobian as J = -df(x)/dx, we will not multiply
        // the Jacobian with -1.
        sys.matrix->zero();
        sys.matrix->add(1., M_avg);
        sys.matrix->add(1., J_int); // see explanation above for positive sign.
        sys.matrix->close();
        
        // The sensitivity problem is linear
        // Our iteration counts and residuals will be sums of the individual
        // results
        std::pair<unsigned int, Real>
        solver_params = sys.get_linear_solve_parameters();
        
        // Solve the linear system.
        libMesh::SparseMatrix<Real> * pc = sys.request_matrix("Preconditioner");
        
        std::pair<unsigned int, Real> rval =
        sys.linear_solver->solve (*sys.matrix, pc,
                                  dsol1,
                                  rhs,
                                  solver_params.second,
                                  solver_params.first);
        
        // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
        sys.get_dof_map().enforce_constraints_exactly (sys, &dsol1, /* homogeneous = */ true);
#endif

        M_avg.vector_mult(*vec1, dsol1);

        // estimate the amplification factor
        amp = _compute_amplification_factor(rhs, *vec1);
        
        if (amp <= this->max_amp ||  _index1 >= max_index) {
            
            continue_it = false;
            
            // we may have reached here due to end of time integration.
            // In this case, the final Dmag may be unstable. If this is the
            // case, we should throw that solution out.
            if (amp <= this->max_amp) {
                
                libMesh::out << "accepting sol" << std::endl;
                dsol0.zero();
                dsol0.add(1., dsol1);
                dsol0.close();
            }
            else {
                // present sol is same as previous sol
                dsol1.zero();
                dsol1.add(1., dsol0);
                dsol1.zero();
                
                libMesh::out << "reached final step. Rjecting solution" << std::endl;
            }
        }
        else {
            
            _index1++;
            sys.time    += this->dt;
        }
    }
    
    
    assembly.clear_elem_operation_object();
}


Real
MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
evaluate_q_sens_for_previous_interval(MAST::AssemblyBase& assembly,
                                      const MAST::FunctionBase& p,
                                      MAST::OutputAssemblyElemOperations& output) {
    
    
    MAST::NonlinearSystem& sys = assembly.system();
    
    Real
    eta         = 0.,
    t1          = sys.time,
    q           = 0.;
    
    libMesh::NumericVector<Real>
    &sol         = this->solution(),
    &dx0         = this->solution_sensitivity(1),
    &dx1         = this->solution_sensitivity();
    
    std::unique_ptr<libMesh::NumericVector<Real>>
    dx(dx0.zero_clone().release());
    
    libmesh_assert_greater(_index1, _index0);
    
    for (unsigned int i=_index0; i<_index1; i++) {
        
        // update the nonlinear solution
        std::ostringstream oss;
        oss << "output_sol_t_" << _index1;
        sys.read_in_vector(sol, "data_M0p25_forced", oss.str(), true);
        
        sys.time  = _t0 + this->dt * ( i - _index0);
        
        eta = (sys.time-_t0)/(t1-_t0);
        
        dx->zero();
        dx->add(1.-eta, dx0);
        dx->add(   eta, dx1);
        dx->close();
        
        output.zero_for_analysis();
        assembly.calculate_output_direct_sensitivity(sol, *dx, p, output);
        
        q  += output.output_sensitivity_total(p) * this->dt;
    }
    
    return q;
}



void
MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
set_element_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                 const std::vector<libMesh::NumericVector<Real>*>& sols) {
    
    libmesh_assert_equal_to(sols.size(), 2);
    
    const unsigned int n_dofs = (unsigned int)dof_indices.size();
    
    // get the current state and velocity estimates
    // also get the current discrete velocity replacement
    RealVectorX
    sol          = RealVectorX::Zero(n_dofs),
    vel          = RealVectorX::Zero(n_dofs);
    
    
    const libMesh::NumericVector<Real>
    &sol_global = *sols[0],
    &vel_global = *sols[1];
    
    // get the references to current and previous sol and velocity
    
    for (unsigned int i=0; i<n_dofs; i++) {
        
        sol(i)          = sol_global(dof_indices[i]);
        vel(i)          = vel_global(dof_indices[i]);
    }
    
    _assembly_ops->set_elem_solution(sol);
    _assembly_ops->set_elem_velocity(vel);
}



void
MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
extract_element_sensitivity_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                                 const std::vector<libMesh::NumericVector<Real>*>& sols,
                                 std::vector<RealVectorX>& local_sols) {
    
    libmesh_assert_equal_to(sols.size(), 2);
    
    const unsigned int n_dofs = (unsigned int)dof_indices.size();
    
    local_sols.resize(2);
    
    RealVectorX
    &sol         = local_sols[0],
    &vel         = local_sols[1];
    
    sol          = RealVectorX::Zero(n_dofs);
    vel          = RealVectorX::Zero(n_dofs);
    
    const libMesh::NumericVector<Real>
    &sol_global = *sols[0],
    &vel_global = *sols[1];
    
    // get the references to current and previous sol and velocity
    
    for (unsigned int i=0; i<n_dofs; i++) {
        
        sol(i)          = sol_global(dof_indices[i]);
        vel(i)          = vel_global(dof_indices[i]);
    }
    
}


void
MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
update_velocity(libMesh::NumericVector<Real>&       vec,
                const libMesh::NumericVector<Real>& sol) {
    
    const libMesh::NumericVector<Real>
    &prev_sol = this->solution(1),
    &prev_vel = this->velocity(1);
    
    vec.zero();
    vec.add( 1.,      sol);
    vec.add(-1., prev_sol);
    vec.scale(1./beta/dt);
    vec.close();
    vec.add(-(1.-beta)/beta, prev_vel);
    
    vec.close();
}


void
MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
update_sensitivity_velocity(libMesh::NumericVector<Real>&       vec,
                            const libMesh::NumericVector<Real>& sol) {
    
    const libMesh::NumericVector<Real>
    &prev_sol = this->solution_sensitivity(1),
    &prev_vel = this->velocity_sensitivity(1);
    
    vec.zero();
    vec.add( 1.,      sol);
    vec.add(-1., prev_sol);
    vec.scale(1./beta/dt);
    vec.close();
    vec.add(-(1.-beta)/beta, prev_vel);
    
    vec.close();
}



void
MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
elem_calculations(bool if_jac,
                  RealVectorX& vec,
                  RealMatrixX& mat) {
    // make sure that the assembly object is provided
    libmesh_assert(_assembly_ops);
    unsigned int n_dofs = (unsigned int)vec.size();
    
    RealVectorX
    f_x     = RealVectorX::Zero(n_dofs),
    f_m     = RealVectorX::Zero(n_dofs);
    
    RealMatrixX
    f_m_jac_xdot  = RealMatrixX::Zero(n_dofs, n_dofs),
    f_m_jac       = RealMatrixX::Zero(n_dofs, n_dofs),
    f_x_jac       = RealMatrixX::Zero(n_dofs, n_dofs);
    
    // perform the element assembly
    _assembly_ops->elem_calculations(if_jac,
                                     f_m,           // mass vector
                                     f_x,           // forcing vector
                                     f_m_jac_xdot,  // Jac of mass wrt x_dot
                                     f_m_jac,       // Jac of mass wrt x
                                     f_x_jac);      // Jac of forcing vector wrt x
    
    if (_assemble_mass)
        mat = f_m_jac_xdot;
    else
        mat = f_m_jac + f_x_jac;
}



void
MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
elem_sensitivity_calculations(const MAST::FunctionBase& f,
                              RealVectorX& vec) {
    
    // make sure that the assembly object is provided
    libmesh_assert(_assembly_ops);
    unsigned int n_dofs = (unsigned int)vec.size();
    
    RealVectorX
    f_x     = RealVectorX::Zero(n_dofs),
    f_m     = RealVectorX::Zero(n_dofs);
    
    // perform the element assembly
    _assembly_ops->elem_sensitivity_calculations(f,
                                                 f_m,           // mass vector
                                                 f_x);          // forcing vector
    
    // sensitivity of system residual involving mass and internal force term
    vec  = (f_m + f_x);
}




void
MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
elem_sensitivity_contribution_previous_timestep(const std::vector<RealVectorX>& prev_sols,
                                                RealVectorX& vec) {
    
    vec.setZero();
}


Real
MAST::StabilizedFirstOrderNewmarkTransientSensitivitySolver::
_compute_amplification_factor(const libMesh::NumericVector<Real>& sol0,
                              const libMesh::NumericVector<Real>& sol1) {
    
    libmesh_assert(_system);
    
    MAST::NonlinearSystem&
    sys = _system->system();
    
    Real
    norm0  = 0.,
    norm1  = 0.;
    
    unsigned int
    n_vars = _system->n_vars();
    
    RealVectorX
    var_norm    =  RealVectorX::Zero(n_vars);
    
    for (unsigned int i=0; i<n_vars; i++) {
        
        norm0 = sys.calculate_norm(sol0, i, libMesh::L2);
        norm1 = sys.calculate_norm(sol1, i, libMesh::L2);

        if (norm0 > 0.)
            var_norm(i) = norm1/norm0;
    }
    
    libMesh::out << "amp coeffs: " << var_norm.transpose() << std::endl;
    return var_norm.maxCoeff();
}
