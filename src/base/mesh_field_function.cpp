
// MAST includes
#include "base/mesh_field_function.h"
#include "base/system_initialization.h"

// libMesh includes
#include "libmesh/dof_map.h"


template <typename ValType>
MAST::MeshFieldFunction<ValType>::
MeshFieldFunction(const std::string& nm):
MAST::FieldFunction<ValType>(nm),
_use_qp_sol(false),
_qp_sol(),
_system(NULL),
_sol(NULL),
_mesh_function(NULL)
{ }



template <typename ValType>
MAST::MeshFieldFunction<ValType>::
MeshFieldFunction(const MAST::MeshFieldFunction<ValType>& f):
MAST::FieldFunction<ValType>(f),
_system(NULL),
_sol(NULL),
_mesh_function(NULL) {
    
    // a cloned function still works with the master function's
    // solution and mesh funciton data structures
    // this is not thread-safe.
    
}




template <typename ValType>
MAST::MeshFieldFunction<ValType>::~MeshFieldFunction() {
 
    this->clear();
}




template <typename ValType>
std::auto_ptr<MAST::FieldFunction<ValType> >
MAST::MeshFieldFunction<ValType>::clone() const {
    
    MAST::MeshFieldFunction<ValType>* rval =
    new MAST::MeshFieldFunction<ValType>(*this);
    
    return std::auto_ptr<MAST::FieldFunction<ValType> >(rval);
}




template <>
void
MAST::MeshFieldFunction<RealVectorX>::operator() (const libMesh::Point& p,
                                                  const Real t,
                                                  RealVectorX& v) const {
    
    // get the pointer to the mesh function from the master
    const MAST::MeshFieldFunction<RealVectorX>* master =
    dynamic_cast<const MAST::MeshFieldFunction<RealVectorX>*>(this->master());
    
    // if the element has provided a quadrature point solution,
    // then use it
    if (master->_use_qp_sol) {
        v = master->_qp_sol;
        return;
    }
    
    // make sure that the object was initialized
    libmesh_assert(master->_mesh_function);
    
    DenseRealVector v1;
    (*master->_mesh_function)(p, t, v1);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert(v1.size());
    
    // now copy this to the output vector
    v = RealVectorX::Zero(v1.size());
    for (unsigned int i=0; i<v1.size(); i++)
        v(i) = v1(i);
}




template <>
void
MAST::MeshFieldFunction<RealVectorX>::derivative (const MAST::DerivativeType d,
                                                  const MAST::FunctionBase& f,
                                                  const libMesh::Point& p,
                                                  const Real t,
                                                  RealVectorX& v) const {
    
    
    libmesh_error_msg("To be implemented.");
}






template <typename ValType>
void
MAST::MeshFieldFunction<ValType>::
init_for_system_and_solution (MAST::SystemInitialization& sys,
                              const libMesh::NumericVector<Real>& sol) {
    
    // this is to be done for the master only
    const MAST::MeshFieldFunction<ValType>* master_const =
    dynamic_cast<const MAST::MeshFieldFunction<ValType>*>(this->master());
    MAST::MeshFieldFunction<ValType>* master =
    const_cast<MAST::MeshFieldFunction<ValType>*>(master_const);
    
    
    // first make sure that the object is not already initialized
    libmesh_assert(!master->_mesh_function);
    
    // attach the system
    master->_system = &sys;
    
    libMesh::System& system = sys.system();
    
    // next, clone this solution and localize to the sendlist
    master->_sol = libMesh::NumericVector<Real>::build(system.comm()).release();

    const std::vector<libMesh::dof_id_type>& send_list =
    system.get_dof_map().get_send_list();
    
    // initialize and then localize the vector with the provided solution
    master->_sol->init(system.n_dofs(),
                       system.n_local_dofs(),
                       send_list,
                       false,
                       libMesh::GHOSTED);
    sol.localize(*master->_sol, send_list);
    
    // finally, create the mesh interpolation function
    master->_mesh_function = new libMesh::MeshFunction(system.get_equation_systems(),
                                                       *master->_sol,
                                                       system.get_dof_map(),
                                                       sys.vars());
    master->_mesh_function->init();
}




template <typename ValType>
void
MAST::MeshFieldFunction<ValType>::clear() {
    
    //only the master function will call this on itself
    if (this->master() != this)
        return;
    
    // if a pointer has been attached, then delete it and the
    // associated vector, and clear the associated system
    if (_mesh_function) {
        delete _mesh_function;
        _mesh_function = NULL;
        
        delete _sol;
        _sol = NULL;
        
        _system = NULL;
    }
    
    // clear flags for quadrature point solution
    _use_qp_sol = false;
}




template <typename ValType>
void
MAST::MeshFieldFunction<ValType>::
set_element_quadrature_point_solution(RealVectorX& sol) {
    
    // this is to be done for the master only
    const MAST::MeshFieldFunction<ValType>* master_const =
    dynamic_cast<const MAST::MeshFieldFunction<ValType>*>(this->master());
    MAST::MeshFieldFunction<ValType>* master =
    const_cast<MAST::MeshFieldFunction<ValType>*>(master_const);

    master->_use_qp_sol = true;
    master->_qp_sol     = sol;
}



template <typename ValType>
void
MAST::MeshFieldFunction<ValType>::
clear_element_quadrature_point_solution() {

    // this is to be done for the master only
    const MAST::MeshFieldFunction<ValType>* master_const =
    dynamic_cast<const MAST::MeshFieldFunction<ValType>*>(this->master());
    MAST::MeshFieldFunction<ValType>* master =
    const_cast<MAST::MeshFieldFunction<ValType>*>(master_const);
    
    master->_use_qp_sol = false;
    master->_qp_sol.setZero();
}




// explicit instantiations
template class MAST::MeshFieldFunction<RealVectorX>;

