
/*
 * the objective of this code is to design a stiffened panel, thermally stressed that would avoid static instability.
 * First bifurcation is a static instability and we want to design the panel such that it does not get to it.
 * The thermal load is applied uniformly accross the panel.
 * Von Karman strain is used for the calculation of strain.
 */

// C++ includes
#include <memory>
#include <iomanip>
#include <map>
#include <iostream>
#include <math.h>

// MAST includes
#include "examples/old/structural/stiffened_plate_optimization/stiffened_plate_optimization_base.h"
#include "examples/base/multilinear_interpolation.h"
#include "examples/old/structural/stiffened_plate_optimization/stiffened_plate_optimization_base.h"
#include "examples/structural/base/thermal_stress_jacobian_scaling_function.h"
#include "examples/base/input_wrapper.h"
#include "examples/fluid/meshing/cylinder.h"
#include "examples/base/multilinear_interpolation.h"
#include "examples/structural/base/blade_stiffened_panel_mesh.h"
#include "examples/base/plot_results.h"


#include "base/field_function_base.h"
#include "base/physics_discipline_base.h"
#include "level_set/level_set_discipline.h"
#include "level_set/level_set_system_initialization.h"
#include "level_set/level_set_eigenproblem_assembly.h"
#include "level_set/level_set_transient_assembly.h"
#include "level_set/level_set_nonlinear_implicit_assembly.h"
#include "level_set/level_set_reinitialization_transient_assembly.h"
#include "level_set/level_set_volume_output.h"
#include "level_set/level_set_perimeter_output.h"
#include "level_set/level_set_boundary_velocity.h"
#include "level_set/indicator_function_constrain_dofs.h"
#include "level_set/level_set_constrain_dofs.h"
#include "level_set/level_set_intersection.h"
#include "level_set/filter_base.h"
#include "level_set/level_set_parameter.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/stress_output_base.h"
#include "elasticity/level_set_stress_assembly.h"
#include "elasticity/structural_system_initialization.h"
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "base/constant_field_function.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "solver/slepc_eigen_solver.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "optimization/gcmma_optimization_interface.h"
#include "optimization/npsol_optimization_interface.h"

#include "base/nonlinear_system.h"
#include "elasticity/stress_output_base.h"
#include "optimization/optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/structural_fluid_interaction_assembly.h"
#include "elasticity/structural_near_null_vector_space.h"
#include "aeroelasticity/time_domain_flutter_solver.h"
#include "aeroelasticity/time_domain_flutter_root.h"
#include "solver/slepc_eigen_solver.h"
#include "solver/pseudo_arclength_continuation_solver.h"

//#include "src/elasticity/structural_discipline.h"
#include "elasticity/structural_system_initialization.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "optimization/gcmma_optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "boundary_condition/dirichlet_boundary_condition.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/fe_type.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh_function.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/fe_type.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/dof_map.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/nonlinear_solver.h"

#include "libmesh/parallel.h"

void
stiffened_plate_thermally_stressed_piston_theory_flutter_optim_obj(int*    mode,
                                                                   int*    n,
                                                                   double* x,
                                                                   double* f,
                                                                   double* g,
                                                                   int*    nstate);

void
stiffened_plate_thermally_stressed_piston_theory_flutter_optim_con(int*    mode,
                                                                   int*    ncnln,
                                                                   int*    n,
                                                                   int*    ldJ,
                                                                   int*    needc,
                                                                   double* x,
                                                                   double* c,
                                                                   double* cJac,
                                                                   int*    nstate);


//extern
//libMesh::LibMeshInit     *__init;
//extern
MAST::FunctionEvaluation *__my_func_eval;





class StiffenedPlateThermallyStressedPistonTheorySizingOptimization:
        public MAST::FunctionEvaluation{ // this class inherits from Function evaluation

protected:      // protected member variables

    MAST::Examples::GetPotWrapper &_input;
    // number of elements and number of stations at which DVs are defined
    unsigned int
            _n_eig,
            _n_divs_x,
            _n_divs_between_stiff,
            _n_stiff, // assumed along x-axis with equidistance spacing
            _n_plate_elems,
            _n_elems_per_stiff,
            _n_elems,
            _n_dv_stations_x,
            _n_load_steps;

    libMesh::SerialMesh*                        _mesh;// create the mesh

    // section offset for plate
    MAST::SectionOffset *_hoff_plate_f;
    // Weight function to calculate the weight of the structure
    MAST::StiffenedPlateWeight *_weight;
    // create the material property card
    MAST::IsotropicMaterialPropertyCard *_m_card;

    // create the Dirichlet boundary condition on left edge
    MAST::DirichletBoundaryCondition *_dirichlet_left;
    // create the Dirichlet boundary condition on right edge
    MAST::DirichletBoundaryCondition *_dirichlet_right;
    // create the Dirichlet boundary condition on bottom edge
    MAST::DirichletBoundaryCondition *_dirichlet_bottom;
    // create the Dirichlet boundary condition on top edge
    MAST::DirichletBoundaryCondition *_dirichlet_top;
    // create the temperature load
    MAST::BoundaryConditionBase
            *_T_load,
            *_p_load;
    // piston theory boundary condition for the whole domain
    MAST::PistonTheoryBoundaryCondition *_piston_bc;
    //   piston theory boundary condition for the whole domain

    // -------------------------------------------------------------------
    // create the element property card, one for each stiffener
    std::vector<MAST::Solid1DSectionElementPropertyCard *> _p_card_stiff;
    // section offset for each stiffener
    std::vector<MAST::SectionOffset *>
            _hyoff_stiff_f,
            _hzoff_stiff_f;
    // length of domain
    Real _length;
    // width of domain
    Real _width;
    // maximum stress limit
    Real _stress_limit;
    // minimum flutter speed
    Real _V0_flutter;
    //   scaling parameters for design optimization problem
    std::vector<Real>
            _dv_scaling,
            _dv_low,
            _dv_init;
    //   interpolates thickness between stations
    std::vector<MAST::MultilinearInterpolation *>
            _thy_stiff_f,
            _thz_stiff_f; // one per stiffener
    // vector of basis vectors from modal analysis
    std::vector<libMesh::NumericVector<Real> *> _basis;
    // stationwise function objects for thickness
    std::vector<MAST::ConstantFieldFunction *>
            _th_station_functions_plate,
            _thy_station_functions_stiff,
            _thz_station_functions_stiff; // N_stiff*N_stations functions
    // stationwise parameter definitions
    std::vector<MAST::Parameter *>
            _th_station_parameters_plate,
            _thy_station_parameters_stiff,
            _thz_station_parameters_stiff; // N_stiff*N_stations params
    std::vector<MAST::Parameter*>                _problem_parameters;
    MAST::MultilinearInterpolation* _th_plate_f;
    // number of velocity divisions from V=0 to V=2*V0_flutter
    unsigned int _n_V_divs_flutter;
    std::vector<MAST::StressStrainOutputBase *> output_vec;
    //type of shape fuctions
    libMesh::FEType                           _fetype;
    Real _dx;
    bool _if_vk;
    MAST::Parameter* _zero;
    MAST::ConstantFieldFunction* _thzoff_stiff_f;
    std::map<Real, MAST::FieldFunction<Real> *> _thy_station_vals,
            _thz_station_vals;
    MAST::StructuralNearNullVectorSpace*       _nsp;

    // create the property functions and add them to the
    MAST::Parameter
            *_mach,
            *_rho_air,
            *_gamma_air;

    MAST::ConstantFieldFunction
            *_velocity_f,
            *_mach_f,
            *_rho_air_f,
            *_gamma_air_f;

    MAST::NonlinearSystem*                      _sys;// create the libmesh system
    MAST::Parameter
            *_temp,
            *_p_cav,
            *_velocity;
    bool _initialized;

    // create the element property card for the plate
    MAST::Solid2DSectionElementPropertyCard *_p_card_plate;
    MAST::StructuralSystemInitialization*       _structural_sys;// initialize the system to the right set of variables
    MAST::PhysicsDisciplineBase*                _discipline;
    MAST::FieldFunction<Real>*                  _jac_scaling;

    MAST::NonlinearImplicitAssembly*            _nonlinear_assembly;// nonlinear assembly object
   MAST::EigenproblemAssembly*                 _modal_assembly;// nonlinear assembly object
    MAST::StructuralFluidInteractionAssembly*   _fsi_assembly;// nonlinear assembly object

    MAST::TimeDomainFlutterSolver*              _flutter_solver;

    libMesh::EquationSystems*                   _eq_sys;// create the equation system
    MAST::StructuralNonlinearAssemblyElemOperations*          _nonlinear_elem_ops;

    MAST::StressAssembly*                             _stress_assembly;
    MAST::StructuralModalEigenproblemAssemblyElemOperations*  _modal_elem_ops;
    Real                                      _p_val, _vm_rho;
    bool _if_continuation_solver;
    MAST::StressStrainOutputBase*               _stress_elem;
    // output quantity objects to evaluate stress
    std::vector<MAST::StressStrainOutputBase *> _outputs;
    bool _if_neg_eig;
public:  // parametric constructor
    StiffenedPlateThermallyStressedPistonTheorySizingOptimization(const libMesh::Parallel::Communicator &comm,
                                                                  MAST::Examples::GetPotWrapper& input) :
            MAST::FunctionEvaluation (comm),
            _input(input),
            _n_eig(0),
            _n_divs_x(0),
            _n_divs_between_stiff(0),
            _n_stiff(0),
            _n_plate_elems(0),
            _n_elems_per_stiff(0),
            _n_elems(0),
            _n_dv_stations_x(0),
            _n_load_steps(0),
            _mesh(nullptr),
            
            _hoff_plate_f(nullptr),
            
            
            _weight(nullptr),
            _m_card(nullptr),
            
            _dirichlet_left(nullptr),
            _dirichlet_right(nullptr),
            _dirichlet_bottom(nullptr),
            _dirichlet_top(nullptr),
            
            _T_load(nullptr),
            _p_load(nullptr),
            _piston_bc(nullptr),
            
            
            
            _length(0.),
            _width(0.),
            
            _stress_limit(0.),
            _V0_flutter(0.),
            _dx(0.),
            
            
            _if_vk(false),
            _zero(nullptr),
            _thzoff_stiff_f(nullptr),
            _nsp(nullptr),
            _mach(nullptr),
            _rho_air(nullptr),
            _gamma_air(nullptr),
            
            
            _velocity_f(nullptr),
            
            _mach_f(nullptr),
            _rho_air_f(nullptr),
            _gamma_air_f(nullptr),
            _sys(nullptr),
            _temp(nullptr),
            _p_cav(nullptr),
            
            _velocity(nullptr),
            
            _initialized(false),
            
            _p_card_plate(nullptr),
            _structural_sys(nullptr),
            _discipline(nullptr),
            _jac_scaling(nullptr),
            _nonlinear_assembly(nullptr),
            _modal_assembly(nullptr),
            
            _fsi_assembly(nullptr),
            _flutter_solver(nullptr),
            _eq_sys(nullptr),
            _nonlinear_elem_ops(nullptr),
            _stress_assembly(nullptr),
            _modal_elem_ops(nullptr),
            
            _p_val(0.),
            _vm_rho(0.),
            _if_continuation_solver(false),
            _stress_elem(nullptr),
            _if_neg_eig(false)

//        _thy_station_vals(nullptr),
//        _thz_station_vals(nullptr)
    {
        libmesh_assert(!_initialized);

        _if_vk = _input("nonlinear", "linear or nonlinear", false);
        _if_continuation_solver = _input( "if_continuation_solver", "continuation solver on/off flag",  false);

        libMesh::out
                << "//////////////////////////////////////////////////////////////////" << std::endl
                << " stiffened_plate_thermally_stressed_piston_theory_optimization" << std::endl
                <<  std::endl
                << " input.in should be provided in the working directory with"
                << " desired parameter values."
                << " In absence of a parameter value, its default value will be used." << std::endl
                << std::endl
                << " Output per iteration is written to optimization_output.txt."      << std::endl
                << "//////////////////////////////////////////////////////////////////" << std::endl
                << std::endl;

        // number of load steps
        _n_load_steps = _input("n_load_steps", " ", 10);
        // number of elements
        _n_divs_x = _input("n_divs_x", " ", 40);
        _n_divs_between_stiff = _input("n_divs_between_stiff", " ", 10);
        _n_stiff = _input("n_stiffeners", " ", 3);

        // number of stations
        _n_dv_stations_x = _input("n_stations", " ", 3);

        _init_fetype(); // we initialize the type of fem (1st ord lagrange default)
        _init_mesh();

        // number of eigenvalues
        _n_eig = 20;
        // now setup the optimization data
        _n_vars = _n_dv_stations_x + 2 * _n_dv_stations_x * _n_stiff; // for thickness variable
        _n_eq = 0;
        _n_ineq = _n_eig + 1 +
                  _n_elems;// constraint that each eigenvalue > 0 flutter constraint + one element stress functional per elem

        _max_iters = 1000;

        // limit stress
        _stress_limit = _input("max_stress", " ", 4.00e8);

        // velocity constraint for flutter analysis
        _V0_flutter = _input("V0_flutter", "", 410.);
        _n_V_divs_flutter = _input("n_V_divs", "", 10);

        _p_val   = _input("constraint_aggregation_p_val", "value of p in p-norm stress aggregation", 2.0);
        _vm_rho  = _input("constraint_aggregation_rho_val", "value of rho in p-norm stress aggregation", 2.0);


        // call the initialization routines for each component

        _init_system_and_discipline();

        _init_dirichlet_conditions();

        _init_eq_sys();

        _init_dv_vector();

        _init_material(); // create the property functions for the plate

        _init_loads();

        _init_piston_vals();

        _init_thickness_variables_plate();

        _init_thickness_variables_stiff();

        _init_nullspace(); // initialize the null space object and assign it to the structural module

        _init_outputs();

        _initialized = true;

        std::ofstream output;
        output.open("optimization_output.txt", std::ofstream::out);

        // move to constructor
        _nonlinear_assembly = new MAST::NonlinearImplicitAssembly;// nonlinear assembly object
        _nonlinear_elem_ops = new MAST::StructuralNonlinearAssemblyElemOperations;

        _modal_assembly     = new MAST::EigenproblemAssembly;// nonlinear assembly object
        _modal_elem_ops     = new MAST::StructuralModalEigenproblemAssemblyElemOperations;

        _stress_assembly    = new MAST::StressAssembly;
        _stress_elem        = new MAST::StressStrainOutputBase;


        _fsi_assembly       = new MAST::StructuralFluidInteractionAssembly;// nonlinear assembly object

        _flutter_solver = new MAST::TimeDomainFlutterSolver;
        std::string nm("flutter_output.txt");
        if (this->comm().rank() == 0)
            _flutter_solver->set_output_file(nm);



        // create the function to calculate weight
        _weight = new MAST::StiffenedPlateWeight(*_discipline);


    }

    ~StiffenedPlateThermallyStressedPistonTheorySizingOptimization() {

        if (_initialized) {

            delete _m_card;
            delete _p_card_plate;
            for (unsigned int i = 0; i < _n_stiff; i++) delete _p_card_stiff[i];

            delete _T_load;
            delete _p_load;

            delete _dirichlet_bottom;
            delete _dirichlet_right;
            delete _dirichlet_top;
            delete _dirichlet_left;


            delete _hoff_plate_f;

            delete _thzoff_stiff_f;
            for (unsigned int i = 0; i < _n_stiff; i++) delete _hyoff_stiff_f[i];
            for (unsigned int i = 0; i < _n_stiff; i++) delete _thy_stiff_f[i];
            for (unsigned int i = 0; i < _n_stiff; i++) delete _thz_stiff_f[i];

            delete _jac_scaling;

            delete _zero;

            delete _weight;

            delete _nonlinear_assembly;
            //delete _modal_assembly;
            delete _fsi_assembly;

            // delete the basis vectors
            if (_basis.size())
                for (unsigned int i = 0; i < _basis.size(); i++)
                    delete _basis[i];

            delete _eq_sys;
            delete _mesh;

            delete _discipline;
            delete _structural_sys;


            //delete _flutter_solver;
            delete _piston_bc;

            delete _nsp;


            // iterate over the output quantities and delete them
            {
                std::vector<MAST::StressStrainOutputBase *>::iterator
                        it = _outputs.begin(),
                        end = _outputs.end();
                for (; it != end; it++) delete *it;

                _outputs.clear();
            }


            // delete the h_y station functions
            {
                std::vector<MAST::ConstantFieldFunction *>::iterator
                        it = _th_station_functions_plate.begin(),
                        end = _th_station_functions_plate.end();
                for (; it != end; it++) delete *it;
            }


            {
                std::vector<MAST::ConstantFieldFunction *>::iterator
                        it = _thy_station_functions_stiff.begin(),
                        end = _thy_station_functions_stiff.end();
                for (; it != end; it++) delete *it;
            }

            {
                std::vector<MAST::ConstantFieldFunction *>::iterator
                        it = _thz_station_functions_stiff.begin(),
                        end = _thz_station_functions_stiff.end();
                for (; it != end; it++) delete *it;
            }


            // delete the h_y station parameters
            {
                std::vector<MAST::Parameter *>::iterator
                        it = _th_station_parameters_plate.begin(),
                        end = _th_station_parameters_plate.end();
                for (; it != end; it++) delete *it;
            }


            {
                std::vector<MAST::Parameter *>::iterator
                        it = _thy_station_parameters_stiff.begin(),
                        end = _thy_station_parameters_stiff.end();
                for (; it != end; it++) delete *it;
            }

            {
                std::vector<MAST::Parameter *>::iterator
                        it = _thz_station_parameters_stiff.begin(),
                        end = _thz_station_parameters_stiff.end();
                for (; it != end; it++) delete *it;
            }
        }
    }





    void _init_fetype() {

        // FEType to initialize the system. Get the order and type of element.
        std::string
                order_str   = _input("fe_order", "order of finite element shape basis functions",    "first"),
                family_str  = _input("fe_family",      "family of finite element shape functions", "lagrange");

        libMesh::Order
                o  = libMesh::Utility::string_to_enum<libMesh::Order>(order_str);
        libMesh::FEFamily
                fe = libMesh::Utility::string_to_enum<libMesh::FEFamily>(family_str);
        _fetype = libMesh::FEType(o, fe);
    }

    void _init_mesh() {

        // create the mesh

        _mesh = new libMesh::SerialMesh(this->comm());

        // length of domain
        _length = _input("length", " ", 0.50);
        _width = _input("width", " ", 0.25);

        // identify the element type from the input file or from the order
        // of the element
        std::string
                t = _input("elem_type", "type of geometric element in the mesh", "quad4");
        libMesh::ElemType
                e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);




        _n_plate_elems = _n_divs_x * (_n_stiff + 1) * _n_divs_between_stiff;

        if (e_type == libMesh::TRI3)
            _n_plate_elems *= 2;

        _n_elems_per_stiff = _n_divs_x;
        _n_elems = _n_plate_elems + _n_stiff * _n_elems_per_stiff;



        // initialize the mesh with one element
        MAST::StiffenedPanelMesh panel_mesh;
        panel_mesh.init(_n_stiff,
                        _n_divs_x,
                        _n_divs_between_stiff,
                        _length,
                        _width,
                        *_mesh,
                        e_type,
                        true);
        libMesh::out
                << "//////////////////////////////////////////////////////////////////" << std::endl
                                                                                        << std::endl;
        _mesh->print_info();
        libMesh::out
                << "//////////////////////////////////////////////////////////////////" << std::endl
                                                                                        << std::endl;

    }

    void _init_system_and_discipline() {
        //
        // make sure that the mesh has been initialized
        //
        libmesh_assert(_mesh);

        // create the equation system
        _eq_sys = new libMesh::EquationSystems(*_mesh);

        // create the libmesh system
        _sys = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));

        // specifying the type of eigenproblem we'd like to solve
        _sys->set_eigenproblem_type(libMesh::GHEP);

        // initialize the system to the right set of variables
        _structural_sys = new MAST::StructuralSystemInitialization(*_sys,
                                                                   _sys->name(),
                                                                   _fetype);

        _discipline = new MAST::PhysicsDisciplineBase(*_eq_sys);



    }

    void _init_dirichlet_conditions() {

        // create and add the boundary condition and loads
        // 0 1 2 3 stands for the edges of the plate
        // _structural_sys->vars() gives the vars indexes (u,v,w,thx.thy.thz)



//        std::vector<unsigned int> constrained_vars(4);
//        // not constraning ty, tz will keep it simply supported
//        constrained_vars[0] = 0;  // u
//        constrained_vars[1] = 1;  // v
//        constrained_vars[2] = 2;  // w
//        constrained_vars[3] = 5;  // tz

        //_dirichlet_bottom->init (0, constrained_vars);
        //_dirichlet_right->init  (1, constrained_vars);
        //_dirichlet_top->init    (2, constrained_vars);
        //_dirichlet_left->init   (3, constrained_vars);


        // create and add the boundary condition
        _dirichlet_bottom = new MAST::DirichletBoundaryCondition;
        _dirichlet_right = new MAST::DirichletBoundaryCondition;
        _dirichlet_top = new MAST::DirichletBoundaryCondition;
        _dirichlet_left = new MAST::DirichletBoundaryCondition;

        _dirichlet_bottom->init(0, _structural_sys->vars());
        _dirichlet_right->init(1, _structural_sys->vars());
        _dirichlet_top->init(2, _structural_sys->vars());
        _dirichlet_left->init(3, _structural_sys->vars());

        _discipline->add_dirichlet_bc(0, *_dirichlet_bottom);
        _discipline->add_dirichlet_bc(1, *_dirichlet_right);
        _discipline->add_dirichlet_bc(2, *_dirichlet_top);
        _discipline->add_dirichlet_bc(3, *_dirichlet_left);

        _discipline->init_system_dirichlet_bc(*_sys);
    }

    void _init_eq_sys() {

        // initialize the equation system
        _eq_sys->init();

        //The EigenSolver, definig which interface, i.e solver package to use.
        _sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);

        //sets the flag to exchange the A and B matrices for a generalized eigenvalue problem.
        //This is needed typically when the B matrix is not positive semi-definite.
        _sys->set_exchange_A_and_B(true);

        //sets the number of eigenvalues requested
        _sys->set_n_requested_eigenvalues(_n_eig);

        //Loop over the dofs on each processor to initialize the list of non-condensed dofs.
        //These are the dofs in the system that are not contained in global_dirichlet_dofs_set.
               _sys->initialize_condensed_dofs(*_discipline);



    }

    void _init_dv_vector() {
        // initialize the dv vector data (dv = design variables)
        const Real
                th_l     = _input("thickness_lower", "", 0.001),
                th_u     = _input("thickness_upper", "", 0.2),
                th       = _input("thickness", "", 0.2),
		        th_stiff = _input("thickness_stiff","",0.2);

        //distance btw stations
        _dx = _length / (_n_dv_stations_x - 1);

        _dv_init.resize(_n_vars);
        _dv_scaling.resize(_n_vars);
        _dv_low.resize(_n_vars);

        _problem_parameters.resize(_n_vars);

        // design variables for the thickness values
        for (unsigned int i = 0; i < _n_vars; i++) {

            //_dv_init[i] = _input("dv_init", "", th / th_u, i);
            _dv_low[i] = th_l / th_u;
            _dv_scaling[i] = th_u;
        }

        for (unsigned int i = 0; i < _n_dv_stations_x; i++) {
            _dv_init[i] = _input("dv_init", "", th / th_u, i);
        }


        for (unsigned int i = _n_dv_stations_x; i < _n_vars; i++) {
            _dv_init[i] = _input("dv_init", "", th_stiff / th_u, i);
        }

    }

    void _init_thickness_variables_plate(){

        // create the thickness variables
        _th_station_parameters_plate.resize(_n_dv_stations_x);
        _th_station_functions_plate.resize(_n_dv_stations_x);


        for (unsigned int i = 0; i < _n_dv_stations_x; i++) {
            std::ostringstream oss;
            oss << "h_" << i;

            // now we need a parameter that defines the thickness at the
            // specified station and a constant function that defines the
            // field function at that location.
            MAST::Parameter *h =
                    new MAST::Parameter(oss.str(), _input("thickness", "", 0.002));

            MAST::ConstantFieldFunction *h_f =
                    new MAST::ConstantFieldFunction(oss.str(), *h);

            // add this to the thickness map
            _thy_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real> *>
                                             (i * _dx, h_f));

            // add the function to the parameter set
            _th_station_parameters_plate[i] = h;
            _th_station_functions_plate[i] = h_f;

            // tell the assembly system about the sensitvity parameter

            //_discipline->add_parameter(*h); no longer needed
            _problem_parameters[i] = h;
        }

        // now create the h_y function and give it to the property card
        _th_plate_f = new MAST::MultilinearInterpolation("h", _thy_station_vals);
        _thy_station_vals.clear();

        _hoff_plate_f   = new MAST::SectionOffset("off",
                                                  *_th_plate_f,
                                                  0.);

        // create the element property card
        _p_card_plate = new MAST::Solid2DSectionElementPropertyCard;

        // add the section properties to the card
        _p_card_plate->add(*_th_plate_f);
        _p_card_plate->add(*_hoff_plate_f);

        // tell the section property about the material property
        _p_card_plate->set_material(*_m_card);

        if (_if_vk) _p_card_plate->set_strain(MAST::NONLINEAR_STRAIN);

        _discipline->set_property_for_subdomain(0, *_p_card_plate);

    }

    void _init_material() {
        // create the property functions and add them to the card


        // create the property functions and add them to the
        MAST::Parameter
                *E,
                *alpha,
                *nu,
                *kappa,
                *rho;
        MAST::ConstantFieldFunction
                *E_f,
                *nu_f,
                *alpha_f,
                *kappa_f,
                *rho_f;


        E = new MAST::Parameter("E", _input("E", "", 72.e9));
        nu = new MAST::Parameter("nu", _input("nu", "", 0.33));
        kappa = new MAST::Parameter("kappa", _input("kappa", "", 5. / 6.));
        alpha = new MAST::Parameter("alpha", _input("alpha", "", 2.5e-5));
        rho = new MAST::Parameter("rho", _input("rho", "", 2700.0));

        E_f = new MAST::ConstantFieldFunction("E", *E);
        nu_f = new MAST::ConstantFieldFunction("nu", *nu);
        kappa_f = new MAST::ConstantFieldFunction("kappa", *kappa);
        alpha_f = new MAST::ConstantFieldFunction("alpha_expansion", *alpha);
        rho_f = new MAST::ConstantFieldFunction("rho", *rho);

        // create the material property card
        _m_card = new MAST::IsotropicMaterialPropertyCard;

        // add the material properties to the card
        _m_card->add(*E_f);
        _m_card->add(*nu_f);
        _m_card->add(*rho_f);
        _m_card->add(*kappa_f);
        _m_card->add(*alpha_f);

    }

    void _init_thickness_variables_stiff() {

        // thickness per stiffener station
        _thy_station_parameters_stiff.resize(_n_dv_stations_x * _n_stiff);
        _thy_station_functions_stiff.resize(_n_dv_stations_x * _n_stiff);
        _thz_station_parameters_stiff.resize(_n_dv_stations_x * _n_stiff);
        _thz_station_functions_stiff.resize(_n_dv_stations_x * _n_stiff);
        _thy_stiff_f.resize(_n_stiff);
        _thz_stiff_f.resize(_n_stiff);
        _hyoff_stiff_f.resize(_n_stiff);

        for (unsigned int i = 0; i < _n_stiff; i++) {

            // this map is used to store the thickness parameter along length
            _thy_station_vals.clear();
            _thz_station_vals.clear();

            // first define the thickness station parameters and the thickness
            // field function
            for (unsigned int j = 0; j < _n_dv_stations_x; j++) {
                std::ostringstream ossy, ossz;
                ossy << "h_y_" << j << "_stiff_" << i;
                ossz << "h_y_" << j << "_stiff_" << i;

                // now we need a parameter that defines the thickness at the
                // specified station and a constant function that defines the
                // field function at that location.
                MAST::Parameter
                        *h_y = new MAST::Parameter(ossy.str(), _input("thickness_stiff", "", 0.002)),
                        *h_z = new MAST::Parameter(ossz.str(), _input("thickness_stiff", "", 0.002));

                MAST::ConstantFieldFunction
                        *h_y_f = new MAST::ConstantFieldFunction(ossy.str(), *h_y),
                        *h_z_f = new MAST::ConstantFieldFunction(ossy.str(), *h_z);

                // add this to the thickness map
                _thy_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real> *>
                                                 (j * _dx, h_y_f));
                _thz_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real> *>
                                                 (j * _dx, h_z_f));

                // add the function to the parameter set
                _thy_station_parameters_stiff[i * _n_dv_stations_x + j] = h_y;
                _thy_station_functions_stiff[i * _n_dv_stations_x + j] = h_y_f;
                _thz_station_parameters_stiff[i * _n_dv_stations_x + j] = h_z;
                _thz_station_functions_stiff[i * _n_dv_stations_x + j] = h_z_f;


                // tell the assembly system about the sensitvity parameter
                //_discipline->add_parameter(*h_y);
                //_discipline->add_parameter(*h_z);
                _problem_parameters[(2 * i + 1) * _n_dv_stations_x + j] = h_y;
                _problem_parameters[(2 * i + 2) * _n_dv_stations_x + j] = h_z;
            }


            _thzoff_stiff_f = new MAST::ConstantFieldFunction("hz_off", *_zero);
            RealVectorX orientation = RealVectorX::Zero(3);
            orientation(1) = 1.;

            // property card per stiffener
            _p_card_stiff.resize(_n_stiff);


                // now create the h_y function and give it to the property card
                _thy_stiff_f[i] = new MAST::MultilinearInterpolation("hy", _thy_station_vals);
                _thz_stiff_f[i] = new MAST::MultilinearInterpolation("hz", _thz_station_vals);
                _hyoff_stiff_f[i] = new MAST::SectionOffset("hy_off",
                                                            *_thy_stiff_f[i],
                                                            -1.);

                _p_card_stiff[i] = new MAST::Solid1DSectionElementPropertyCard;


                // add the section properties to the card
                _p_card_stiff[i]->add(*_thy_stiff_f[i]);
                _p_card_stiff[i]->add(*_thz_stiff_f[i]);
                _p_card_stiff[i]->add(*_hyoff_stiff_f[i]);
                _p_card_stiff[i]->add(*_thzoff_stiff_f);
                _p_card_stiff[i]->y_vector() = orientation;

                // tell the section property about the material property
                _p_card_stiff[i]->set_material(*_m_card);

                //_p_card_stiff[i]->set_bending_model(MAST::TIMOSHENKO);
                //_p_card_stiff[i]->set_bending_model(MAST::BERNOULLI);

                if (_if_vk) _p_card_stiff[i]->set_strain(MAST::NONLINEAR_STRAIN);

                _p_card_stiff[i]->init();

                // the domain ID of the stiffener is 1 plus the stiff number
                _discipline->set_property_for_subdomain(i + 1, *_p_card_stiff[i]);
            }
    }

    void _init_piston_vals() {
        // create the property functions and add them to the

        _velocity = new MAST::Parameter("V", 0.);
        _mach = new MAST::Parameter("mach", _input("mach", "", 4.5));
        _rho_air = new MAST::Parameter("rho", _input("rho_f", "", 1.05));
        _gamma_air = new MAST::Parameter("gamma", _input("gamma", "", 1.4));


        _velocity_f = new MAST::ConstantFieldFunction("V", *_velocity);
        _mach_f = new MAST::ConstantFieldFunction("mach", *_mach);
        _rho_air_f = new MAST::ConstantFieldFunction("rho", *_rho_air);
        _gamma_air_f = new MAST::ConstantFieldFunction("gamma", *_gamma_air);

        // now initialize the piston theory boundary conditions
        RealVectorX vel = RealVectorX::Zero(3);
        vel(0) = 1.;  // flow along the x-axis
        _piston_bc = new MAST::PistonTheoryBoundaryCondition(1,     // order
                                                             vel);  // vel vector
        _piston_bc->add(*_velocity_f);
        _piston_bc->add(*_mach_f);
        _piston_bc->add(*_rho_air_f);
        _piston_bc->add(*_gamma_air_f);
        //_discipline->add_volume_load(0, *_piston_bc);
        //_discipline->add_parameter(*_velocity); no longer needed// //
    }

    void _init_loads() {

        // create the temperature load
        // create the property functions and add them to the

        MAST::ConstantFieldFunction
                *ref_temp_f,
                *temp_f,
                *p_cav_f;



        _p_cav = new MAST::Parameter("p_cav", _input("p_cav", "", -300.));
        p_cav_f = new MAST::ConstantFieldFunction("pressure", *_p_cav);

        _zero = new MAST::Parameter("zero", 0.);

        ref_temp_f = new MAST::ConstantFieldFunction("ref_temperature", *_zero);
        _temp       = new MAST::Parameter("temperature", _input("temp", "", 10.));
        temp_f     = new MAST::ConstantFieldFunction("temperature", *_temp);


            // initialize the load
            _jac_scaling = new MAST::Examples::ThermalJacobianScaling;


        _T_load = new MAST::BoundaryConditionBase(MAST::TEMPERATURE);
        _T_load->add(*temp_f);
        _T_load->add(*ref_temp_f);

        if(!_if_continuation_solver)
            _T_load->add(*_jac_scaling);


        _discipline->add_volume_load(0, *_T_load);          // for the panel
        for (unsigned int i = 0; i < _n_stiff; i++)
            _discipline->add_volume_load(i + 1, *_T_load);    // for the stiffeners

        // pressure load
        _p_load = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
        _p_load->add(*p_cav_f);
        _discipline->add_volume_load(0, *_p_load);          // for the panel
    }

    void _init_nullspace(){

        _nsp = new MAST::StructuralNearNullVectorSpace;
        _sys->nonlinear_solver->nearnullspace_object = _nsp;
    }

    void _init_outputs(){
        // create the output objects, one for each element

            _outputs.resize(_mesh->n_elem(), nullptr);

            for (int i = 0; i < _mesh->n_elem(); i++) {
                MAST::StressStrainOutputBase *output = new MAST::StressStrainOutputBase;
                output->set_discipline_and_system(*_discipline, *_structural_sys);
                output->set_aggregation_coefficients(_p_val, 1.0, _vm_rho, _stress_limit);
                _outputs[i] = output;
        }


        std::vector<unsigned int> n_size(comm().size(), 0);

        for (int i = 0; i < comm().size(); i++) {
            n_size[i] = _mesh->n_elem_on_proc(i);
        }

          //  std::cout << "number of elems on proc " << 0 << " is  " << n_size[0] << std::endl;
           // std::cout << "number of elems on proc " << 1 << " is  " << n_size[1] << std::endl;
        //       n[comm().rank()] = _mesh->n_local_elem();
        //     comm().sum(n);

        std::vector<unsigned int> prev_elems(comm().size() , 0);

        if (comm().rank() != 0) {
            for (int j = 0; j < comm().rank() ; j++)
                prev_elems[comm().rank()] += n_size[j];
        }

        comm().sum(prev_elems);


        libMesh::MeshBase::const_element_iterator
                e_it = _mesh->local_elements_begin(),
                e_end = _mesh->local_elements_end();

        int i = 0;
        for (; e_it != e_end; e_it++) {
            _outputs[ (i + prev_elems[comm().rank()]) ]->set_participating_elements({*e_it});

            i += 1;
        }

        if (_outputs.size() != _mesh->n_elem())
            libMesh::out << "_outputs is not the correct size " << std::endl;

            libmesh_assert_equal_to(_outputs.size(), _mesh->n_elem());


    }

    virtual void init_dvar(std::vector<Real>& x,
                           std::vector<Real>& xmin,
                           std::vector<Real>& xmax) {

        // one DV for each element

       x.resize(_n_vars);
       xmin.resize(_n_vars);
       xmax.resize(_n_vars);

       xmin    = _dv_low;
       xmax    = _dv_scaling;

        //
        // now, check if the user asked to initialize dvs from a previous file
        //
        std::string
                nm    =  _input("restart_optimization_file", "filename with optimization history for restart", "");

        if (nm.length()) {

            unsigned int
                    iter = _input("restart_optimization_iter", "restart iteration number from file", 0);
            this->initialize_dv_from_output_file(nm, iter, x);

            libmesh_assert_equal_to(x.size(),_n_vars);
        }
        else {
                x       = _dv_init;
        }
    }

    virtual void evaluate(const std::vector<Real> &dvars, // design variable
                          Real &obj,                      // objective function
                          bool eval_obj_grad,             // flag to evaluate or not grad of obj func
                          std::vector<Real> &obj_grad,    // the gradient of the obj func
                          std::vector<Real> &fvals,       // constraint functions
                          std::vector<bool> &eval_grads,  // vector of flags to evaluate gradient of constraint func
                          std::vector<Real> &grads) {     // gradients of constraint functions

        libmesh_assert(_initialized);
        libmesh_assert_equal_to(dvars.size(), _n_vars);



        //MAST::EigenproblemAssembly                               modal_assembly;
        //MAST::StructuralModalEigenproblemAssemblyElemOperations  modal_elem_ops;


        // set the parameter values equal to the DV value
        // first the plate thickness values
        for (unsigned int i = 0; i < _n_vars; i++)
            (*_problem_parameters[i]) = dvars[i] * _dv_scaling[i];


        // DO NOT zero out the gradient vector, since GCMMA needs it for the
        // subproblem solution
        // zero the function evaluations
        std::fill(fvals.begin(), fvals.end(), 0.);


        libMesh::Point pt; // dummy point object

        // print the values of all design variables
        libMesh::out << "//////////////////////////////////////////////////////////////////////"<< std::endl
                     << " New Eval " << std::endl;
        for (unsigned int i = 0; i < _n_vars; i++)
            libMesh::out
                    << "th     [ " << std::setw(10) << i << " ] = "
                    << std::setw(20) << (*_problem_parameters[i])() << std::endl;

        bool
        if_write_output = _input("if_write_output", "print outputs", false);

        
        //////////////////////////////////////////////////////////////////////
        libMesh::out << " calculation of wheight " << std::endl;
        // the optimization problem is defined as
        // min weight, subject to constraints on displacement and stresses
        Real
                wt = 0.;
        // calculate weight
        (*_weight)(pt, 0., wt);
        libMesh::out << " wheight calculated " << std::endl;


        // solve for the steady state at zero velocity
        (*_velocity) = 0.;
        //////////////////////////////////////////////////////////////////////
        // steady state solution
        StiffenedPlateSteadySolverInterface steady_solve(*this,
                                                         if_write_output,
                                                         false,
                                                         _n_load_steps);

        if (_if_vk) {
            if (!_if_continuation_solver) {
                libMesh::out << "** Steady state solution w/ Newton raphson solver **" << std::endl;
                _if_neg_eig = false;
            }
            else {
                libMesh::out << "** Steady state solution w/ continuation solver **" << std::endl;
                _if_neg_eig = false;
            }
        } else{
            libMesh::out << "** Steady state solution w/ Newton raphson solver **" << std::endl;
        }

        steady_solve.solve();


        // us this solution as the base solution later if no flutter is found.
        libMesh::NumericVector<Real> &
                steady_sol_wo_aero = _sys->add_vector("steady_sol_wo_aero");
        steady_sol_wo_aero = _sys->get_vector("base_solution");


        // now that we have the solution under the influence of thermal loads,
        // we will use a small number of load steps to get to the steady-state
        // solution

        //steady_solve.set_n_load_steps(50);
        //steady_solve.set_modify_only_aero_load(true);

        if (_if_neg_eig) {
            obj = 1.e10;
            for (unsigned int i=0; i<_n_ineq; i++)
                fvals[i] = 1.e10;
            return;
        }

        //////////////////////////////////////////////////////////////////////
        // perform the modal analysis

        // modal analysis is about the base state exclusing aerodynamic loads.
        // So, they will act as generalized coordinates that will not provide
        // diagonal reduced order mass/stiffness operator. The eigenvalues
        // will also be independent of velocity.

        libMesh::out << "** modal analysis **" << std::endl;

        _modal_assembly->set_discipline_and_system(*_discipline, *_structural_sys); // modf_w
        _modal_assembly->set_base_solution(steady_sol_wo_aero);
        _modal_elem_ops->set_discipline_and_system(*_discipline, *_structural_sys);
        _sys->eigen_solver->set_position_of_spectrum( libMesh::LARGEST_MAGNITUDE);
        _sys->eigenproblem_solve( *_modal_elem_ops, *_modal_assembly);
        _modal_assembly->clear_base_solution();
        _modal_assembly->clear_discipline_and_system();
        _modal_elem_ops->clear_discipline_and_system();


        unsigned int
                nconv = std::min(_sys->get_n_converged_eigenvalues(),
                                 _sys->get_n_requested_eigenvalues());
        if (_basis.size() > 0)
            libmesh_assert(_basis.size() == nconv);
        else {
            _basis.resize(nconv);
            for (unsigned int i = 0; i < _basis.size(); i++)
                _basis[i] = nullptr;
        }

        // vector of eigenvalues
        std::vector<Real> eig_vals(nconv);

        bool if_all_eig_positive = true;

        libMesh::ExodusII_IO *
                writer = nullptr;

        if (if_write_output)
            writer = new libMesh::ExodusII_IO(*_mesh);

        for (unsigned int i = 0; i < nconv; i++) {

            // create a vector to store the basis
            if (_basis[i] == nullptr)
                _basis[i] = _sys->solution->zero_clone().release(); // what happens in this line ?

            // now write the eigenvalue
            Real
                    re = 0.,
                    im = 0.;
            _sys->get_eigenpair(i, re, im, *_basis[i]);

            libMesh::out
                    << std::setw(35) << std::fixed << std::setprecision(15)
                    << re << std::endl;

            eig_vals[i] = re;

            // keep the flag true or change to false
            if_all_eig_positive = (if_all_eig_positive && (re > 0.)) ? true : false;

            if (if_write_output) {

                // copy the solution for output
                (*_sys->solution) = *_basis[i];

                // We write the file in the ExodusII format.
                std::set<std::string> nm;
                nm.insert(_sys->name());
                writer->write_timestep("modes.exo",
                                       *_eq_sys,
                                       i + 1, //  time step
                                       i);    //  time
            }
        }

        if (if_write_output)
            *_sys->solution = steady_sol_wo_aero;
        
        //////////////////////////////////////////////////////////////////////
        // perform the flutter analysis
        //////////////////////////////////////////////////////////////////////

        // velocity should be set to zero for all residual calculations
        (*_velocity) = 0.;
        steady_solve.solution() = steady_sol_wo_aero;
        *_sys->solution = steady_sol_wo_aero;
        

        //////////////////////////////////////////////////////////////////////
        //  plot stress solution
        //////////////////////////////////////////////////////////////////////
        if (if_write_output) {

            _stress_elem->set_aggregation_coefficients(_p_val,1.,_vm_rho,_stress_limit);
            _stress_elem->set_participating_elements_to_all();
            _stress_elem->set_discipline_and_system(*_discipline,*_structural_sys);
            _stress_assembly->set_discipline_and_system(*_discipline,*_structural_sys);
            _stress_assembly->update_stress_strain_data(*_stress_elem, steady_sol_wo_aero);

            libMesh::out << "Writing output to : output.exo" << std::endl;

            //std::set<std::string> nm;
            //nm.insert(_sys->name());
            // write the solution for visualization
            libMesh::ExodusII_IO(*_mesh).write_equation_systems("output.exo",
                                                                *_eq_sys);//,&nm);

            _stress_elem->clear_discipline_and_system();
            _stress_assembly->clear_discipline_and_system();
        }

        //////////////////////////////////////////////////////////////////////
        // now calculate the stress output based on the velocity output
        //////////////////////////////////////////////////////////////////////

        _nonlinear_assembly->set_discipline_and_system(*_discipline, *_structural_sys);
        for (int i=0; i < _mesh->n_elem() ; i++){
            _nonlinear_assembly->calculate_output(steady_sol_wo_aero,*_outputs[i]);
        }
        _nonlinear_assembly->clear_discipline_and_system();


        //////////////////////////////////////////////////////////////////////
        // get the objective
        //////////////////////////////////////////////////////////////////////

        // set the function and objective values
        obj = wt;

        // parallel sum of the weight
        this->comm().sum(obj);


        //////////////////////////////////////////////////////////////////////
        // stress constraints
        //////////////////////////////////////////////////////////////////////

        // copy the element von Mises stress values as the functions
        for (unsigned int i = 0; i < _mesh->n_elem(); i++)
            fvals[_n_eig + 1 + i ] = -1. +   _outputs[i]->output_total() /
                                             _stress_limit;

        //////////////////////////////////////////////////////////////////////
        // evaluate the eigenvalue constraint
        //////////////////////////////////////////////////////////////////////

        if (nconv) {
            // set the eigenvalue constraints  -eig <= 0. scale
            // by an arbitrary 1/1.e7 factor
            for (unsigned int i = 0; i < nconv; i++)
                fvals[i] = -eig_vals[i] / 1.e7;
        }

        //////////////////////////////////////////////////////////////////////
        // evaluate the flutter constraint
        //////////////////////////////////////////////////////////////////////

            fvals[_n_eig+0]  =  -100.;

        //////////////////////////////////////////////////////////////////
        //   evaluate sensitivity if needed
        //////////////////////////////////////////////////////////////////

        // sensitivity of the objective function

        if (eval_obj_grad) {

            Real w_sens = 0.; // needs to be adress ?

            // set gradient of weight
            for (unsigned int i = 0; i < _n_vars; i++) {

                _weight->derivative(*_problem_parameters[i],
                                    pt,
                                    0.,
                                    w_sens);

                obj_grad[i] = w_sens * _dv_scaling[i];
            }

            // parallel sum
            this->comm().sum(obj_grad);
        }


        //////////////////////////////////////////////////////////////////////
        // check to see if the sensitivity of constraint is requested
        //////////////////////////////////////////////////////////////////////
        bool if_sens = false;
        for (unsigned int i = 0; i < eval_grads.size(); i++)
            if_sens = (if_sens || eval_grads[i]);


        if (if_sens) {

            libMesh::out << "** Sensitivity analysis **" << std::endl;
            // first initialize the gradient vector to zero
            std::fill(grads.begin(), grads.end(), 0.);

            //////////////////////////////////////////////////////////////////
            // indices used by GCMMA follow this rule:
            // grad_k = dfi/dxj  ,  where k = j*NFunc + i
            //////////////////////////////////////////////////////////////////


            // copy the solution to be used for sensitivity
            // If a flutter solution was found, then this depends on velocity.
            // Otherwise, it is independent of velocity
            *_sys->solution = steady_solve.solution();

            // first do sensitivity analysis wrt velocity, which is necessary for
            // flutter sensitivity.
            libMesh::NumericVector<Real> &dXdV = _sys->add_vector("sol_V_sens");
            std::vector<Real> dsigma_dV(_outputs.size(), 0.);

            // no flutter solution
            dXdV.zero();
            
            _modal_assembly->set_discipline_and_system(*_discipline, *_structural_sys); // modf_w
            _modal_assembly->set_base_solution(steady_sol_wo_aero);
            _modal_elem_ops->set_discipline_and_system(*_discipline, *_structural_sys);
            _nonlinear_assembly->set_discipline_and_system(*_discipline,*_structural_sys);
            _nonlinear_elem_ops->set_discipline_and_system(*_discipline,*_structural_sys);

            // we are going to choose to use one parametric sensitivity at a time
            for (unsigned int i = 0; i < _n_vars; i++) {
                libMesh::out << "design variable " << i << std::endl;

                // iterate over each dv and calculate the sensitivity
                libMesh::NumericVector<Real> &dXdp = _sys->add_sensitivity_solution(0);
                dXdp.zero();

                // sensitivity analysis
                _sys->sensitivity_solve(*_nonlinear_elem_ops,
                                        *_nonlinear_assembly,
                                        *_problem_parameters[i],
                                        true);

                dXdp = _sys->get_sensitivity_solution(0);

                for (unsigned int j = 0 ; j < _mesh->n_elem(); j++){
                    // evaluate sensitivity of the outputs
                    _nonlinear_assembly->calculate_output_direct_sensitivity(steady_sol_wo_aero,
                                                                             &dXdp,
                                                                             *(_problem_parameters[i]),
                                                                             *(_outputs[j])  );
                }

                // copy the sensitivity values in the output. This accounts for the
                // sensitivity of state wrt parameter. However, if a flutter root
                // was found, the state depends on velocity, which depends on the
                // parameter. Hence, the total sensitivity of stress constraint
                // would need to include the latter component, which was added
                // above.
                for (unsigned int j = 0; j < _mesh->n_elem(); j++) {
                    grads[(i * _n_ineq) + (j + _n_eig + 1 )] = _dv_scaling[i] / _stress_limit *
                                    _outputs[j]->output_sensitivity_total(*(_problem_parameters[i]));
                }


                // if all eigenvalues are positive, calculate at the sensitivity of
                // flutter velocity
                // if no root was found, then set the sensitivity to a zero value
                
                grads[(i * _n_ineq) + (_n_eig + 0)] = 0.;
                

                // calculate the sensitivity of the eigenvalues
                std::vector<Real> eig_sens(nconv,0.);

                if (nconv) {
                    _modal_assembly->set_base_solution(dXdp, true);
                    _sys->eigenproblem_sensitivity_solve(*_modal_elem_ops,
                                                         *_modal_assembly,
                                                         *_problem_parameters[i],
                                                         eig_sens);


                    _modal_assembly->clear_base_solution(true);

                    for (unsigned int j = 0; j < nconv; j++) {
                        grads[(i * _n_ineq) + j] = -_dv_scaling[i] * eig_sens[j] / 1.e7;
                    }
                }
            }

            _nonlinear_assembly->clear_discipline_and_system();
            _nonlinear_elem_ops->clear_discipline_and_system();
            _modal_assembly->clear_discipline_and_system();
            _modal_assembly->clear_base_solution();
            _modal_elem_ops->clear_discipline_and_system();

            libMesh::out << "** sensitivity analysis DONE **" << std::endl;

        }
    }

    void clear_stresss() {

        // iterate over the output quantities and delete them
        std::vector<MAST::StressStrainOutputBase *>::iterator
                it = _outputs.begin(),
                end = _outputs.end();

        for (; it != end; it++)
            (*it)->clear();
    }


    virtual void output(unsigned int iter,
                        const std::vector<Real> &x,
                        Real obj,
                        const std::vector<Real> &fval,
                        bool if_write_to_optim_file)  {

        libmesh_assert_equal_to(x.size(), _n_vars);

        // write the DVs in the physical dimension
        for (unsigned int i = 0; i < _n_vars; i++)
            libMesh::out
                    << "th     [ " << std::setw(10) << i << " ] = "
                    << std::setw(20) << (*_problem_parameters[i])() << std::endl;




        // write the solution for visualization
        *_sys->solution = _sys->get_vector("base_solution");

        _stress_elem->set_discipline_and_system(*_discipline,*_structural_sys);
        _stress_assembly->set_discipline_and_system(*_discipline,*_structural_sys);
        _stress_assembly->update_stress_strain_data(*_stress_elem, *_sys->solution);

        libMesh::out << "Writing output to : output.exo" << std::endl;
        libMesh::ExodusII_IO(*_mesh).write_equation_systems("output.exo",
                                                            *_eq_sys);

        _stress_elem->clear_discipline_and_system();
        _stress_assembly->clear_discipline_and_system();

        MAST::FunctionEvaluation::output(iter, x, obj, fval, if_write_to_optim_file);
    }



    class StiffenedPlateSteadySolverInterface:
            public MAST::FlutterSolverBase::SteadySolver {
    protected:
        /*!
         *   pointer to the object that hold all the solution data
         */
        StiffenedPlateThermallyStressedPistonTheorySizingOptimization& _obj;

        /*!
         *   flag to toggle output
         */
        bool _if_write_output;

        /*!
         *   number of nonliear load increment steps
         */
        unsigned int _n_steps;

        /*!
         *   only aero load is modified during nonlinear iterations
         */
        bool _if_only_aero_load_steps;


        /*!
         *   deletes the solution vector from system when the class is
         *   destructed unless this flag is false.
         */
        bool _if_clear_vector_on_exit;
        libMesh::ExodusII_IO* writer;
        unsigned int inc;

    public:
        StiffenedPlateSteadySolverInterface(StiffenedPlateThermallyStressedPistonTheorySizingOptimization& obj,
                                            bool if_output,
                                            bool if_clear_vector_on_exit,
                                            unsigned int n_steps):
                MAST::FlutterSolverBase::SteadySolver(),
                _obj(obj),
                _if_write_output(if_output),
                _n_steps(n_steps),
                _if_only_aero_load_steps(false),
                _if_clear_vector_on_exit(if_clear_vector_on_exit),
                inc (0) {
            // add vector in sys call it bas_solution will store disps
            _obj._sys->add_vector("base_solution");
            writer = new libMesh::ExodusII_IO(_obj._sys->get_mesh());
        }

        virtual ~StiffenedPlateSteadySolverInterface() {

            if (_if_clear_vector_on_exit)
                _obj._sys->remove_vector("base_solution");
            delete writer;
        }


        /*!
         *  solves for the steady state solution, and @returns
         *  a const-reference to the solution.
         */
        virtual const libMesh::NumericVector<Real>&
        solve() {

            // create a vector to store the solution
            libMesh::NumericVector<Real>& sol = _obj._sys->get_vector("base_solution");
            *_obj._sys->solution = sol;


            libmesh_assert(_obj._initialized);

            // check if the system solved is linear or nonlinear
            bool if_vk = _obj._if_vk;

            ///////////////////////////////////////////////////////////////
            // first, solve the quasi-steady problem
            ///////////////////////////////////////////////////////////////
            // set the number of load steps
            unsigned int
                    n_steps = 1;
            if (if_vk) n_steps = _n_steps;

            Real
                    T0      = (*_obj._temp)(),
                    V0      = (*_obj._velocity)(),
                    p0      = (*_obj._p_cav)();


            //////////////////////////////////////////////////////////////////////
            libMesh::out << " clear the solution " << std::endl;
            // first zero the solution
            _obj._sys->solution->zero();
            libMesh::out << " solution cleared" << std::endl;
            //////////////////////////////////////////////////////////////////////
            libMesh::out << " clear stress " << std::endl;
            _obj.clear_stresss();
            libMesh::out << " stress cleared" << std::endl;


            _obj._nonlinear_assembly->set_discipline_and_system(*_obj._discipline,
                                                                *_obj._structural_sys);

            _obj._nonlinear_elem_ops->set_discipline_and_system(*_obj._discipline,
                                                                *_obj._structural_sys);


            _obj._modal_assembly->set_discipline_and_system(*_obj._discipline,
                                                                *_obj._structural_sys);

            _obj._modal_elem_ops->set_discipline_and_system(*_obj._discipline,
                                                                *_obj._structural_sys);

            bool if_continuation_solver = _obj._if_continuation_solver;

            if (!if_continuation_solver) {

                // apply scaling on thermal jaobian
                MAST::Examples::ThermalJacobianScaling
                        *jac_scaling =  dynamic_cast<MAST::Examples::ThermalJacobianScaling*>(_obj._jac_scaling);
                jac_scaling->set_assembly(*_obj._nonlinear_assembly);

                // now iterate over the load steps
                for (unsigned int i = 0; i < n_steps; i++) {

                    jac_scaling->set_enable(true);
                    _obj._nonlinear_assembly->reset_residual_norm_history();


                    // modify aero component
                    (*_obj._velocity)() = V0 * (i + 1.) / (1. * n_steps);

                    // modify the thermal load if specified by the user
                    if (!_if_only_aero_load_steps) {

                        (*_obj._temp)() = T0 * (i + 1.) / (1. * n_steps);
                        (*_obj._p_cav)() = p0;//*(i+1.)/(1.*n_steps);
                    }

                    libMesh::out
                            << "Load step: " << i
                            << "  : T = " << (*_obj._temp)()
                            << "  : p = " << (*_obj._p_cav)()
                            << "  : V = " << (*_obj._velocity)()
                            << std::endl;


                    // we can solve the problem using solve fiunctionality in sys by giving it
                    // assembly and elem_ops
                    _obj._sys->solve(*_obj._nonlinear_elem_ops, *_obj._nonlinear_assembly);



                    //  we can then write the solution into the .exo file which will contain all variables
                    writer->write_timestep("sol_n_r_solver.exo",
                                           *_obj._eq_sys,
                                           inc + 1,
                                           inc + 1);
                    inc++;

                }

                // disable scaling of the thermal jacobian cuz we need the true jacobian for sens
                jac_scaling->set_enable(false);
                jac_scaling->clear_assembly();
            }
            else {

                // Solve the system and print displacement degrees-of-freedom to screen.
                libMesh::Point
                        pt((_obj._length)/2., (_obj._width)/2., 0.), // location of mid-point before shift
                        pt0,
                        R = 0, //radius of the circle where the circumference defines the curved plate
                        dr1, dr2;
                const libMesh::Node
                *nd = nullptr;

                // if a finite radius is defined, change the mesh to a circular arc of specified radius
                libMesh::MeshBase::node_iterator
                        n_it   = _obj._mesh->nodes_begin(),
                        n_end  = _obj._mesh->nodes_end();

                // initialize the pointer to a node
                nd   = *n_it;
                pt0  = *nd;

                for ( ; n_it != n_end; n_it++) {

                    dr1  = pt0;
                    dr1 -= pt;

                    dr2  = **n_it;
                    dr2 -= pt;

                    if (dr2.norm() < dr1.norm()) {

                        nd  = *n_it;
                        pt0 = *nd;
                    }
                }

                // std::cout << *nd << std::endl;

                Real
                        max_temp = (*_obj._temp)(),
                        dt =1;

                Real init_step      = _obj._input("init_step", "init_temperature  for c-s",  (*_obj._temp)()/50);
                (*_obj._temp)() = init_step;

                const unsigned int
                        dof_num = nd->dof_number(0, 2, 0);

                unsigned int
                n_temp_steps  = _obj._input( "n_temp_steps", "number of load steps for temperature increase",  1000);
                        

                // write the header to the load.txt file

                std::ofstream out;   // text file for nl solution
                std::ofstream out_eig;  // text file for eigenvalues

                if (_obj.comm().rank() == 0) {

                    out.open("continuation_solver_load.txt", std::ofstream::out);
                    out
                            << std::setw(10) << "iter"
                            << std::setw(25) << "temperature"
                            << std::setw(25) << "pressure"
                            << std::setw(25) << "displ" << std::endl;

                    out_eig.open("continuation_solver_eig.txt", std::ofstream::out);
                    out_eig
                            << std::setw(10) << "iter"
                            << std::setw(25) << "temperature"
                            << std::setw(25) << "pressure";

                    for (int di = 0; di < _obj._n_eig; di++)
                        out_eig  << std::setw(25) << "Re_of_eigenvalue" << di+1;

                    out_eig << std::endl;
                }

                // first solve the the temperature increments
                std::vector<Real> vec1;
                std::vector<unsigned int> vec2 = {dof_num};


                if (n_temp_steps) {

                    MAST::PseudoArclengthContinuationSolver solver;
                    solver.schur_factorization = _obj._input("if_schur_factorization", "use Schur-factorization in continuation solver", true);
                    solver.min_step            = _obj._input("min_step", "minimum arc-length step-size for continuation solver",          10.);
                    solver.max_it              = _obj._input("max_it", "max nr iterations",          10);
                    solver.max_step            = _obj._input("max_step", "maximum arc-length step-size for continuation solver",   20.);
                    solver.step_desired_iters  = _obj._input("step_desired_iters", "maximum arc-length step-size for continuation solver",5);



                    // specify temperature as the load parameter to be changed per
                    // load step
                    solver.set_assembly_and_load_parameter(*_obj._nonlinear_elem_ops,
                                                           *_obj._nonlinear_assembly,
                                                           *_obj._temp);

                    // the initial deformation direction is identified with a
                    // unit change in temperature.
                    solver.initialize((*_obj._temp)());
                    // with the search direction defined, we define the arc length
                    // per load step to be a factor of 2 greater than the initial step.

                    solver.arc_length *= 2;

                    std::vector<Real> eig_vec(_obj._n_eig, 0);

                    for (unsigned int i=0; i<n_temp_steps; i++) {

                        solver.solve();
                        libMesh::out
                                << "  iter: " << i
                                << "  temperature: " << (*_obj._temp)()
                                << "  pressure: "    << (*_obj._p_cav)() << std::endl;


                        // get the value of the node at the center of the plate for output
                        _obj._sys->solution->localize(vec1, vec2);

                        // write the value to the load.txt file
                        if (_obj.comm().rank() == 0) {
                            out
                                    << std::setw(10) << i
                                    << std::setw(25) << (*_obj._temp)()
                                    << std::setw(25) << (*_obj._p_cav)()
                                    << std::setw(25) << vec1[0] << std::endl;
                        }


                        if (true //i%2== 0
                                ){
                            _obj._modal_assembly->set_base_solution(*_obj._sys->solution);
                            _obj._sys->eigenproblem_solve( *_obj._modal_elem_ops, *_obj._modal_assembly);
                            unsigned int
                                    nconv = std::min(_obj._sys->get_n_converged_eigenvalues(),
                                                     _obj._sys->get_n_requested_eigenvalues());

                            if (_obj.comm().rank() == 0) {
                                out_eig
                                        << std::setw(10) << i
                                        << std::setw(25) << (*_obj._temp)()
                                        << std::setw(25) << (*_obj._p_cav)();
                            }

                            eig_vec.resize(nconv);
                            for (int dj =0 ; dj < nconv; dj++){
                                // now write the eigenvalue
                                Real
                                        re = 0.,
                                        im = 0.;
                                _obj._sys->get_eigenvalue(dj, re, im);
                                eig_vec[dj] = re;

                                if (_obj.comm().rank() == 0) {
                                    out_eig  << std::setw(25) << re  ;
                                    }
                            }

                            if (nconv < _obj._n_eig) {
                                int diff_eigs = _obj._n_eig - nconv ;
                                for (int di = 0; di < diff_eigs; di++)
                                    out_eig << std::setw(25) << "N/A";
                            }

                            out_eig << std::endl;
                            _obj._modal_assembly->clear_base_solution();
                        }


                        _obj._sys->time += dt;
//                         write the current solution to the exodus file for
//                         visualization
                        writer->write_timestep("sol_continuation_solver.exo",
                                               *_obj._eq_sys,
                                                 i+1,
                                                 i+1);//,
                                              // _obj._sys->time);


                        // if a negative eigenvalue is detected
                        // change flag to true to increase obj and fvals
                        // and solve the system one last time and exit

                        auto min_eig  = std::min_element(eig_vec.begin(),eig_vec.end());
                        if ( (*min_eig < 100.0) && ( (*_obj._temp)() < max_temp) )  {
                            _obj._if_neg_eig = true;
                            libMesh::out << " negative eigenvalue found" << std::endl;
                            (*_obj._temp)() = max_temp;
                            break;
                        }

                        if (  (*_obj._temp)() < 0.0 )   {
                            _obj._if_neg_eig = true;
                            libMesh::out << " Continuation solver diverged" << std::endl;
                            (*_obj._temp)() = max_temp;
                            break;
                        }
                        
                        // if the temperature given by the solver is bigger than tmax
                        // go back to tmax and solve the system one more time and exit
                        if ((*_obj._temp)() > max_temp) {
                            (*_obj._temp)() = max_temp;
                            _obj._sys->solve(*_obj._nonlinear_elem_ops,
                                             *_obj._nonlinear_assembly);
                            break;
                        }
                    }
                    // clear assembly and loading parameter from continuation solver
                    solver.clear_assembly_and_load_parameters();
                }
            }


            // copy the solution to the base solution vector
            sol = *_obj._sys->solution;

            _obj._nonlinear_assembly->clear_discipline_and_system();
            _obj._nonlinear_elem_ops->clear_discipline_and_system();

            _obj._modal_assembly->clear_discipline_and_system();
            _obj._modal_elem_ops->clear_discipline_and_system();
            return sol;
        }


        /*!
         *   sets the number of steps to be used for nonlinaer steady analysis.
         *   The default is 25 for noninear and 1 for linear.
         */
        void set_n_load_steps(unsigned int n) {
            _n_steps = n;
        }


        /*!
         *   Tells the solver to increment only the aero loads during the
         *   nonlinear stepping. \p false by default.
         */
        void set_modify_only_aero_load(bool f) {
            _if_only_aero_load_steps = f;
        }


        /*!
         * @returns  a non-const-reference to the solution.
         */
        virtual libMesh::NumericVector<Real>&
        solution() { return _obj._sys->get_vector("base_solution"); }



        /*!
         * @returns  a const-reference to the solution.
         */
        virtual const libMesh::NumericVector<Real>&
        solution() const { return _obj._sys->get_vector("base_solution"); }

    };
};



//StiffenedPlateThermallyStressedPistonTheorySizingOptimization* _my_func_eval=nullptr;


int main(int argc, char* argv[]) {

    libMesh::LibMeshInit init(argc, argv);

    MAST::Examples::GetPotWrapper
            _input(argc, argv, "input");

    bool         verify_grads  = _input("verify_grads",  "verify the gradients", false);

    //verify_grads = true;

    // create and attach sizing optimization object
    StiffenedPlateThermallyStressedPistonTheorySizingOptimization func_eval(init.comm(), _input);

    if (init.comm().rank() == 0)
        func_eval.set_output_file("optimization_output.txt");

    __my_func_eval = &func_eval;


    MAST::GCMMAOptimizationInterface optimizer;
    

    unsigned int
            max_inner_iters        = _input("max_inner_iters", "maximum inner iterations in GCMMA", 15);

    Real
            constr_penalty         = _input("constraint_penalty", "constraint penalty in GCMMA", 50.),
            initial_rel_step       = _input("initial_rel_step", "initial step size in GCMMA", 1.e-2),
            asymptote_reduction    = _input("asymptote_reduction", "reduction of aymptote in GCMMA", 0.7),
            asymptote_expansion    = _input("asymptote_expansion", "expansion of asymptote in GCMMA", 1.2);

    optimizer.set_real_parameter   ("constraint_penalty",  constr_penalty);
    optimizer.set_real_parameter   ("initial_rel_step",  initial_rel_step);
    optimizer.set_real_parameter   ("asymptote_reduction",  asymptote_reduction);
    optimizer.set_real_parameter   ("asymptote_expansion",  asymptote_expansion);
    optimizer.set_integer_parameter(   "max_inner_iters", max_inner_iters);
    
    if (verify_grads) {
        std::vector<Real>
                dvals(func_eval.n_vars()),
                dummy(func_eval.n_vars());
        func_eval.init_dvar(dvals, dummy, dummy);

        libMesh::out << "******* Begin: Verifying gradients ***********" << std::endl;
        func_eval.verify_gradients(dvals);
        libMesh::out << "******* End: Verifying gradients ***********" << std::endl;
    }
    else {
        // attach and optimize
        optimizer.attach_function_evaluation_object(func_eval);
        optimizer.optimize();
    }


    //output.close();
    
    // END_TRANSLATE
    return 0;
}





//void
//stiffened_plate_thermally_stressed_piston_theory_flutter_optim_obj(int*    mode,
//                                                                   int*    n,
//                                                                   double* x,
//                                                                   double* f,
//                                                                   double* g,
//                                                                   int*    nstate) {
//
//
//    // make sure that the global variable has been setup
//    libmesh_assert(__my_func_eval);
//
//    // initialize the local variables
//    Real
//            obj = 0.;
//
//    unsigned int
//            n_vars  =  __my_func_eval->n_vars(),
//            n_con   =  __my_func_eval->n_eq()+__my_func_eval->n_ineq();
//
//    libmesh_assert_equal_to(*n, n_vars);
//
//    std::vector<Real>
//            dvars   (*n,    0.),
//            obj_grad(*n,    0.),
//            fvals   (n_con, 0.),
//            grads   (0);
//
//    std::vector<bool>
//            eval_grads(n_con);
//    std::fill(eval_grads.begin(), eval_grads.end(), false);
//
//    // copy the dvars
//    for (unsigned int i=0; i<n_vars; i++)
//        dvars[i] = x[i];
//
//
//    __my_func_eval->evaluate(dvars,
//                             obj,
//                             true,       // request the derivatives of obj
//                             obj_grad,
//                             fvals,
//                             eval_grads,
//                             grads);
//
//
//    // now copy them back as necessary
//    *f  = obj;
//    for (unsigned int i=0; i<n_vars; i++)
//        g[i] = obj_grad[i];
//}
//
//
//void
//stiffened_plate_thermally_stressed_piston_theory_flutter_optim_con(int*    mode,
//                                                                   int*    ncnln,
//                                                                   int*    n,
//                                                                   int*    ldJ,
//                                                                   int*    needc,
//                                                                   double* x,
//                                                                   double* c,
//                                                                   double* cJac,
//                                                                   int*    nstate) {
//
//
//    // make sure that the global variable has been setup
//    libmesh_assert(__my_func_eval);
//
//    // initialize the local variables
//    Real
//            obj = 0.;
//
//    unsigned int
//            n_vars  =  __my_func_eval->n_vars(),
//            n_con   =  __my_func_eval->n_eq()+__my_func_eval->n_ineq();
//
//    libmesh_assert_equal_to(    *n, n_vars);
//    libmesh_assert_equal_to(*ncnln, n_con);
//
//    std::vector<Real>
//            dvars   (*n,    0.),
//            obj_grad(*n,    0.),
//            fvals   (n_con, 0.),
//            grads   (n_vars*n_con, 0.);
//
//    std::vector<bool>
//            eval_grads(n_con);
//    std::fill(eval_grads.begin(), eval_grads.end(), true);
//
//    // copy the dvars
//    for (unsigned int i=0; i<n_vars; i++)
//        dvars[i] = x[i];
//
//
//    __my_func_eval->evaluate(dvars,
//                             obj,
//                             true,       // request the derivatives of obj
//                             obj_grad,
//                             fvals,
//                             eval_grads,
//                             grads);
//
//
//    // now copy them back as necessary
//
//    // first the constraint functions
//    for (unsigned int i=0; i<n_con; i++)
//        c[i] = fvals[i];
//
//    // next, the constraint gradients
//    for (unsigned int i=0; i<n_con*n_vars; i++)
//        cJac[i] = grads[i];
//
//
//}




