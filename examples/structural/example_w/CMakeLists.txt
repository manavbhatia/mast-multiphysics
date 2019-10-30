# since an optimizer is required for optimization we turn this on only if
# NLOPT was enabled for MAST
if(ENABLE_NLOPT OR ENABLE_SNOPT OR ENABLE_GCMMA)
    add_executable(structural_example_w
                    example_w.cpp
                    ${MAST_ROOT_DIR}/examples/fluid/meshing/cylinder.cpp
                    ${MAST_ROOT_DIR}/examples/old/structural/stiffened_plate_optimization/stiffened_plate_optimization_base.cpp
                    ${MAST_ROOT_DIR}/examples/base/multilinear_interpolation.cpp
                    ${MAST_ROOT_DIR}/examples/old/structural/stiffened_plate_optimization/stiffened_plate_optimization_base.cpp
                    ${MAST_ROOT_DIR}/examples/structural/base/thermal_stress_jacobian_scaling_function.cpp
                    ${MAST_ROOT_DIR}/examples/fluid/meshing/cylinder.cpp
                    ${MAST_ROOT_DIR}/examples/base/multilinear_interpolation.cpp
                    ${MAST_ROOT_DIR}/examples/structural/base/blade_stiffened_panel_mesh.cpp
                    ${MAST_ROOT_DIR}/examples/base/plot_results.cpp
                    ${MAST_ROOT_DIR}/examples/old/structural/plate_optimization/plate_optimization_base.cpp
                     )

    target_include_directories(structural_example_w PRIVATE
                                ${CMAKE_CURRENT_LIST_DIR}
                                ${PROJECT_SOURCE_DIR}
                                ${PROJECT_SOURCE_DIR}/src)
    target_link_libraries(structural_example_w mast)

    install(TARGETS structural_example_w
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
endif()
