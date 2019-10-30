# Install script for directory: /Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/libmast.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmast.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmast.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libmast.dylib")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/aeroelasticity/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/base/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/boundary_condition/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/coordinates/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/elasticity/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/fluid/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/heat_conduction/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/level_set/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/mesh/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/numerics/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/optimization/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/property_cards/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/solver/cmake_install.cmake")
  include("/Users/walidarsalane/Documents/computation/MAST_multiphysics_fork/mast-multiphysics/cmake-build-release/src/utility/cmake_install.cmake")

endif()

