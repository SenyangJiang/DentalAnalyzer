# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( compare )


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()
set(CGAL_libs CGAL::CGAL)

# Boost and its components
find_package(Boost REQUIRED COMPONENTS
    thread
    filesystem
    system
    program_options
    date_time
    chrono)

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  

endif()

# Eigen 3.1 or greater
find_package( Eigen3 3.1.0 )
include(CGAL_Eigen_support)
if ( TARGET CGAL::Eigen_support )

  set(CGAL_libs ${CGAL_libs} CGAL::Eigen_support)

else()
  
  message(STATUS "This project requires Eigen 3.1 (or greater), and will not be compiled.")
  return()
  
endif()

# OpenGR
find_package(OpenGR QUIET)
include(CGAL_OpenGR_support)
if (TARGET CGAL::OpenGR_support)
  
  set(CGAL_libs ${CGAL_libs} CGAL::OpenGR_support)
  
else()
  
  message(STATUS "This project requires OpenGR, and will not be compiled.")
  
endif()

# libpointmatcher
find_package(libpointmatcher QUIET)
include(CGAL_pointmatcher_support)
if (TARGET CGAL::pointmatcher_support)

  set(CGAL_libs ${CGAL_libs} CGAL::pointmatcher_support)

else()

  message(STATUS "This project requires libpointmatcher, and will not be compiled.")

endif()

# Creating entries for target: compare
# ############################

add_executable( compare  compare.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS compare )

# Link the executable to CGAL and third-party libraries
target_link_libraries(compare PRIVATE ${CGAL_libs})

