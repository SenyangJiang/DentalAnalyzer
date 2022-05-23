cmake_minimum_required(VERSION 3.5)

project(DentalAnalyzer LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CGAL_DIR "C:/dev/CGAL-5.1.1")
set(GMP_LIBRARIES "C:/dev/CGAL-5.1.1/auxiliary/gmp/lib/libgmp-10.lib")
set(GMP_INCLUDE_DIR "C:/dev/CGAL-5.1.1/auxiliary/gmp/include")
set(MPFR_LIBRARIES "C:/dev/CGAL-5.1.1/auxiliary/gmp/lib/libmpfr-4.lib")
set(MPFR_INCLUDE_DIR "C:/dev/CGAL-5.1.1/auxiliary/gmp/include")
set(Boost_INCLUDE_DIR "C:/dev/boost_1_75_0")
set(Boost_DIR "C:/dev/boost_1_75_0/stage/lib/cmake/Boost-1.75.0")
set(Boost_DEBUG OFF)
set(CGAL_Boost_USE_STATIC_LIBS ON)

set(EIGEN3_INCLUDE_DIR "C:/eigen-3.3.7")
set(OpenGR_DIR "C:/Program Files (x86)/OpenGR/lib/cmake/opengr")
set(libnabo_DIR "C:/dev/libnabo/build")
set(libpointmatcher_DIR "C:/dev/libpointmatcher/build/install/share/libpointmatcher/cmake")
set(CMAKE_BUILD_TYPE "Release")
find_package(QT NAMES Qt6 Qt5 COMPONENTS Widgets Multimedia REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets Multimedia REQUIRED)

# CGAL and its components
find_package( CGAL COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()
set(CGAL_libs CGAL::CGAL)
#
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
find_package(OpenGR)
include(CGAL_OpenGR_support)
if (TARGET CGAL::OpenGR_support)

  set(CGAL_libs ${CGAL_libs} CGAL::OpenGR_support)

else()

  message(STATUS "This project requires OpenGR, and will not be compiled.")
  return()

endif()

# libpointmatcher
find_package(libpointmatcher)
include(CGAL_pointmatcher_support)
if (TARGET CGAL::pointmatcher_support)

  set(CGAL_libs ${CGAL_libs} CGAL::pointmatcher_support)

else()

  message(STATUS "This project requires libpointmatcher, and will not be compiled.")
  return()

endif()

set(PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        mainwindow.h
        mainwindow.ui
        Analyzer.cpp
        Analyzer.h
        objects.h
        utility.cpp
        utility.h
        extract_Vpropres.cpp
        extract_Vpropres.h
        normal_cycle.h
        polyhedron_enrichment_define.h
        Roughness.h
        Roughness_Items.h
        Polyhedron/polyhedron.h
        Polyhedron/Correct_CGAL_Structure.h
        Polyhedron/Polyhedron_Builder.h
        Polyhedron/Polyhedron_Copy.h
        Polyhedron/polyhedron_enriched_polyhedron.h
        Polyhedron/Polyhedron_OBJ_CGALImporter.h
        Polyhedron/Polyhedron_OFF_CGALImporter.h
        Polyhedron/polyhedron_shared_items.h
        audiolevel.cpp
        audiolevel.h
        audiorecorder.cpp
        audiorecorder.h
        audiorecorder.ui
        result.h
        parameter.h
        result.cpp
)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(DentalAnalyzer
        ${PROJECT_SOURCES}
    )
else()
    if(ANDROID)
        add_library(DentalAnalyzer SHARED
            ${PROJECT_SOURCES}
        )
    else()
        add_executable(DentalAnalyzer
            ${PROJECT_SOURCES}
        )
    endif()
endif()

target_link_libraries(DentalAnalyzer PRIVATE Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::Multimedia ${CGAL_libs})
