cmake_minimum_required(VERSION 3.11)

# set(CMAKE_SYSTEM_NAME Generic)

# set(CMAKE_BUILD_TYPE Debug)

project(
    "mySPH"
    LANGUAGES "Fortran"
    VERSION "0.1.0"    
)

# file(GLOB_RECURSE srcs0 src/0/*.f90)
# file(GLOB_RECURSE srcs1 src/1/*.f90)
# file(GLOB_RECURSE srcs2 src/2/*.f90)
# file(GLOB_RECURSE srcs3 src/3/*.f90)
set(srcs)
add_subdirectory("src")

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/module")
    make_directory("${PROJECT_SOURCE_DIR}/module")
endif()

set(PROF OFF CACHE BOOL "Enable gprof profiling")
if(PROF)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pg") # Compiler ptions for gprof profiling
    message(STATUS "Gprof profiling enabled")
endif()

set(OMP OFF CACHE BOOL "Enable OpenMP")
if(OMP)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    message(STATUS "OpenMP enabled")
else(OMP)
    message(STATUS "OpenMP disabled")
endif()

# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -fbacktrace")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

add_library(
    mySPH-lib
    ${srcs}
)

find_package(Toml-f)

if(LOCAL_TOMLF)
    if(Toml-f_FOUND)
        message(STATUS "Local toml-f library found")
        target_include_directories(
            mySPH-lib
            PRIVATE
            ${Toml-f_INCLUDE_DIR}
        )
        target_link_libraries(
            mySPH-lib
            ${Toml-f_LIBRARY}
        )
    elseif(NOT Toml-f_FOUND)
        message(STATUS "Local toml-f library not found, downloading and building it")
        include(toml-f)
        target_link_libraries(
            mySPH-lib
            toml-f
        )
    endif(Toml-f_FOUND)
elseif(NOT LOCAL_TOMLF)
    message(STATUS "Downloading and building toml-f")
    include(toml-f)
    target_link_libraries(
        mySPH-lib
        toml-f
    )
endif(LOCAL_TOMLF)


add_executable(
    mySPH
    src/main.f90
)

target_link_libraries(
    ${PROJECT_NAME}
    mySPH-lib
)

set_target_properties(
    mySPH-lib
    PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    OUTPUT_NAME "${PROJECT_NAME}"
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION_MAJOR}"
    Fortran_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/module"
)

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION_MAJOR}"
    Fortran_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/module"
)

enable_testing()
add_subdirectory(test)
