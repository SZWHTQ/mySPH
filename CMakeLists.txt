cmake_minimum_required(VERSION 3.11)

# set(CMAKE_SYSTEM_NAME Linux)

# set(CMAKE_BUILD_TYPE Debug)

project(
    "mySPH"
    LANGUAGES "Fortran"
    VERSION "0.2.0"
)
# set(CMAKE_GENERATOR "Ninja")

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    add_definitions("-std=f2018 -cpp")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    if("{CMAKE_SYSTEM_NAME}" MATCHES "Windows")
        add_definitions("/fpp /stand=f18")
    else()
        add_definitions("-fpp -stand=f18")
    endif()
endif()

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

set(OMP ON CACHE BOOL "Enable OpenMP")
if(OMP)
    if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lgomp -fopenmp")
    elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
        if(CMAKE_SYSTEM_NAME MATCHES "Windows")
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /Qopenmp")
        else()
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
        endif()
    endif()
    message(STATUS "OpenMP enabled")
else(OMP)
    message(STATUS "OpenMP disabled")
endif()

# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -fbacktrace")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

add_library(
    ${PROJECT_NAME}-lib
    ${srcs}
)

find_package(Toml-f)

if(LOCAL_TOMLF)
    if(Toml-f_FOUND)
        message(STATUS "Local toml-f library found")
        target_include_directories(
            ${PROJECT_NAME}-lib
            PRIVATE
            ${Toml-f_INCLUDE_DIR}
        )
        target_link_libraries(
            ${PROJECT_NAME}-lib
            ${Toml-f_LIBRARY}
        )
    elseif(NOT Toml-f_FOUND)
        message(STATUS "Local toml-f library not found, downloading and building it")
        include(toml-f)
        target_link_libraries(
            ${PROJECT_NAME}-lib
            toml-f
        )
    endif(Toml-f_FOUND)
elseif(NOT LOCAL_TOMLF)
    message(STATUS "Downloading and building toml-f")
    include(toml-f)
    target_link_libraries(
        ${PROJECT_NAME}-lib
        toml-f
    )
endif(LOCAL_TOMLF)

# target_link_libraries(
#     ${PROJECT_NAME}-lib
#     blas
#     lapack
# )

add_executable(
    ${PROJECT_NAME}
    "src/main.f90"
)

target_link_libraries(
    ${PROJECT_NAME}
    ${PROJECT_NAME}-lib
)

set_target_properties(
    ${PROJECT_NAME}-lib
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

# enable_testing()
# add_subdirectory(test)
