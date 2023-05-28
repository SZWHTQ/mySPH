include(FetchContent)

if(NOT EXISTS "${CMAKE_BINARY_DIR}/dependencies/toml-f")
    make_directory("${CMAKE_BINARY_DIR}/dependencies/toml-f")
endif()
# set(ENV{http_proxy} "http://172.30.176.1:7890")
# set(ENV{https_proxy} "http://172.30.176.1:7890")
FetchContent_Declare(
    toml-f
    GIT_REPOSITORY https://github.com/toml-f/toml-f.git
    GIT_TAG v0.4.1
    # GIT_REPOSITORY https://gitee.com/fortran-base/toml-f.git
    # GIT_TAG v0.3.0
    SOURCE_DIR "${CMAKE_BINARY_DIR}/dependencies/toml-f"
)
FetchContent_MakeAvailable(
    toml-f
)