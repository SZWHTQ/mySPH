cmake_minimum_required(VERSION 3.11)


add_executable(check "check.f90")
add_test("${PROJECT_NAME}/chech" check)


add_executable(toml-fTest "toml-fTest.f90")
if(Toml-f_FOUND)
    target_include_directories(
        toml-fTest
        PRIVATE
        ${Toml-f_INCLUDE_DIR}
    )
    target_link_libraries(
        toml-fTest
        PRIVATE
        ${Toml-f_LIBRARY}
        mySPH-lib
    )
else(Toml-f_FOUND)
    target_link_libraries(
        toml-fTest
        PRIVATE
        mySPH-lib
        toml-f
    )
endif(Toml-f_FOUND)
set_target_properties(
    toml-fTest
    PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/module"
)
add_test("${PROJECT_NAME}/toml-f" toml-fTest)


add_executable(progressBarTest "progressBarTest.f90")
target_link_libraries(
    progressBarTest
    PRIVATE
    mySPH-lib
)
set_target_properties(
    progressBarTest
    PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/module"
)
add_test("${PROJECT_NAME}/progressBar" progressBarTest)


add_executable(toStringTest "toStringTest.f90")
target_link_libraries(
    toStringTest
    PRIVATE
    mySPH-lib
)
set_target_properties(
    toStringTest
    PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/module"
)
add_test("${PROJECT_NAME}/toString" toStringTest)


add_executable(roundTest "roundTest.f90")
target_link_libraries(
    roundTest
    PRIVATE
    mySPH-lib
)
set_target_properties(
    roundTest
    PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/module"
)
add_test("${PROJECT_NAME}/round" roundTest)


add_executable(exceptionTest "exceptionTest.f90")
target_link_libraries(
    exceptionTest
    PRIVATE
    mySPH-lib
)
set_target_properties(
    exceptionTest
    PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/module"
)
add_test("${PROJECT_NAME}/exception" exceptionTest)

add_executable(kernelTest "kernelTest.f90")
target_link_libraries(
    kernelTest
    PRIVATE
    mySPH-lib
)
set_target_properties(
    kernelTest
    PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/module"
)
add_test("${PROJECT_NAME}/kernel" kernelTest)
