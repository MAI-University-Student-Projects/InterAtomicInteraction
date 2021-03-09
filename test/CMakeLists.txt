find_package(Catch2 REQUIRED) # catch2 as system package requirement
set(_TEST_FILES test-atomic.cpp)

# add executable app from source files
add_executable(inter_atomic_test test-main.cpp ${_TEST_FILES})

target_include_directories(inter_atomic_test PRIVATE ../src)
# link target (executable) with catch2 by private flag
# target_link_libraries(inter_atomic_test PRIVATE Catch2::Catch2WithMain)

target_link_libraries(inter_atomic_test PRIVATE Catch2::Catch2WithMain atomic_lib)

#create test named 'atomic_ctest' invokable through ctest
add_test(NAME atomic_ctest COMMAND inter_atomic_test)
