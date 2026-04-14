add_executable(distool_tests
    main.cpp
    test_statistics.cpp
    test_descriptors.cpp
    test_classifier.cpp
    test_reader.cpp
)

target_include_directories(distool_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}   # test_harness.h
)

target_link_libraries(distool_tests PRIVATE distool_lib)

# Inject path to the test data (dump files live in the source root).
target_compile_definitions(distool_tests PRIVATE
    TEST_DATA_DIR="${CMAKE_SOURCE_DIR}"
)

add_test(NAME distool_tests COMMAND distool_tests)
