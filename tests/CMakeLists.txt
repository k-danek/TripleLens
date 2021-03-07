set(BINARY triplelens_test)

set(TEST_SOURCES "../tests/main.cc")

list(APPEND TEST_SOURCES "../tests/laguerretest.cc")

set(SOURCES ${TEST_SOURCES})

add_executable(${BINARY} ${TEST_SOURCES})

add_test(NAME ${BINARY} COMMAND ${BINARY})

target_link_libraries(${BINARY} PUBLIC gtest ccc)

