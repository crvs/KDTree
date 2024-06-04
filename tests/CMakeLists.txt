add_executable(construction_time construction_time.cpp)
target_link_libraries(construction_time KDTree)

add_executable(error_test error_test.cpp)
target_link_libraries(error_test KDTree)

add_executable(knn_error_test knn_error_test.cpp)
target_link_libraries(knn_error_test KDTree)

add_executable(toy_test toy_test.cpp)
target_link_libraries(toy_test KDTree)

add_test(NAME construction_time COMMAND $<TARGET_FILE:construction_time>)
add_test(NAME error_test COMMAND $<TARGET_FILE:error_test>)
add_test(NAME knn_error_test COMMAND $<TARGET_FILE:knn_error_test>)
add_test(NAME toy_test COMMAND $<TARGET_FILE:toy_test>)

add_custom_target(check ${CMAKE_CTEST_COMMAND} DEPENDS construction_time error_test knn_error_test toy_test)
