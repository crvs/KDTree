
add_executable(construction_time construction_time.cpp)
target_link_libraries(construction_time KDTree)

add_executable(error_test error_test.cpp)
target_link_libraries(error_test KDTree)

add_executable(toy_test toy_test.cpp)
target_link_libraries(toy_test KDTree)
