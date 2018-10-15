message("this is to test the systems random number generators")

add_executable(test_std test_std.cpp)

find_package(GSL REQUIRED)
add_executable(test_gsl test_gsl.cpp)
target_link_libraries(test_gsl ${GSL_LIBRARIES})

find_package(Boost REQUIRED COMPONENTS "random")
add_executable(test_boost test_boost.cpp)
message("Boost_LIBRARIES: ${Boost_LIBRARIES}")
target_link_libraries(test_boost ${Boost_LIBRARIES})
