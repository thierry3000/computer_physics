
add_executable(ising IS.cpp )
target_link_libraries(ising boost_program_options)
add_executable(LJ LJ_basic.cpp)
target_link_libraries(LJ boost_program_options)
