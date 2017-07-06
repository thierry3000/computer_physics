#find_package(Boost REQUIRED COMPONENTS "program_options")
#find_package(Boost REQUIRED COMPONENTS "filesystem")
#find_package(Boost REQUIRED COMPONENTS "system")

#not working on cip pool
#find_package(Boost REQUIRED COMPONENTS "system")

add_executable(Lanczos Lanczos.cpp)

#not working on cip pool
#target_link_libraries(bosehubbard gsl gslcblas boost_program_options boost_filesystem boost_system)


#target_link_libraries(bosehubbard boost_program_options boost_filesystem boost_system)
