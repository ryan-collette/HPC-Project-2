
add_library(simulator psystem.c)
target_link_libraries(simulator PRIVATE MPI::MPI_C)

add_executable(nbody nbody.c)
target_link_libraries(nbody simulator)
