
add_executable(nbody nbody.c)
target_link_libraries(nbody PRIVATE MPI::MPI_C)
