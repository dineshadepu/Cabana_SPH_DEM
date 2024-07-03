add_executable(01DEM 01_dem_two_head_on_collision.cpp)
target_link_libraries(01DEM LINK_PUBLIC Cabana::cabanacore Cabana::Cajita)

add_executable(02DEMManyParticles 02_dem_many_particles_collision.cpp)
target_link_libraries(02DEMManyParticles LINK_PUBLIC Cabana::cabanacore Cabana::Cajita)

install(TARGETS 02DEMManyParticles 01DEM DESTINATION ${CMAKE_INSTALL_BINDIR})



add_definitions(-w)
