add_executable(TwoBlocksColliding two_block_colliding.cpp)
target_link_libraries(TwoBlocksColliding LINK_PUBLIC Cabana::cabanacore Cabana::Cajita)

add_executable(DB2d db_2d.cpp)
target_link_libraries(DB2d LINK_PUBLIC Cabana::cabanacore Cabana::Cajita)

add_executable(DB2dBrute db_2d_without_nnps.cpp)
target_link_libraries(DB2dBrute LINK_PUBLIC Cabana::cabanacore Cabana::Cajita)

add_executable(01RBFreelyRotatingBodies 01_rb_freely_rotating_square_bodies.cpp)
target_link_libraries(01RBFreelyRotatingBodies LINK_PUBLIC Cabana::cabanacore Cabana::Cajita)

install(TARGETS 01RBFreelyRotatingBodies TwoBlocksColliding DB2d DB2dBrute DESTINATION ${CMAKE_INSTALL_BINDIR})

add_definitions(-w)
