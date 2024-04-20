add_executable(TwoBlocksColliding two_block_colliding.cpp)
target_link_libraries(TwoBlocksColliding LINK_PUBLIC Cabana::cabanacore Cabana::Cajita)

add_executable(DB2d db_2d.cpp)
target_link_libraries(DB2d LINK_PUBLIC Cabana::cabanacore Cabana::Cajita)

add_executable(DB2dBrute db_2d_without_nnps.cpp)
target_link_libraries(DB2dBrute LINK_PUBLIC Cabana::cabanacore Cabana::Cajita)

install(TARGETS TwoBlocksColliding DB2d DB2dBrute DESTINATION ${CMAKE_INSTALL_BINDIR})
