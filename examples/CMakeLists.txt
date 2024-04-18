add_executable(DB2D db_2d.cpp)
target_link_libraries(DB2d LINK_PUBLIC CabanaPD)

# add_executable(KalthoffWinkler kalthoff_winkler.cpp)
# target_link_libraries(KalthoffWinkler LINK_PUBLIC CabanaPD)

# add_executable(CrackBranching crack_branching.cpp)
# target_link_libraries(CrackBranching LINK_PUBLIC CabanaPD)

install(TARGETS DB2d DESTINATION ${CMAKE_INSTALL_BINDIR})
# install(TARGETS ElasticWave KalthoffWinkler CrackBranching DESTINATION ${CMAKE_INSTALL_BINDIR})
