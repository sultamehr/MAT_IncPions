add_library(util SafeROOTName.cpp)
target_link_libraries(util ${ROOT_LIBRARIES})
install(TARGETS util DESTINATION lib)
