# Sequential root CMakeLists
add_subdirectory(directMethodSSA)
add_subdirectory(boundedTauLeaping)
add_subdirectory(common)

add_executable(stochemsim_direct driverSSA.cpp)
add_executable(stochemsim_btl driverBTL.cpp)

target_include_directories(stochemsim_direct PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(stochemsim_btl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link libraries necessary for compiling BTL and direct method
target_link_libraries(stochemsim_direct directMethod common)
target_link_libraries(stochemsim_btl btl common)