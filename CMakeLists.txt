cmake_minimum_required (VERSION 3.9)

project(ProgettoPCS2025 LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

# Add dependencies
################################################################################

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads
find_package(Threads REQUIRED)
list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest
find_package(GTest REQUIRED)
list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PRIVATE GTest::gtest)
list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PRIVATE GTest::gmock)
list(APPEND ${PROJECT_NAME}_LINKED_LIBRARIES PRIVATE GTest::gmock_main)

# Insert Sources
################################################################################

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src_test)
include_directories(${CMAKE_SOURCE_DIR}/ExportParaview)

add_subdirectory(src)

list(APPEND ${PROJECT_NAME}_sources ${src_sources})
list(APPEND ${PROJECT_NAME}_headers ${src_headers})
list(APPEND ${PROJECT_NAME}_includes ${src_includes})

add_subdirectory(ExportParaview)

list(APPEND ${PROJECT_NAME}_sources ${src_paraview_sources})
list(APPEND ${PROJECT_NAME}_headers ${src_paraview_headers})
list(APPEND ${PROJECT_NAME}_includes ${src_paraview_includes})

add_subdirectory(src_test)

list(APPEND ${PROJECT_NAME}_test_sources ${${PROJECT_NAME}_sources} ${src_test_sources})
list(APPEND ${PROJECT_NAME}_test_headers ${${PROJECT_NAME}_headers} ${src_test_headers})
list(APPEND ${PROJECT_NAME}_test_includes ${${PROJECT_NAME}_includes} ${src_test_includes})

# Copy file
################################################################################

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/tetrahedron/Cell0Ds.txt  ${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/tetrahedron/Cell1Ds.txt  ${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/tetrahedron/Cell2Ds.txt ${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/tetrahedron/Cell3Ds.txt 
	${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/octahedron/Cell0Ds.txt  ${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/octahedron/Cell1Ds.txt  ${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/octahedron/Cell2Ds.txt ${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/octahedron/Cell3Ds.txt 
	${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/icosahedron/Cell0Ds.txt  ${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/icosahedron/Cell1Ds.txt  ${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/icosahedron/Cell2Ds.txt ${CMAKE_CURRENT_SOURCE_DIR}/PlatonicSolid/icosahedron/Cell3Ds.txt
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Create executable
################################################################################

add_executable(${PROJECT_NAME}
      main.cpp
      ${${PROJECT_NAME}_sources}
      ${${PROJECT_NAME}_headers})

target_link_libraries(${PROJECT_NAME} ${${PROJECT_NAME}_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${${PROJECT_NAME}_includes})
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

add_executable(${PROJECT_NAME}_TEST
      main_test.cpp
	  ${${PROJECT_NAME}_test_sources}
	  ${${PROJECT_NAME}_test_headers})

target_link_libraries(${PROJECT_NAME}_TEST ${${PROJECT_NAME}_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME}_TEST PRIVATE ${${PROJECT_NAME}_test_includes})
target_compile_options(${PROJECT_NAME}_TEST PUBLIC -fPIC)
