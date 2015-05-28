project(EOA)  # Elimination Ordering Agrum
cmake_minimum_required(VERSION 2.8)

# do not forget to change this line if needed ("act install -d...")
set(AGRUM_INSTALLATION_DIRECTORY "/usr")
set(aGrUM_DIR "${AGRUM_INSTALLATION_DIRECTORY}/lib/cmake/aGrUM/")

find_package(aGrUM)

if (aGrUM_FOUND)
    include(${AGRUM_USE_FILE})
  # include_directories(${AGRUM_INCLUDE_DIR})
  # link_directories(${AGRUM_LIB_DIR})
else (aGrUM_FOUND)
  message(FATAL_ERROR "Please install aGrUM")
endif (aGrUM_FOUND)

# cmake -DCMAKE_BUILD_TYPE=DEBUG
# or
# cmake -DCMAKE_BUILD_TYPE=RELEASE
#     RELEASE is the default option (thanks to the next 3 lines)
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release)
endif()

file(GLOB_RECURSE EOA_SOURCE ${EOA_SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE EOA_INCLUDE ${EOA_SOURCE_DIR}/*.h)

add_executable (eoa ${EOA_SOURCE})

target_link_libraries(eoa ${AGRUM_LIBRARIES})
target_link_libraries(eoa ${LIBAGRUM} "stdc++")
