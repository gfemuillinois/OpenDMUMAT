cmake_minimum_required(VERSION 3.16)

project(OpenDMUmat)

# Get OS and define macro for testing purposes
if (CMAKE_SYSTEM_NAME MATCHES "Windows")
  add_compile_definitions(WINDOWS_OS=1)
else()
  add_compile_definitions(WINDOWS_OS=0)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()

# add files to project
add_subdirectory(OpenDMLib)

# Get cpp file in Source Dir
list(APPEND OpenDMUmat_SRCS "./umat_opendm.cpp")

# Make an executable
option(OpenDM_Exec "Make an executable out of OpenDM" OFF)
# Make a lib
option(OpenDM_Lib "Make a shared lib out of OpenDM" ON)
if( OpenDM_Exec AND OpenDM_Lib)
  message(FATAL_ERROR "Please Make only an executable or shared lib!!!")
endif()

# Option to use numerical tangent or analytical
option(OpenDM_USE_NUM_TANG "Use numerical tangent stiffness" ON)

if (${OpenDM_Exec})
  list(APPEND OpenDMUmat_SRCS "./test_opendm.cpp")
  add_executable(openDMUmat ${OpenDMUmat_SRCS})
endif()

if (${OpenDM_Lib})
  add_library(openDMUmat SHARED ${OpenDMUmat_SRCS})
endif()



target_link_libraries(openDMUmat openDMLib)
