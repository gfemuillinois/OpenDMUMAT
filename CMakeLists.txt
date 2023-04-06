cmake_minimum_required(VERSION 3.20)

project(OpenDMUmat)

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


if (${OpenDM_Exec})
  add_executable(openDMUmat ${OpenDMUmat_SRCS})
endif()

if (${OpenDM_Lib})
  add_library(openDMUmat SHARED ${OpenDMUmat_SRCS})
endif()



target_link_libraries(openDMUmat openDMLib)
