# File to make library out of files here

# Gather source and header files
file(GLOB_RECURSE OpenDMLib_SRCS CONFIGURE_DEPENDS "./*cpp")
file(GLOB_RECURSE OpenDMLib_HEADER CONFIGURE_DEPENDS "./*hpp")

set(OpenDMUmat_INCLUDE_DIRS "")
foreach(_headerFile ${foo_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list(APPEND OpenDMLib_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES OpenDMLib_INCLUDE_DIRS)
# Add eigen headers
find_path(EIGEN_HEADERS NAMES /Eigen/Dense PATHS /Users/bryce/eigen/eigen-3.4.0 "C:/Users/Bryce/eigen-3.4.0/eigen-3.4.0") 
list(APPEND OpenDMLib_INCLUDE_DIRS ${EIGEN_HEADERS})
list(APPEND OpenDMLib_HEADER  "${EIGEN_HEADERS}/Eigen/Dense")

# Print things
message("OpenDMLib_SRCS=${OpenDMLib_SRCS}")
message("OpenDMLib_HEADER=${OpenDMLib_HEADER}")
message("OpenDMLib_INCLUDE_DIRS=${OpenDMLib_INCLUDE_DIRS}")
message("EIGEN_HEADERS=${EIGEN_HEADERS}")

# Make lib
add_library(openDMLib STATIC ${OpenDMLib_SRCS})
target_include_directories(openDMLib PRIVATE ${OpenDMLib_INCLUDE_DIRS})
