cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

############ Make your changes here ############

set(name   PriorFilter)
set(WITH_VIS true)
############ What's below should not need to be changed ############

# http://www.cmake.org/cmake/help/cmake_tutorial.html
# http://www.cmake.org/cmake/help/cmake2.6docs.html
# http://www.cmake.org/Wiki/CMake_FAQ

set(mymain ${name}.cc )
set(myexe  ${name}    )
set(myproj ${name}_   )  # any name as long as it is different from myexe

project(${myproj})
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include)
        EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(${myexe} ${mymain} ${sources} ${headers})
EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LD_FLAGS})
target_link_libraries(${myexe} ${ROOT_LIBRARIES})

install(TARGETS ${myexe} DESTINATION bin)
