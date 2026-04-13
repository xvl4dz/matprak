cmake_minimum_required(VERSION 3.11)
project(torus_steps)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel IOXML)

if(VTK_FOUND)
    if(${VTK_VERSION} VERSION_LESS "8.90")
        include(${VTK_USE_FILE})
    endif()
endif()

set(STEPS 
    step1_simple_torus
    step2_static_heat
    step3_moving 
    step4_strobe 
    step5_combined 
    step6_quad_static 
    step7_quad_dynamic
)

foreach(STEP ${STEPS})
    add_executable(${STEP} ${STEP}.cpp)
    target_link_libraries(${STEP} dolfin ${VTK_LIBRARIES})
endforeach()
