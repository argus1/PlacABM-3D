add_library(PhysiMeSS OBJECT)

file(GLOB SOURCES "*.cpp")
target_sources(PhysiMeSS PRIVATE ${SOURCES})

find_package(OpenMP 4)
if(OpenMP_CXX_FOUND)
  target_link_libraries(PhysiMeSS PUBLIC OpenMP::OpenMP_CXX)
endif()
