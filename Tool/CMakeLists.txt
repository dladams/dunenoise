# myproj/Tool/CMakeLists.txt

install(FILES ApaWirePlot.h TYPE INCLUDE)

# Build a tool.
add_library(ApaWirePlot_tool SHARED ApaWirePlot_tool.cc)
find_library(libDuneCommon NAMES dune_DuneCommon)
find_library(libDuneGeometry NAMES dune_Geometry)
target_include_directories(ApaWirePlot_tool BEFORE PRIVATE ../Interface)
target_link_libraries(ApaWirePlot_tool ${tool_lib_list})
target_link_libraries(ApaWirePlot_tool ${libDuneCommon})
target_link_libraries(ApaWirePlot_tool ${libDuneGeometry})
target_link_libraries(ApaWirePlot_tool ${ROOT_LIBRARIES})
install(TARGETS ApaWirePlot_tool DESTINATION lib)

#add_subdirectory(test)
