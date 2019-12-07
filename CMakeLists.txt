cmake_minimum_required(VERSION 3.1)
project(final_proj)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "-Wno-format -Wno-deprecated-declarations -Wno-unused-result -O3")

include_directories(dependencies/ dependencies/eigen/ source/ )

add_executable(final_proj 
	main.cpp 
	source/VCFrontEnd.cpp
	source/vertexClustering.cpp
	source/meshSimplify.cpp
	source/config.cpp
	source/edge.cpp
	source/edgeHeap.cpp
	source/matrix.cpp
	source/solve.cpp
	source/vector3.cpp
	source/vector4.cpp
	source/vertex.cpp
	source/vertexGroup.cpp
)

target_link_libraries(final_proj -pthread)

