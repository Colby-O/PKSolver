message("-- PARAMETER_EXTRACTION - Loading Eigen3...")

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIRS})

if(NOT EIGEN3_FOUND)
	message(ERROR " EIGEN not found!")
else()
	add_definitions(-DWITH_EIGEN)
endif() 
