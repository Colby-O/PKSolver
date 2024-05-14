message(STATUS "PARAMETER_EXTRACTION - Loading ITK...")

find_package(ITK REQUIRED)

#--- Error Handling
if(NOT ITK_FOUND)
	messgae(ERROR "ITK was not found!")
else()
	include(${ITK_USE_FILE})
endif()
