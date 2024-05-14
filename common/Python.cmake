message("-- PARAMETER_EXTRACTION - Loading Python...")

find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)

include_directories(${PYTHON_INCLUDE_DIRS})

#--- Error Handling
if(NOT PYTHONINTERP_FOUND)
	message(ERROR " Python Interpreter not found!")
else()
	add_definitions(-DWITH_PYTHONINTERP)
endif()

if(NOT PYTHONLIBS_FOUND)
	message(ERROR " Python libs not found!")
else()
	add_definitions(-DWITH_PYTHONLIBS)
endif()
