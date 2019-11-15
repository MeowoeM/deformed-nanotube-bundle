cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# this is where we will place the Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# compiler flags
SET( FOX_FLAGS "-I/Users/Meow/local/FoX-4.1.2/objs/finclude /Users/Meow/local/FoX-4.1.2/objs/lib/libFoX_wkml.a /Users/Meow/local/FoX-4.1.2/objs/lib/libFoX_dom.a /Users/Meow/local/FoX-4.1.2/objs/lib/libFoX_sax.a /Users/Meow/local/FoX-4.1.2/objs/lib/libFoX_wcml.a /Users/Meow/local/FoX-4.1.2/objs/lib/libFoX_wxml.a /Users/Meow/local/FoX-4.1.2/objs/lib/libFoX_common.a /Users/Meow/local/FoX-4.1.2/objs/lib/libFoX_utils.a /Users/Meow/local/FoX-4.1.2/objs/lib/libFoX_fsys.a" )
SET( MKL_FLAGS  "-mkl" )
SET( IFORT_FLAGS "-fpp" )
SET( CMAKE_EXE_LINKER_FLAGS "${FOX_FLAGS} ${MKL_FLAGS}" )

# Define some directories
# SET(LIB ${CMAKE_SOURCE_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

add_library(Lattice2d Lattice2d.f90)
set_target_properties(
	Lattice2d PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(LatticeVector LatticeVector.f90)
set_target_properties(
	LatticeVector PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(ClassVectorDoubleArray1d ClassVectorDoubleArray1d.f90)
set_target_properties(
	ClassVectorDoubleArray1d PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(ClassVectorComplex8 ClassVectorComplex8.f90)
set_target_properties(
	ClassVectorComplex8 PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(ClassVectorDouble ClassVectorDouble.f90)
set_target_properties(
	ClassVectorDouble PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(ClassVectorInteger ClassVectorInteger.f90)
set_target_properties(
	ClassVectorInteger PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(ClassVectorIntegerArray1d ClassVectorIntegerArray1d.f90)
set_target_properties(
	ClassVectorIntegerArray1d PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(ClassInterlayerAttenuation ClassInterlayerAttenuation.f90)
set_target_properties(
	ClassInterlayerAttenuation PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(ModuleHermitianMatrixEigenvalueSolver ModuleHermitianMatrixEigenvalueSolver.f90)
set_target_properties(
	ModuleHermitianMatrixEigenvalueSolver PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS} ${MKL_FLAGS}"
)

add_library(Class_kBase Class_kBase.f90)
target_link_libraries(Class_kBase ClassVectorInteger ClassVectorIntegerArray1d)
set_target_properties(
	Class_kBase PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(ClassSingleWallNanotube ClassSingleWallNanotube.f90)
target_link_libraries(ClassSingleWallNanotube Class_kBase Lattice2d LatticeVector)
set_target_properties(
	ClassSingleWallNanotube PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(ClassNanotubeBundle ClassNanotubeBundle.f90)
target_link_libraries(ClassNanotubeBundle ClassSingleWallNanotube ClassInterlayerAttenuation Lattice2d LatticeVector)
set_target_properties(
	ClassNanotubeBundle PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(ClassHamiltonianConstructor ClassHamiltonianConstructor.f90)
target_link_libraries(ClassHamiltonianConstructor Class_kBase ClassNanotubeBundle ClassInterlayerAttenuation ClassVectorComplex8)
set_target_properties(
	ClassHamiltonianConstructor PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_library(ClassSpectralFunction ClassSpectralFunction.f90)
target_link_libraries(ClassSpectralFunction Class_kBase ClassVectorInteger ClassVectorDouble ClassVectorDoubleArray1d)
set_target_properties(
	ClassSpectralFunction PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)

add_executable(NanotubeBundleDeformed main.f90)
target_link_libraries(NanotubeBundleDeformed 
	Lattice2d 
	LatticeVector
	ClassSingleWallNanotube
	ClassHamiltonianConstructor
	ModuleHermitianMatrixEigenvalueSolver
	Class_kBase 
	ClassNanotubeBundle 
	ClassSpectralFunction
	ClassVectorDouble 
	ClassVectorDoubleArray1d)
set_target_properties(
	NanotubeBundleDeformed PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS} ${FOX_FLAGS}"
)

add_executable(TestClassInterlayerAttenuation TestClassInterlayerAttenuation.f90)
target_link_libraries(TestClassInterlayerAttenuation ClassInterlayerAttenuation)
set_target_properties(
	TestClassInterlayerAttenuation PROPERTIES
		COMPILE_FLAGS
			"${IFORT_FLAGS}"
)