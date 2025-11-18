#!/bin/bash
#===============================================================================
# Simple OpenMPI compilation script for SST program
# Uses mpifort compiler
#
# Usage:
#   export kLength=144
#   bash compile_mpi.sh
#===============================================================================

# Check if kLength is set
if [ -z "$kLength" ]; then
    echo "ERROR: kLength environment variable not set!"
    echo ""
    echo "Usage:"
    echo "  export kLength=144"
    echo "  bash compile_mpi.sh"
    exit 1
fi

echo "========================================="
echo "Compiling SST with OpenMPI"
echo "kLength = $kLength"
echo "========================================="

# Clean previous build
echo "Cleaning old build files..."
rm -f *.o *.mod mainSSTProgram.x mainSST_debug.x

# Compiler and flags
FC=mpifort
FFLAGS="-O3 -DkLength=$kLength"

echo "Compiler: $FC"
echo "Flags: $FFLAGS"
echo ""

# Source files in correct dependency order
echo "Step 1: Compiling modules..."

echo "  Compiling array_dimensions.F90..."
$FC $FFLAGS -c array_dimensions.F90
if [ $? -ne 0 ]; then echo "ERROR compiling array_dimensions.F90"; exit 1; fi

echo "  Compiling integrationModules.F90..."
$FC $FFLAGS -c integrationModules.F90
if [ $? -ne 0 ]; then echo "ERROR compiling integrationModules.F90"; exit 1; fi

echo "  Compiling kernel_interp.F90..."
$FC $FFLAGS -c kernel_interp.F90
if [ $? -ne 0 ]; then echo "ERROR compiling kernel_interp.F90"; exit 1; fi

echo "  Compiling EDQNMstratifiedModules.F90..."
$FC $FFLAGS -c EDQNMstratifiedModules.F90
if [ $? -ne 0 ]; then echo "ERROR compiling EDQNMstratifiedModules.F90"; exit 1; fi

echo "  Compiling job_parameters.F90..."
$FC $FFLAGS -c job_parameters.F90
if [ $? -ne 0 ]; then echo "ERROR compiling job_parameters.F90"; exit 1; fi

echo "  Compiling data_loader_all.F90..."
$FC $FFLAGS -c data_loader_all.F90
if [ $? -ne 0 ]; then echo "ERROR compiling data_loader_all.F90"; exit 1; fi

echo "  Compiling getForcing.F90..."
$FC $FFLAGS -c getForcing.F90
if [ $? -ne 0 ]; then echo "ERROR compiling getForcing.F90"; exit 1; fi

echo "  Compiling timeRoutines.F90..."
$FC $FFLAGS -c timeRoutines.F90
if [ $? -ne 0 ]; then echo "ERROR compiling timeRoutines.F90"; exit 1; fi

echo ""
echo "Step 2: Compiling main program..."
echo "  Compiling mainSST.F90..."
$FC $FFLAGS -c mainSST.F90
if [ $? -ne 0 ]; then echo "ERROR compiling mainSST.F90"; exit 1; fi

echo ""
echo "Step 3: Linking executable..."
$FC $FFLAGS -o mainSSTProgram.x \
    array_dimensions.o \
    integrationModules.o \
    kernel_interp.o \
    EDQNMstratifiedModules.o \
    job_parameters.o \
    data_loader_all.o \
    getForcing.o \
    timeRoutines.o \
    mainSST.o

if [ $? -ne 0 ]; then
    echo "ERROR: Linking failed!"
    exit 1
fi

echo ""
echo "========================================="
echo "SUCCESS!"
echo "Executable created: mainSSTProgram.x"
echo "Compiled with kLength = $kLength"
echo "========================================="
echo ""
echo "To run with MPI:"
echo "  mpirun -np 16 ./mainSSTProgram.x"
echo ""
