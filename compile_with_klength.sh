#!/bin/bash

#===============================================================================
# Compilation script for SST program with fixed array dimensions
# Usage:
#   1. Set kLength environment variable: export kLength=144
#   2. Run: ./compile_with_klength.sh
#
# Or in a SLURM script:
#   export kLength=144
#   ./compile_with_klength.sh
#===============================================================================

# Check if kLength is set in environment
if [ -z "$kLength" ]; then
    echo "ERROR: kLength environment variable is not set!"
    echo "Please set it before compiling, e.g.:"
    echo "  export kLength=144"
    echo "  $0"
    exit 1
fi

echo "========================================="
echo "Compiling SST program with kLength=$kLength"
echo "========================================="

# Compiler settings
FC=mpif90
FFLAGS="-O3 -DkLength=$kLength"
FFLAGS_DEBUG="-g -O0 -DkLength=$kLength -fcheck=all -fbacktrace"

# Choose compilation mode (comment/uncomment as needed)
MODE="optimized"
# MODE="debug"

if [ "$MODE" == "debug" ]; then
    CURRENT_FLAGS="$FFLAGS_DEBUG"
    OUTPUT="mainSST_debug.x"
else
    CURRENT_FLAGS="$FFLAGS"
    OUTPUT="mainSSTProgram.x"
fi

echo "Compilation mode: $MODE"
echo "Compiler flags: $CURRENT_FLAGS"
echo "Output executable: $OUTPUT"
echo ""

# Source files in compilation order
SOURCES="
array_dimensions.F90
integrationModules.F90
kernel_interp.F90
EDQNMstratifiedModules.F90
job_parameters.F90
data_loader_all.F90
getForcing.F90
timeRoutines.F90
mainSST.F90
"

# Clean old module files
echo "Cleaning old .mod and .o files..."
rm -f *.mod *.o

# Compile
echo "Compiling source files..."
$FC $CURRENT_FLAGS -c $SOURCES

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Compilation failed!"
    exit 1
fi

# Link
echo "Linking executable..."
OBJECT_FILES=$(echo $SOURCES | sed 's/\.F90/\.o/g')
$FC $CURRENT_FLAGS -o $OUTPUT $OBJECT_FILES

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Linking failed!"
    exit 1
fi

echo ""
echo "========================================="
echo "SUCCESS! Compilation completed."
echo "Executable: $OUTPUT"
echo "Compiled with kLength = $kLength"
echo "========================================="
echo ""
echo "To run:"
echo "  mpirun -np <num_procs> ./$OUTPUT"
