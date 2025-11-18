#!/bin/bash
# Simple compilation wrapper - no execution permissions needed
# Can be called with: bash compile.sh

if [ -z "$kLength" ]; then
    echo "ERROR: kLength not set!"
    echo "Usage: export kLength=144 && bash compile.sh"
    exit 1
fi

echo "Compiling with kLength=$kLength"

mpif90 -O3 -DkLength=$kLength -c \
    array_dimensions.F90 \
    integrationModules.F90 \
    kernel_interp.F90 \
    EDQNMstratifiedModules.F90 \
    job_parameters.F90 \
    data_loader_all.F90 \
    getForcing.F90 \
    timeRoutines.F90 \
    mainSST.F90

if [ $? -eq 0 ]; then
    mpif90 -O3 -DkLength=$kLength -o mainSSTProgram.x \
        array_dimensions.o \
        integrationModules.o \
        kernel_interp.o \
        EDQNMstratifiedModules.o \
        job_parameters.o \
        data_loader_all.o \
        getForcing.o \
        timeRoutines.o \
        mainSST.o

    if [ $? -eq 0 ]; then
        echo "SUCCESS: Executable created: mainSSTProgram.x"
    else
        echo "ERROR: Linking failed"
        exit 1
    fi
else
    echo "ERROR: Compilation failed"
    exit 1
fi
