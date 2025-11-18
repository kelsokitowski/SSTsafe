#!/bin/bash
#SBATCH --job-name=SST_sim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --output=sst_%j.out
#SBATCH --error=sst_%j.err

#===============================================================================
# SLURM script for compiling and running SST simulation
# with fixed array dimensions based on kLength
#
# This script:
#   1. Sets kLength from environment or uses default
#   2. Compiles the program with that kLength
#   3. Runs the simulation
#===============================================================================

# Set kLength - you can override this in your environment
# Example: sbatch --export=kLength=200 slurm_compile_and_run.sh
if [ -z "$kLength" ]; then
    # Default value if not set
    export kLength=144
    echo "kLength not set in environment, using default: $kLength"
else
    echo "Using kLength from environment: $kLength"
fi

# Load required modules (adjust for your cluster)
# module load openmpi
# module load gcc

# Print environment info
echo "========================================="
echo "SLURM Job Information"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Tasks: $SLURM_NTASKS"
echo "kLength: $kLength"
echo "========================================="
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Compile the program
echo "Step 1: Compiling program with kLength=$kLength..."
echo "Current directory: $(pwd)"

# Try using the simpler compile.sh first, fallback to compile_with_klength.sh
if [ -f compile.sh ]; then
    bash compile.sh
else
    bash compile_with_klength.sh
fi

if [ $? -ne 0 ]; then
    echo "Compilation failed! Exiting."
    exit 1
fi

echo ""
echo "Step 2: Running simulation..."
echo "========================================="

# Run the program
mpirun -np $SLURM_NTASKS ./mainSSTProgram.x

echo ""
echo "========================================="
echo "Simulation complete!"
echo "========================================="
