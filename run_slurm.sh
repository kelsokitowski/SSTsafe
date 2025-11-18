#!/bin/bash
#SBATCH --job-name=SST_sim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --output=sst_%j.out
#SBATCH --error=sst_%j.err

#===============================================================================
# SLURM script for SST simulation with fixed arrays
#
# Set kLength before submitting:
#   export kLength=144 && sbatch run_slurm.sh
# Or:
#   sbatch --export=kLength=144 run_slurm.sh
#===============================================================================

# Set default kLength if not provided
if [ -z "$kLength" ]; then
    export kLength=144
    echo "Using default kLength=144"
fi

# Print job info
echo "========================================="
echo "SLURM Job Information"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Tasks: $SLURM_NTASKS"
echo "kLength: $kLength"
echo "Working directory: $(pwd)"
echo "========================================="
echo ""

# Load modules if needed (uncomment and adjust for your cluster)
# module load openmpi
# module load gcc

# Compile the program
echo "STEP 1: Compiling program..."
echo "========================================="
bash compile_mpi.sh

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Compilation failed!"
    echo "Check the error messages above."
    exit 1
fi

# Run the simulation
echo ""
echo "STEP 2: Running simulation..."
echo "========================================="
mpirun -np $SLURM_NTASKS ./mainSSTProgram.x

EXIT_CODE=$?

echo ""
echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "Simulation completed successfully!"
else
    echo "Simulation failed with exit code $EXIT_CODE"
fi
echo "========================================="

exit $EXIT_CODE
