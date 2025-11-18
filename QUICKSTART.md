# SST Fixed Arrays - Quick Start Guide

## Simple Compilation and Running

### 1. Compile the Program

```bash
export kLength=144
bash compile_mpi.sh
```

### 2. Run Locally

```bash
mpirun -np 16 ./mainSSTProgram.x
```

### 3. Run on SLURM

```bash
# Default kLength=144
sbatch run_slurm.sh

# Custom kLength
sbatch --export=kLength=200 run_slurm.sh
```

## That's It!

The program now uses **fixed-size arrays** determined at compile time.

## Important Notes

1. **kLength must match your data files** - The program will check and error if they don't match
2. **Recompile when changing kLength** - Just run `compile_mpi.sh` again with new kLength
3. **Check your data's kLength** - Look at the error message if you get a mismatch

## Scripts

- `compile_mpi.sh` - Single compilation script (uses mpifort)
- `run_slurm.sh` - SLURM job submission script

## Example Workflow

```bash
# 1. Set kLength
export kLength=144

# 2. Compile
bash compile_mpi.sh

# 3. Test locally (optional)
mpirun -np 4 ./mainSSTProgram.x

# 4. Submit to SLURM
sbatch run_slurm.sh
```

## Troubleshooting

**Error: "kLength mismatch"**
- Your data files have different kLength than compiled
- Solution: Check data files and recompile with correct kLength

**Error: "kLength not set"**
- You forgot to set the environment variable
- Solution: `export kLength=144` before compiling

**Error: "mpifort not found"**
- OpenMPI not loaded/installed
- Solution: Load MPI module: `module load openmpi`

For detailed documentation, see `README_FIXED_ARRAYS.md`
