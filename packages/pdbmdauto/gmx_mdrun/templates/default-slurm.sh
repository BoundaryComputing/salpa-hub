#!/bin/bash -l
#=============================================================================
# GROMACS MD Run - Example SLURM Script Template
#
# This is a template file showing common SLURM directives and GROMACS commands.
# Customize this script for your specific HPC cluster configuration.
#
# Template Variables (automatically replaced by BoCoFlow):
#   {{JOB_NAME}}         - Auto-generated job name
#   {{REMOTE_WORK_DIR}}  - Remote working directory
#   {{OUTPUT_DIR}}       - Output directory (same as REMOTE_WORK_DIR)
#   {{WORKING_DIR}}      - Working directory (same as REMOTE_WORK_DIR)
#   {{NODE_ID}}          - BoCoFlow node ID
#   {{RUN_LABEL}}        - User-specified run label (e.g., md, nvt, npt)
#   {{INPUT_TOP_FILE}}   - Topology file basename
#   {{INPUT_GRO_FILE}}   - Structure file basename
#   {{INPUT_MDP_FILE}}   - Parameters file basename
#   {{INPUT_NDX_FILE}}   - Index file basename (may be empty)
#
# Usage:
#   1. Copy this template to your SLURM job script field in BoCoFlow
#   2. Modify the SLURM directives for your cluster
#   3. Adjust module loads to match your environment
#   4. Keep the {{VARIABLE}} placeholders - they will be replaced at runtime
#=============================================================================

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{OUTPUT_DIR}}/slurm-%j.out
#SBATCH --error={{OUTPUT_DIR}}/slurm-%j.err

# Resource allocation - customize for your cluster
#SBATCH --partition=gpu          # Change to your partition name
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8      # Number of MPI tasks per node
#SBATCH --cpus-per-task=2        # OpenMP threads per task
#SBATCH --gres=gpu:1             # Number of GPUs per node
#SBATCH --time=24:00:00          # Maximum wall time
#SBATCH --mem=32G                # Memory per node

# Email notifications (optional)
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=your.email@example.com

#-----------------------------------------------------------------------------
# Environment Setup - Customize for your cluster
#-----------------------------------------------------------------------------

# Load GROMACS module - adjust for your cluster's module system
# Examples for different clusters:
module load GROMACS/2023.3              # Generic
# module load gromacs/2023.3-cuda       # With CUDA support
# module load gromacs/2023/intel-mpi    # With Intel MPI

# Set OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#-----------------------------------------------------------------------------
# Change to working directory
#-----------------------------------------------------------------------------
cd {{WORKING_DIR}}

echo "============================================"
echo "GROMACS MD Simulation - {{RUN_LABEL}}"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"
echo "============================================"

#-----------------------------------------------------------------------------
# Run grompp to generate TPR file
#-----------------------------------------------------------------------------
echo "Running grompp..."

# Build grompp command
GROMPP_CMD="gmx grompp \
    -f {{INPUT_MDP_FILE}} \
    -c {{INPUT_GRO_FILE}} \
    -r {{INPUT_GRO_FILE}} \
    -p {{INPUT_TOP_FILE}} \
    -o {{RUN_LABEL}}.tpr \
    -maxwarn 10"

# Add index file if provided
if [ -n "{{INPUT_NDX_FILE}}" ]; then
    GROMPP_CMD="$GROMPP_CMD -n {{INPUT_NDX_FILE}}"
fi

# Execute grompp
$GROMPP_CMD

if [ $? -ne 0 ]; then
    echo "ERROR: grompp failed"
    exit 1
fi

#-----------------------------------------------------------------------------
# Run mdrun
#-----------------------------------------------------------------------------
echo "Running mdrun..."

# For MPI-enabled GROMACS (recommended for multi-node jobs):
srun gmx_mpi mdrun -deffnm {{RUN_LABEL}} -v

# Alternative for single-node jobs without MPI:
# gmx mdrun -deffnm {{RUN_LABEL}} -v -ntmpi 8 -ntomp $OMP_NUM_THREADS

# Alternative with GPU acceleration:
# srun gmx_mpi mdrun -deffnm {{RUN_LABEL}} -v -nb gpu -bonded gpu -pme gpu

if [ $? -ne 0 ]; then
    echo "ERROR: mdrun failed"
    exit 1
fi

#-----------------------------------------------------------------------------
# Completion
#-----------------------------------------------------------------------------
echo "============================================"
echo "Simulation {{RUN_LABEL}} completed"
echo "End time: $(date)"
echo "============================================"

# List output files
ls -la {{RUN_LABEL}}.*
