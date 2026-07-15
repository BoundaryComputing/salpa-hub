#!/bin/bash -l
#=============================================================================
# EasyParm — ORCA Run — SLURM job script (formal default)
#
# Prefilled into the "SLURM Job Script" field of the ORCA Run node. The node
# substitutes {{VARIABLE}} placeholders at submit time, so customise the
# cluster-specific bits (partition, module tree, memory, scratch) for your
# site and the placeholders will keep working.
#
# Validated on Snellius (SURF) for a 94-atom Sn(IV)-porphyrin DFT
# opt+freq+CHELPG job (B3LYP-D3BJ / def2-SVP + def2-ECP on Sn).
#
# Template variables (substituted by HPCNodeBase at submit time):
#   {{JOB_NAME}}         Auto-generated SLURM job name
#   {{REMOTE_WORK_DIR}}  HPC scratch root from the HPC profile
#   {{OUTPUT_DIR}}       Per-job dir on the remote (== WORKING_DIR)
#   {{WORKING_DIR}}      Per-job dir on the remote
#   {{NODE_ID}}          BoCoFlow node UUID
#   {{RUN_LABEL}}        User-supplied run label (default "orca")
#   {{ORCA_INPUT_FILE}}  {{RUN_LABEL}}.inp
#   {{ORCA_OUTPUT_FILE}} {{RUN_LABEL}}.out
#   {{XYZ_FILE}}         XYZ filename (auto input mode)
#   {{NPROCS}}           %pal nprocs — MUST match --ntasks below
#=============================================================================

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{OUTPUT_DIR}}/{{RUN_LABEL}}-%j.out
#SBATCH --error={{OUTPUT_DIR}}/{{RUN_LABEL}}-%j.err

# ---- Cluster resources (Snellius defaults — adjust for your site) ----------
#SBATCH --partition=genoa            # 128-core EPYC 9654 nodes on Snellius
#SBATCH --nodes=1
#SBATCH --ntasks={{NPROCS}}          # MUST match %pal nprocs in the .inp
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=128G                   # ORCA freq jobs are memory-hungry; the
                                     # default ~1.75 GB/core on genoa OOMs a
                                     # def2-SVP freq on 94 atoms. Don't lower
                                     # without checking %maxcore × NPROCS.
#SBATCH --time=24:00:00

# ----------------------------------------------------------------------------
# Module environment — ORCA via EasyBuild
# ----------------------------------------------------------------------------
# Invoke ORCA via ${EBROOTORCA}/bin/orca (set by EasyBuild) rather than
# $(which orca). If `module load` silently fails, `$(which orca)` expands
# to empty and `$ORCA_BIN file.inp > out` becomes `file.inp > out` — the
# shell tries to run the .inp as a command and the job dies with no useful
# trace.
module purge
module load 2025
module load ORCA/6.1.0-gompi-2025a-avx2

# Pin OpenMP / MKL — ORCA spawns its own MPI; oversubscription tanks throughput.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# ----------------------------------------------------------------------------
# Run ORCA in NODE-LOCAL scratch, not the submit dir
# ----------------------------------------------------------------------------
# The submit dir lives under $HOME, which has a modest per-user quota
# (200 GiB on Snellius). A freq job writes multi-GB intermediates
# (.gbw, .bas0/.bas1, *.densities, integral tmp). If $HOME is near quota
# ORCA fails fast with "Unable to write data in TBasis::WriteElement!".
# $TMPDIR is node-local SSD, multi-TB, purged at job end — copy only the
# small artifacts downstream needs (.out/.hess/.xyz/CHELPG) back.
SUBMIT_DIR="{{WORKING_DIR}}"
SCRATCH="${TMPDIR:-/scratch-local/$USER/$SLURM_JOB_ID}"
mkdir -p "$SCRATCH"

cd "$SUBMIT_DIR"

echo "============================================"
echo "ORCA job: {{RUN_LABEL}}"
echo "Host    : $(hostname)"
echo "Start   : $(date)"
echo "Submit  : $SUBMIT_DIR"
echo "Scratch : $SCRATCH"
echo "ORCA bin: ${EBROOTORCA}/bin/orca"
echo "============================================"

cp -p "$SUBMIT_DIR"/{{ORCA_INPUT_FILE}} "$SUBMIT_DIR"/*.xyz "$SCRATCH"/ 2>/dev/null
cd "$SCRATCH"

# Do NOT wrap in mpirun/srun — ORCA spawns its own MPI from %pal nprocs.
# `tee` mirrors stdout+stderr into the .out file AND the SLURM stdout
# stream so the live log is visible during the run.
${EBROOTORCA}/bin/orca {{ORCA_INPUT_FILE}} 2>&1 | tee {{ORCA_OUTPUT_FILE}}
RC=${PIPESTATUS[0]}

# Copy back only the small artifacts the downstream nodes consume.
# Leave the multi-GB .gbw, .bas*, *.densities, integral tmp in scratch.
for f in {{RUN_LABEL}}.out {{RUN_LABEL}}.hess {{RUN_LABEL}}.xyz \
         {{RUN_LABEL}}.property.txt {{RUN_LABEL}}.chelpg.xyz \
         {{RUN_LABEL}}_atom*.out {{RUN_LABEL}}_atom*.property.txt; do
  [ -f "$f" ] && cp -p "$f" "$SUBMIT_DIR"/ 2>/dev/null
done
cd "$SUBMIT_DIR"

echo "ORCA returncode: $RC"
echo "End   : $(date)"
ls -la {{RUN_LABEL}}.*

exit $RC
