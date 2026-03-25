#!/bin/bash
#SBATCH --job-name=cortical_features_retry
#SBATCH --account=rrg-akhanf_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=0-11
#SBATCH --output=/scratch/apooladi/halawork/work/sb_features/logs/%x_%A_%a.out
#SBATCH --error=/scratch/apooladi/halawork/work/sb_features/logs/%x_%A_%a.err

# --- 12 subjects that failed due to multiple T1w files in BIDS ---
SUBJECTS=(
  LHS0001 LHS0005 LHS0008 LHS0012 LHS0013 LHS0020
  LHS5216 LHS5217 LHS5220 LHS5221 LHS5224 LHS5226
)

SUBJECT=${SUBJECTS[$SLURM_ARRAY_TASK_ID]}

echo "=== Running cortical features for sub-${SUBJECT} (array index ${SLURM_ARRAY_TASK_ID}) ==="
echo "Host: $(hostname)"
echo "Start: $(date)"

# --- Load modules ---
module load freesurfer/8.0.0-1
module load connectome-workbench
source $EBROOTFREESURFER/FreeSurferEnv.sh
export FS_LICENSE=/project/ctb-akhanf/akhanf/opt/freesurfer/.license

# --- Activate pixi environment ---
export PATH="/home/apooladi/.pixi/bin:$PATH"
eval "$(pixi shell-hook --manifest-path /home/apooladi/features/pixi.toml)"

# --- Paths ---
BIDS_DIR=/scratch/apooladi/halawork/work/bids
RECON_DIR=/scratch/habudaqq/deepprep/output
FINAL_DIR=/scratch/apooladi/halawork/work/sb_features
WORK_DIR=${FINAL_DIR}/.work/sub-${SUBJECT}
PIPELINE_DIR=/home/apooladi/features

# --- Clean previous results for this subject ---
rm -rf "${FINAL_DIR}/sub-${SUBJECT}"
rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"

# --- Run pipeline ---
cd "${PIPELINE_DIR}"
python run.py \
  "${BIDS_DIR}" \
  "${WORK_DIR}" \
  participant \
  --participant-label "${SUBJECT}" \
  --deepprep_dir "${RECON_DIR}" \
  --skip-bids-validation \
  --force-output \
  --cores "${SLURM_CPUS_PER_TASK}" \
  --freesurfer_layout deepprep

RC=$?
echo "Exit code: ${RC}"

# --- Move results to final output dir ---
if [ ${RC} -eq 0 ] && [ -d "${WORK_DIR}/sub-${SUBJECT}" ]; then
    cp -a "${WORK_DIR}/sub-${SUBJECT}" "${FINAL_DIR}/"
    echo "Copied results to ${FINAL_DIR}/sub-${SUBJECT}"
fi

echo "End: $(date)"
