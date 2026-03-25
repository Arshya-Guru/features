#!/bin/bash
#SBATCH --job-name=cortical_features
#SBATCH --account=rrg-akhanf_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=0-136
#SBATCH --output=/scratch/apooladi/halawork/work/sb_features/logs/%x_%A_%a.out
#SBATCH --error=/scratch/apooladi/halawork/work/sb_features/logs/%x_%A_%a.err

# --- Subject list (137 subjects — 7 excluded due to incomplete deepprep) ---
# Excluded (no Recon output): TWH0003 TWH0018 TWH0020 TWH0023 TWH0049 TWH0053 TWH0070
SUBJECTS=(
  HSC0002 HSC0003 HSC0004 HSC0006 HSC0007 HSC0009 HSC0011 HSC0014
  HSC0015 HSC0016 HSC0017 HSC0018 HSC0019 HSC0023 HSC0025 HSC0026
  HSC0028 HSC0030 HSC0034 HSC0035 HSC0037
  LHS0001 LHS0003 LHS0004 LHS0005 LHS0007 LHS0008 LHS0009 LHS0010
  LHS0011 LHS0012 LHS0013 LHS0014 LHS0015 LHS0016 LHS0018 LHS0019
  LHS0020 LHS0021 LHS0022 LHS0025 LHS0027 LHS0029 LHS0031 LHS0034
  LHS0035 LHS0036 LHS0038 LHS0039 LHS0041 LHS0042
  LHS5201 LHS5202 LHS5203 LHS5204 LHS5205 LHS5206 LHS5207 LHS5208
  LHS5209 LHS5211 LHS5213 LHS5214 LHS5215 LHS5216 LHS5217 LHS5219
  LHS5220 LHS5221 LHS5222 LHS5224 LHS5226
  TWH0001 TWH0002 TWH0004 TWH0005 TWH0006 TWH0007 TWH0008 TWH0009
  TWH0011 TWH0012 TWH0013 TWH0014 TWH0017 TWH0019 TWH0021 TWH0022
  TWH0024 TWH0026 TWH0027 TWH0029 TWH0030 TWH0031 TWH0033 TWH0035
  TWH0036 TWH0037 TWH0038 TWH0040 TWH0042 TWH0043 TWH0045 TWH0048
  TWH0050 TWH0051 TWH0052 TWH0054 TWH0056 TWH0057 TWH0058 TWH0059
  TWH0060 TWH0061 TWH0063 TWH0064 TWH0065 TWH0066 TWH0067 TWH0068
  TWH0069 TWH0071 TWH0072 TWH0073 TWH0074 TWH0075 TWH0076 TWH0077
  TWH0078 TWH0079 TWH0080 TWH0081 TWH0082 TWH0083 TWH0084 TWH0085
  TWH0086
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
TEMPLATE_DIR=/scratch/apooladi/halawork/work/templates
FINAL_DIR=/scratch/apooladi/halawork/work/sb_features
WORK_DIR=${FINAL_DIR}/.work/sub-${SUBJECT}
PIPELINE_DIR=/home/apooladi/features

# --- Seed work dir with existing results so snakemake only runs new rules ---
rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"
if [ -d "${FINAL_DIR}/sub-${SUBJECT}" ]; then
    cp -a "${FINAL_DIR}/sub-${SUBJECT}" "${WORK_DIR}/"
fi

# --- Run pipeline (per-subject work dir to avoid race conditions) ---
cd "${PIPELINE_DIR}"
python run.py \
  "${BIDS_DIR}" \
  "${WORK_DIR}" \
  participant \
  --participant-label "${SUBJECT}" \
  --deepprep_dir "${RECON_DIR}" \
  --template_dir "${TEMPLATE_DIR}" \
  --skip-bids-validation \
  --force-output \
  --cores "${SLURM_CPUS_PER_TASK}" \
  --freesurfer_layout deepprep \
  --rerun-incomplete

RC=$?
echo "Exit code: ${RC}"

# --- Copy all results to final output dir ---
if [ ${RC} -eq 0 ] && [ -d "${WORK_DIR}/sub-${SUBJECT}" ]; then
    cp -a "${WORK_DIR}/sub-${SUBJECT}" "${FINAL_DIR}/"
    echo "Copied results to ${FINAL_DIR}/sub-${SUBJECT}"
fi

echo "End: $(date)"
