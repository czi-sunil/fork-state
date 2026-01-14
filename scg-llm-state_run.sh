#!/bin/bash

#
# This file copied from:
#    https://github.com/czi-ai/scg-vae/blob/giovp/perturb/benchmark/state_run.sh
#

# Re-exec with bash if invoked with sh (to support [[ ]] and ==)
if [ -z "${BASH_VERSION:-}" ]; then
    exec bash "$0" "$@"
fi

# Script to run STATE training
# Usage: state_run.sh <dataset> <fold> <seed>
# Example: state_run.sh replogle 1 12345

DATASET=$1
FOLD=$2
SEED=$3

# Set dataset-specific variables
if [[ "${DATASET}" == "parse1m" ]]; then
  CELL_TYPE_KEY="cell_type"
  PERT_COL="cytokine"
  CONTROL_PERT="PBS"
  BATCH_COL="donor"
  TOML_CONFIG_PATH="/mnt/czi-sci-ai/project-scg-llm-pvc/giovanni/state-reproduce/baselines/state_sets_reproduce/configs/splits/parse1m_state.toml"
else
  CELL_TYPE_KEY="cell_line"
  PERT_COL="gene"
  CONTROL_PERT="non-targeting"
  BATCH_COL="gem_group"
  TOML_CONFIG_PATH="/mnt/czi-sci-ai/project-scg-llm-pvc/giovanni/state-reproduce/baselines/state_sets_reproduce/configs/splits/replogle_rpe1_state.toml"
fi

BASE_DIR="/mnt/czi-sci-ai/project-scg-llm-data-2/state_baselines"

RUN_NAME="state-${DATASET}-f${FOLD}-s${SEED}"
echo "Running ${RUN_NAME}"


source /mnt/czi-sci-ai/project-scg-llm-pvc/giovanni/envstate/bin/activate && \
# Run state tx train with all parameters
state tx train \
  data.kwargs.toml_config_path="${TOML_CONFIG_PATH}" \
  data.kwargs.embed_key=X_hvg \
  data.kwargs.num_workers=12 \
  data.kwargs.pert_col=${PERT_COL} \
  data.kwargs.cell_type_key=${CELL_TYPE_KEY} \
  data.kwargs.control_pert=${CONTROL_PERT} \
  data.kwargs.batch_col=${BATCH_COL} \
  training.max_steps=80000 \
  training.val_freq=1000 \
  training.ckpt_every_n_steps=10000 \
  training.batch_size=64 \
  training.lr=1e-3 \
  training.train_seed=${SEED} \
  model.kwargs.hidden_dim=128 \
  model.kwargs.batch_encoder=True \
  model=state \
  wandb.entity="giovp" \
  wandb.project="state-reproduce" \
  wandb.tags="[${DATASET}-fold${FOLD}]" \
  output_dir="$BASE_DIR/state_${DATASET}_fold${FOLD}_seed${SEED}" \
  name="state-${DATASET}-f${FOLD}-s${SEED}"

echo "Completed ${RUN_NAME}"
