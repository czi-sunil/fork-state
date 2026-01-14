#!/usr/bin/env bash
#
# Reproducing results from SCG-LLM:
# Train the State model on SCG-LLM datasets, then run Prediction.
#
# See also: scg-llm-state_run.sh
# 
#
# The code automatically uses all GPUs, if available.
#
# NOTE: Set your WANDB API Key in the following var first!
#	WANDB_API_KEY
#	WANDB_BASE_URL ... default value is "https://czi.wandb.io"
#
#       Edit the Options below before using!
#

# Exit on any error
set -e

# -- Get path to this script

# Get the absolute full path of the script, resolving any symlinks
SCRIPT_FULL_PATH="$(readlink -f "${BASH_SOURCE[0]}")"

# Extract the directory from the absolute path
SCRIPT_DIR="$(dirname "$SCRIPT_FULL_PATH")"

CMD=`basename -- "${BASH_SOURCE[0]}"`

start_date=`date`
echo "Start: ${start_date}"
echo


# -- Options: command-line

# SCG-LLM: 80000
MAXSTEPS=80000

# SCG-LLM: 1000
VALSTEPS=1000

# Resume from Checkpoint?
RESUME=


# -- SCG-LLM Options: edit here

# Gio overrides this to 128
# Leave unset if using the setting in the MODEL.yaml config
HIDDEN_DIM=

# 'parse1m' or 'replogle'
DATASET=parse1m

FOLD=1

# train_seed
SEED=12345


# -- Other Options: edit here

# All runs are under this
RUNDIR=./Runs


TOML_BASE_DIR=/mnt/vcm-perturbation-v1-10t/giovanni/state-reproduce/baselines/state_sets_reproduce/configs/splits

# Path to this experiment's data
# ... this is just the base-dir.
# The TOML points to the data, but EXPERIMENT requires a valid path.
DATADIR=/mnt/vcm-perturbation-v1-10t/czb_data/

# Experiment name (base name of config file)
EXPERIMENT=scg_llm_state

# Run name
RUNNAME="${EXPERIMENT}-${DATASET}-f${FOLD}-s${SEED}"

# Run sub-dir
RUN_SUBDIR=scg_llm

#
# Note:
#    Output is in: ${OUTPUTDIR}/${RUNNAME}/eval_last.ckpt/"
#


#
# Set dataset-specific variables
#
if [[ "${DATASET}" == "parse1m" ]]; then
  CELL_TYPE_KEY="cell_type"
  PERT_COL="cytokine"
  CONTROL_PERT="PBS"
  BATCH_COL="donor"
  TOML_CONFIG_PATH=${TOML_BASE_DIR}/parse1m_state.toml
else
  CELL_TYPE_KEY="cell_line"
  PERT_COL="gene"
  CONTROL_PERT="non-targeting"
  BATCH_COL="gem_group"
  TOML_CONFIG_PATH=${TOML_BASE_DIR}/replogle_rpe1_state.toml
fi


# -- Options from command line

function ShowOpts {
    echo "SCG-LLM Run with following options:"
    echo
    echo "DATASET = ${DATASET}"
    echo "EXPERIMENT = ${EXPERIMENT}"
    echo
    echo "MAXSTEPS = ${MAXSTEPS}"
    echo "VALSTEPS = ${VALSTEPS}"
    echo
    echo "RUNNAME = ${RUNNAME}"
    echo "RUN_SUBDIR = ${RUN_SUBDIR}"
    echo
}

function Usage {
    echo "Usage: ${CMD} [-h] [-m MAX_STEPS] [-v VALIDATION_FREQ_STEPS] [-r]"
    echo
    echo "Defaults: MAX_STEPS=${MAXSTEPS}, VALIDATION_FREQ_STEPS=${VALSTEPS}"
    echo
    echo "  MAX_STEPS = max nbr batches for training"
    echo "  VALIDATION_FREQ_STEPS = Validation frequency (nbr batches), and checkpoint freq."
    echo
    echo "  -r = Resume from last checkpoint, if possible"
    echo "       otherwise Runs/${RUN_SUBDIR}/${RUNNAME} gets deleted"
    echo
    echo "Edit the script for other options shown below."
    echo

    ShowOpts
}


OPTSTRING=":hm:v:r"

while getopts "$OPTSTRING" opt; do
  case ${opt} in
      h)
      Usage
      exit 0
      ;;
    m)
      MAXSTEPS=${OPTARG}
      ;;
    v)
      VALSTEPS=${OPTARG}
      ;;
    r)
      RESUME=1
      ;;
    \?)
      # Handles invalid options (e.g., -x)
      Usage
      exit 1
      ;;
    :)
      # Handles missing arguments for options that require them (e.g., -m without a value)
      Usage
      exit 1
      ;;
  esac
done

# Shift processed options so that positional arguments remain in $@
shift $((OPTIND - 1))


ShowOpts


# -- wandb

if [[ -z "$WANDB_API_KEY" && "${SCRIPT_DIR}" == /mnt/*/sunil/* ]]; then
    
    MY_INIT_SCR="/mnt/vcm-perturbation-v1/sunil/cluster.sh"

    if [ -e "${MY_INIT_SCR}" ]; then
	source "${MY_INIT_SCR}"
    fi
fi


if [ -z "$WANDB_API_KEY" ]; then
    echo "WANDB_API_KEY is not set!"
fi

if [ -n "$WANDB_API_KEY" ] && [ -z "$WANDB_BASE_URL" ]; then
    export WANDB_BASE_URL=https://czi.wandb.io
fi


# -- Paths

OUTPUTDIR="${RUNDIR}/${RUN_SUBDIR}"

TRNG_LOGFILE="${RUNDIR}/log_${RUNNAME}.txt"


# -- Invoke venv

cd "${SCRIPT_DIR}"

source ${SCRIPT_DIR}/.venv/bin/activate


# -- Check paths

if [ ! -d "${RUNDIR}" ]; then
    echo "Creating dir: ${RUNDIR}"
    mkdir -p $RUNDIR
fi

if [ ! -d "${DATADIR}" ]; then
    echo "Data dir ${DATADIR} not found!"
    exit 1
fi

# If not Resuming from Checkpoint and RUNNAME exists then delete it

if [[ -z "${RESUME}" && -d "${OUTPUTDIR}/${RUNNAME}" ]]; then
    echo "Clearing old RUNNAME dir:  ${OUTPUTDIR}/${RUNNAME}"
    rm -rf "${OUTPUTDIR}/${RUNNAME}"
fi

if [ -n "${RESUME}" ]; then
    echo
    echo "Resuming from checkpoint, if possible"
    echo
fi


# Optional overrides

OVERRIDES=

if [ ! -z "${HIDDEN_DIM}" ]; then
    OVERRIDES="${OVERRIDES}  model.kwargs.hidden_dim=${HIDDEN_DIM}"
fi


# -- Capture all remaining output to TRNG_LOGFILE

echo "Remaining logs captured in: ${TRNG_LOGFILE}"

exec &> "${TRNG_LOGFILE}"

ShowOpts


# --  Train

echo "Starting training ..."
echo

# Checkpoint after every validation

uv run state tx train \
  +experiment=${EXPERIMENT} \
  datadir="${DATADIR}" \
  data.kwargs.toml_config_path=${TOML_CONFIG_PATH} \
  data.kwargs.pert_col=${PERT_COL} \
  data.kwargs.cell_type_key=${CELL_TYPE_KEY} \
  data.kwargs.control_pert=${CONTROL_PERT} \
  data.kwargs.batch_col=${BATCH_COL} \
  training.train_seed=${SEED} \
  training.max_steps=${MAXSTEPS} \
  training.val_freq=${VALSTEPS} \
  training.ckpt_every_n_steps=${VALSTEPS} \
  training.devices=auto \
  ${OVERRIDES} \
  wandb.tags="[${RUNNAME}]" \
  wandb.project=scg-llm-state \
  wandb.entity="" \
  output_dir="${OUTPUTDIR}" \
  name="${RUNNAME}"

if [ $? -ne 0 ]; then
    echo
    echo "Training failed!"
    echo
    exit 1
fi

trng_copmlete_date=`date`

echo
echo "   Training completed"
echo
echo "   Started at:   ${start_date}"
echo "   Completed at: ${trng_copmlete_date}"
echo "-------------------------"
echo


# -- Predict and score

echo "Starting prediction ..."
echo

uv run state tx predict --output-dir "${OUTPUTDIR}/${RUNNAME}" --checkpoint "last.ckpt"

if [ $? -ne 0 ]; then
    echo
    echo "Predict failed!"
    echo
    exit 1
fi

echo
echo "   Predictions and Metrics completed"
echo "   Output is in: ${OUTPUTDIR}/${RUNNAME}/eval_last.ckpt/"
echo "-------------------------"
echo

echo "-- ${CMD} --"
echo "All started at:        ${start_date}"
echo "Prediction started at: ${trng_copmlete_date}"
echo "All Completed at:    " `date`
echo
