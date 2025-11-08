#!/usr/bin/env bash
#
# Train the State model on VCC-curated dataset
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

# Recommended: >= 40000
MAXSTEPS=4000

# Recommended: 20000
VALSTEPS=200


# -- Options: edit here

# Path to this experiment's data
DATADIR=../../Data/Arc/vcc_curated

# Experiment name (base name of config file)
EXPERIMENT=vcc_state_sm

# Run sub-dir
RUN_SUBDIR="vcc_gx"

# Run name
RUNNAME="${EXPERIMENT}-batch-run"


# -- Options from command line

function ShowOpts {
    echo "Running with following options:"
    echo "MAXSTEPS = ${MAXSTEPS}"
    echo "VALSTEPS = ${VALSTEPS}"
    echo
    echo "RUNNAME = ${RUNNAME}"
    echo "RUN_SUBDIR = ${RUN_SUBDIR}"
    echo
    echo "EXPERIMENT = ${EXPERIMENT}"
    echo "DATADIR = ${DATADIR}"
    echo
}

function Usage {
    echo "Usage: ${CMD} [-h] [-m MAX_STEPS] [-v VALIDATION_FREQ_STEPS]"
    echo
    echo "Defaults: MAX_STEPS=${MAXSTEPS}, VALIDATION_FREQ_STEPS=${VALSTEPS}"
    echo
    echo "  MAX_STEPS = max nbr batches for training"
    echo "  VALIDATION_FREQ_STEPS = Validation frequency (nbr batches), and checkpoint freq."
    echo
    echo "Edit the script for other options shown below."
    echo

    ShowOpts
}


OPTSTRING=":hm:v:"

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

RUNDIR=./Runs

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

# If RUNNAME exists then delete it

if [ -d "${OUTPUTDIR}/${RUNNAME}" ]; then
    echo "Clearing old RUNNAME dir:  ${OUTPUTDIR}/${RUNNAME}"
    rm -rf "${OUTPUTDIR}/${RUNNAME}"
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
  training.max_steps=${MAXSTEPS} \
  training.val_freq=${VALSTEPS} \
  training.ckpt_every_n_steps=${VALSTEPS} \
  training.devices=auto \
  wandb.tags="[${RUNNAME}]" \
  wandb.project=pert-state \
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
