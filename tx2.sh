#!/usr/bin/env bash
#
# Train the State model on VCC data curated by CZI.
# NOTE: First use `run_preprocess_trng.sh` to prepare the data from raw counts.
# 
#
# The code automatically uses GPU if available.
#
# NOTE: Set your WANDB API Key in the following var first!
#	WANDB_API_KEY
#	WANDB_BASE_URL ... default value is "https://czi.wandb.io"
#

# -- Get path to this script

# Get the absolute full path of the script, resolving any symlinks
SCRIPT_FULL_PATH="$(readlink -f "${BASH_SOURCE[0]}")"

# Extract the directory from the absolute path
SCRIPT_DIR="$(dirname "$SCRIPT_FULL_PATH")"

CMD=`basename -- "${BASH_SOURCE[0]}"`

start_date=`date`
echo "Start: ${start_date}"
echo


# -- Options

# Recommended: >= 40000
MAXSTEPS=100

# Recommended: 20000
CKPTSTEPS=25

# Run name
RUNNAME="tx2-test"


# -- My init scr, and wandb

if [[ "${SCRIPT_DIR}" == /mnt/* ]]; then
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


# -- Options from command line

OPTSTRING=":hm:c:"

while getopts "$OPTSTRING" opt; do
  case ${opt} in
    h)
      echo "Usage: ${CMD} [-h] [-m MAX_STEPS] [-c CHECKPOINT_STEPS]"
      exit 0
      ;;
    m)
      MAXSTEPS=${OPTARG}
      ;;
    c)
      CKPTSTEPS=${OPTARG}
      ;;
    \?)
      # Handles invalid options (e.g., -x)
      echo "Error: Invalid option -${OPTARG}." >&2
      echo "Usage: ${CMD} [-h] [-m MAX_STEPS] [-c CHECKPOINT_STEPS]" >&2
      exit 1
      ;;
    :)
      # Handles missing arguments for options that require them (e.g., -m without a value)
      echo "Error: Option -${OPTARG} requires an argument." >&2
      echo "Usage: ${CMD} [-h] [-m MAX_STEPS] [-c CHECKPOINT_STEPS]" >&2
      exit 1
      ;;
  esac
done

# Shift processed options so that positional arguments remain in $@
shift $((OPTIND - 1))


function ShowOpts {
    echo "Running with following options:"
    echo "MAXSTEPS = $MAXSTEPS"
    echo "CKPTSTEPS = $CKPTSTEPS"
    echo "RUNNAME = ${RUNNAME}"
    echo
}

ShowOpts

# -- Invoke venv

source ${SCRIPT_DIR}/.venv/bin/activate

# -- Paths

RUNDIR=./Runs

DATADIR=../../Data/Arc/vcc_curated

OUTPUTDIR=${RUNDIR}/vcc_g1

TRNG_LOGFILE="${RUNDIR}/log_${RUNNAME}.txt"


if [ ! -d "${RUNDIR}" ]; then
    mkdir -pv $RUNDIR
fi

if [ ! -d "${DATADIR}" ]; then
    echo "Data dir ${DATADIR} not found!"
    exit 1
fi


# -- Capture all remaining output to TRNG_LOGFILE

# echo
# echo "Output logs sent to: ${TRNG_LOGFILE}"
# echo

# exec &> "${TRNG_LOGFILE}"


ShowOpts


# --  Train

echo "Starting training ..."
echo

#  datadir="${DATADIR}" \

uv run state tx2 train \
  +experiment=vcc_czi \
  datadir="${DATADIR}" \
  training.max_steps=${MAXSTEPS} \
  training.ckpt_every_n_steps=${CKPTSTEPS} \
  training.devices=1 \
  wandb.tags="[${RUNNAME}]" \
  wandb.project=pert-vcc-st \
  wandb.entity="" \
  output_dir="${OUTPUTDIR}" \
  name="${RUNNAME}"

if [ $? -ne 0 ]; then
    echo
    echo "Training failed!"
    echo
    exit 1
fi

echo
echo "   Training completed"
echo "-------------------------"
echo

