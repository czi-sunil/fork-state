#!/usr/bin/env bash
#
# Train the State model on VCC data curated by CZI
# 
#
# The code automatically uses GPU if available.
#
# NOTE: Set your WANDB API Key in the following var first!
#	WANDB_API_KEY
#	WANDB_BASE_URL ... default value is "https://czi.wandb.io"
#

# Exit on any error
# set -e

# --

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

CMD=`basename -- "${BASH_SOURCE[0]}"`

RUNNAME='batch_test'

LOGFILE="${SCRIPT_DIR}/log_${RUNNAME}.txt"


echo "test_batch.sh:  SCRIPT_DIR = ${SCRIPT_DIR}"
echo "test_batch.sh:  LOGFILE = ${LOGFILE}"


# Capture all remaining output to LOGFILE

exec &> "${LOGFILE}"


start_date=`date`
echo "Start: ${start_date}"
echo


# -- cd

cd ${SCRIPT_DIR}

echo "Changed dir to:" `pwd`
echo


# -- Get WANDB params

echo "Setting up env ..."

source /mnt/vcm-perturbation-v1/sunil/cluster.sh

# -- Get the env

source ./.venv/bin/activate

echo "Done env"
echo "SCRIPT_DIR = '${SCRIPT_DIR}'"
echo

# -- Options

# Recommended: >= 40000
MAXSTEPS=4000

# Recommended: 20000
CKPTSTEPS=200

# Run name
RUNNAME="vcc-czi-batch-run"

function ShowOpts {
    echo "Running with following options:"
    echo "MAXSTEPS = $MAXSTEPS"
    echo "CKPTSTEPS = $CKPTSTEPS"
    echo
}


# -- wandb

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

ShowOpts

echo "Args: $*"
echo

nvidia-smi

echo "Done."
