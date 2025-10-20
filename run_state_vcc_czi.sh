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

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

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
RUNNAME="vcc-czi-run"


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


echo "Running with following options:"
echo "MAXSTEPS = $MAXSTEPS"
echo "CKPTSTEPS = $CKPTSTEPS"
echo


# -- Invoke venv

source ${SCRIPT_DIR}/.venv/bin/activate

# -- Paths

RUNDIR=./Runs

DATADIR=../../Data/Arc/vcc_curated

OUTPUTDIR=${RUNDIR}/vcc

TRNG_LOGFILE="${RUNDIR}/log_${RUNNAME}.txt"

if [ ! -d "${RUNDIR}" ]; then
    mkdir -pv $RUNDIR
fi

if [ ! -d "${DATADIR}" ]; then
    echo "Data dir ${DATADIR} not found!"
    exit 1
fi


# --  Train

echo "Starting training ..."

uv run state tx train \
  data.kwargs.toml_config_path="${DATADIR}/statecfg.toml" \
  data.kwargs.num_workers=8 \
  data.kwargs.batch_col="batch" \
  data.kwargs.pert_col="target_gene" \
  data.kwargs.cell_type_key="tissue_ontology_term_id" \
  data.kwargs.control_pert="non-targeting" \
  data.kwargs.perturbation_features_file="${DATADIR}/ESM2_pert_features.pt" \
  training.max_steps=${MAXSTEPS} \
  training.ckpt_every_n_steps=${CKPTSTEPS} \
  model=state_sm \
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


# -- Predict and score

echo "Starting prediction ..."

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
echo "Started at: ${start_date}"
echo "All Completed at:" `date`
echo
