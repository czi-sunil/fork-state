#!/usr/bin/env bash
#
# Train the State model on VCC data provided by the State team
#	Downloaded from 
#	https://storage.googleapis.com/vcc_data_prod/datasets/state/competition_support_set.zip
# 
# and then produce metrics.
#
# The code automatically uses GPU if available.
#
# NOTE: Set your WANDB API Key in the following var first!
#	WANDB_API_KEY
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
RUNNAME="vcc-comp-run"


# -- wandb

if [ -z "$WANDB_API_KEY" ]; then
    echo "WANDB_API_KEY is not set!"
fi

if [ -z "$WANDB_BASE_URL" ]; then
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


# -- Chdir and invoke venv

echo "Changing dir to ${SCRIPT_DIR}"
echo
cd $SCRIPT_DIR

source .venv/bin/activate

# -- Paths

RUNDIR=./Runs

DATADIR=../../Data/Arc/competition_support_set

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
  data.kwargs.toml_config_path="${DATADIR}/starter.toml" \
  data.kwargs.num_workers=8 \
  data.kwargs.batch_col="batch_var" \
  data.kwargs.pert_col="target_gene" \
  data.kwargs.cell_type_key="cell_type" \
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

echo
echo "   Training completed"
echo "-------------------------"
echo


# -- Predict, produce output in VCC format

echo "Starting prediction ..."

uv run state tx infer \
  --output "${OUTPUTDIR}/prediction.h5ad" \
  --model-dir "${OUTPUTDIR}/${RUNNAME}" \
  --checkpoint "${OUTPUTDIR}/${RUNNAME}/checkpoints/last.ckpt" \
  --adata "${DATADIR}/competition_val_template.h5ad" \
  --pert-col "target_gene"


# ... now convert predictions into VCC submission format

echo
echo "Preparing VCC submission format ..."

uv run cell-eval prep -i ${OUTPUTDIR}/prediction.h5ad -g ${DATADIR}/gene_names.csv

echo
echo "   Predictions completed"
echo "   Output is in: ${OUTPUTDIR}/prediction.prep.vcc"
echo "-------------------------"
echo


echo "-- ${CMD} --"
echo "Started at: ${start_date}"
echo "All Completed at:" `date`
echo
