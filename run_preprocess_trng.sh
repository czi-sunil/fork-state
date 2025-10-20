#!/usr/bin/env bash
#
# Preprocess CZI curated data for training
# 
#

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

CMD=`basename -- "${BASH_SOURCE[0]}"`

start_date=`date`
echo "Start: ${start_date}"
echo


Usage () {
    local code=${1:-0}
    echo "Usage: ${CMD} [-h] [-n NUM_HVGS] RAW_DATA_FILE OUTPUT_FILE"
    echo
    exit $code
}


# -- Options

# Recommended: 2000
NUM_HVGS=2000

# -- Options from command line

OPTSTRING=":hn:"

while getopts "$OPTSTRING" opt; do
  case ${opt} in
    h)
      Usage
      ;;
    n)
      NUM_HVGS=${OPTARG}
      ;;
    \?)
      # Handles invalid options (e.g., -x)
      echo "Error: Invalid option -${OPTARG}." >&2
      Usage 1
      ;;
    :)
      # Handles missing arguments for options that require them (e.g., -m without a value)
      echo "Error: Option -${OPTARG} requires an argument." >&2
      Usage 1
      ;;
  esac
done

# Shift processed options so that positional arguments remain in $@
shift $((OPTIND - 1))

if [ ${#*} != 2 ]; then
    Usage 1
fi

RAW_DATA="$1"
OUT_DATA="$2"

echo "Running with following options:"
echo "NUM_HVGS = $NUM_HVGS"
echo


# -- Invoke venv

source ${SCRIPT_DIR}/.venv/bin/activate


# --  Preprocess

echo "Processing data ..."
echo

state tx preprocess_train \
  --adata ${RAW_DATA} \
  --output ${OUT_DATA} \
  --num_hvgs ${NUM_HVGS}

if [ $? -ne 0 ]; then
    echo
    echo "Failed to complete preprocessing!"
    echo
    exit 1
fi

echo
echo "   Preprocessing completed"
echo "   Output is in: ${OUT_DATA}"
echo "-------------------------"
echo
