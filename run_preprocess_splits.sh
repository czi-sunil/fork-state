#!/usr/bin/env bash
#
# Preprocess CZI curated data that has been pre-split into training/test
# See `state._cli._tx._preprocess_splits.run_tx_preprocess_splits()` for details.
# 
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


Usage () {
    local code=${1:-0}
    echo "Usage: ${CMD} [-h] [-n NUM_HVGS] [-t TARGET_SUM] TRAINING_SPLIT_FILE TEST_SPLIT_FILE"
    echo
    exit $code
}


# -- Options

# Recommended: 2000
NUM_HVGS=2000

# Follows PerturBench, which also uses 1e4
TARGET_SUM=10000


# -- Options from command line

OPTSTRING=":hn:t:"

while getopts "$OPTSTRING" opt; do
  case ${opt} in
    h)
      Usage
      ;;
    n)
      NUM_HVGS=${OPTARG}
      ;;
    t)
      TARGET_SUM=${OPTARG}
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

TRAINING_SPLIT_FILE="$1"
TEST_SPLIT_FILE="$2"

echo "Running with following options:"
echo "NUM_HVGS = $NUM_HVGS"
echo "TARGET_SUM = $TARGET_SUM"
echo


# -- Invoke venv

source ${SCRIPT_DIR}/.venv/bin/activate


# --  Preprocess

echo "Processing data ..."
echo

uv run state tx preprocess_splits \
  --train_split ${TRAINING_SPLIT_FILE} \
  --test_split ${TEST_SPLIT_FILE} \
  --num_hvgs ${NUM_HVGS} \
  --target_sum ${TARGET_SUM}

status=$?

if [ $status -ne 0 ]; then
    echo
    echo "Failed to complete preprocessing!"
    echo

    if [ $status -eq 137 ]; then
        echo "Exit status = ${status}, likely OOM error"
        echo
    fi
    exit $status
fi

echo
echo "   Preprocessing completed"
echo "   Output is in: " `dirname ${TRAINING_SPLIT_FILE}`
echo "-------------------------"
echo

echo "-- ${CMD} --"
echo "All started at:      ${start_date}"
echo "All Completed at:    " `date`
echo
