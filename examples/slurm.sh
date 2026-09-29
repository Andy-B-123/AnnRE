#!/usr/bin/env bash
# Supply account, partition, resources, array and log path when submitting.
# Example: sbatch --array=0-71 --cpus-per-task=1 --mem=8G --time=01:00:00 \
#   --account=YOUR_ACCOUNT --partition=YOUR_PARTITION examples/slurm.sh /absolute/path/to/run
set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?Submit this script as a SLURM array}"
RUN=${1:?Provide the prepared run directory}
annre extract --run "$RUN" --library "$SLURM_ARRAY_TASK_ID"
# Run `annre report --run /absolute/path/to/run` only after all tasks succeed.
