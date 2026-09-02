#!/bin/zsh
#
# Local (zsh) helper: self-adjusting memory throttle for macOS. Blocks until free memory
# exceeds a safe threshold (50% of current free, floored at 500 MB). Sourced by the GNU
# parallel runner scripts (wts_gnu.sh) to avoid oversubscribing memory when running many
# R jobs concurrently.
#

get_free_mem_mb() {
    vm_stat | awk '
        /free/        {free=$3}
        /speculative/ {spec=$3}
        END {print (free+spec)*4096/1024/1024}'
}

# Calculate a safe threshold: 50% of currently free memory or 500 MB, whichever is lower
CURRENT_FREE=$(get_free_mem_mb)
MIN_FREE_MB=$(( CURRENT_FREE / 2 ))
(( MIN_FREE_MB < 500 )) && MIN_FREE_MB=500

# Wait until free memory exceeds threshold
while true; do
    free_mb=$(get_free_mem_mb)
    if (( free_mb > MIN_FREE_MB )); then
        break
    fi
    echo "Waiting for free memory: $free_mb MB < $MIN_FREE_MB MB"  # Optional status
    sleep 1
done
