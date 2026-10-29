#!/bin/bash
# usage: queue.sh jobsN.txt   (lines: TAG<TAB>JSON); 4 workers pinned to CPUS (default 6 7 22 23)
C=$(cd "$(dirname "$0")" && pwd); cd $C
CPUS=${CPUS:-"6 7 22 23"}
worker() { cpu=$1; while IFS=$'\t' read -r tag json; do [ -f runs/$tag.json ] && continue; ./launch.sh $cpu $tag "$json"; done; }
split -n r/4 -d "$1" qpart_
i=0; for c in $CPUS; do worker $c < qpart_0$i & i=$((i+1)); done; wait; rm -f qpart_0*
echo QUEUE_DONE
