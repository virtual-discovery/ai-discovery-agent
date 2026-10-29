#!/bin/bash
# Refetch the Zenodo 22031607 archives used by extract.py (about 765 MB) and unpack them.
set -euo pipefail
D=${1:-$HOME/.cache/ai-theorist/her-alloy-workfunction-descriptor}
mkdir -p "$D/x" && cd "$D"
for f in 01-dataset 02-models 03-inference 04-activity; do
  curl -s -L -o $f.tar.gz https://zenodo.org/api/records/22031607/files/$f.tar.gz/content
done
md5sum -c <<'SUMS'
eff9518f36a306a6db04347c5b427032  01-dataset.tar.gz
4f184f3bd3798bf60da36aa27ee0a624  02-models.tar.gz
77eacf62b4e5d107a018649fc55b9ea9  03-inference.tar.gz
1541a46cecb4c58b3aea0852a233079a  04-activity.tar.gz
SUMS
cd x
tar xzf ../01-dataset.tar.gz --exclude='01-dataset/data/raw/*'
tar xzf ../03-inference.tar.gz --exclude='*.pkl'
tar xzf ../04-activity.tar.gz --exclude='04-activity/runs/diagnostics/lolo-crossfit/fold-predictions.csv.gz'
