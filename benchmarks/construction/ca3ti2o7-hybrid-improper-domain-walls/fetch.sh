#!/bin/bash
# Refetch the Zenodo archive (CC BY 4.0) and unpack it into the cache used by the construction scripts.
set -euo pipefail
D=$HOME/.cache/ai-theorist/ca3ti2o7-hybrid-improper-domain-walls
mkdir -p $D && cd $D
curl -sL -o Data.zip https://zenodo.org/api/records/22801877/files/Data.zip/content
echo "4007d33af8f3b2b26055b4fbc711abeb  Data.zip" | md5sum -c -
unzip -q -o Data.zip -d raw
