#!/bin/bash
# Refetch the Materials Cloud archive 238ad-1bg44 (CC BY 4.0). The zip is read in place by extract.py (no unpacking).
set -euo pipefail
D=$HOME/.cache/ai-theorist/ge-pressure-polymorph-pathways
mkdir -p $D && cd $D
curl -sL -o allfiles.zip https://archive.materialscloud.org/api/records/238ad-1bg44/files/allfiles.zip/content
echo "d5faa3ad7e44d8b797ff59d1d7ea66c8  allfiles.zip" | md5sum -c -
