#!/bin/bash
# Refetch the pinned files of Materials Cloud Archive record xd06z-0pp29 (doi 10.24435/materialscloud:ek-fp, CC BY 4.0)
# and unpack the DFT tarball. Retries because the archive's object store (rgw.cscs.ch) is intermittently 503.
set -uo pipefail
DEST="${1:-$HOME/.cache/ai-theorist/cr2o3-magnetoelectric-sign}"
mkdir -p "$DEST"; cd "$DEST"
API=https://archive.materialscloud.org/api/records/xd06z-0pp29/files
declare -A MD5=([README.txt]=9c6dfac9cdfd791755a6bd86dc4d28af [files_description.md]=cde11635fe27af49edafa52a5bc677e9 \
                [DFT_calculations.tar.gz]=dd3a0fa4aa1e7b990d9fad58f5f02ca8)
for f in README.txt files_description.md DFT_calculations.tar.gz; do
  for attempt in $(seq 1 ${TRIES:-120}); do
    if [ -f "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "${MD5[$f]}" ]; then break; fi
    code=$(curl -sSL -C - -o "$f.part" -w '%{http_code}' "$API/$f/content")
    if [ "$code" = 200 ] || [ "$code" = 206 ]; then mv "$f.part" "$f"; else rm -f "$f.part"; sleep 30; fi
  done
  echo "${MD5[$f]}  $f" | md5sum -c - || exit 1
done
mkdir -p raw && tar -xzf DFT_calculations.tar.gz -C raw
