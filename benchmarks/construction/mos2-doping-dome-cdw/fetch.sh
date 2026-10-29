#!/bin/bash
# Refetch the pinned Materials Cloud archive record 51yyy-ksf72 (DOI 10.24435/materialscloud:z1-aw, CC BY 4.0):
# Girotto Erhardt, Berges, Ponce and Novko, "Understanding the origin of superconducting dome in electron-doped
# MoS2 monolayer", npj 2D Mater. Appl. 9, 44 (2025).  Only the six archives below are used (11 MB in total);
# fig6.tar.gz (the elphmod "model III" supercell model, 97 MB) and figS11.tar.gz are deliberately not fetched.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DEST="${1:-$HOME/.cache/ai-theorist/mos2-doping-dome-cdw}"
API="https://archive.materialscloud.org/api/records/51yyy-ksf72/files"
mkdir -p "$DEST"
cd "$DEST"
fetch() {  # name md5
  if [ ! -f "$1" ] || [ "$(md5sum "$1" | cut -d' ' -f1)" != "$2" ]; then
    curl -sSL "$API/$1/content" -o "$1"
  fi
  echo "$2  $1" | md5sum -c -
}
fetch README.md 3b22839753c99e088435274842d0449d
fetch fig2a-d.tar.gz 4cdb4ac5e83aaa4e3067d6304cb5e97e
fetch fig2e.tar.gz 0514e9cf3bd294ef9f49e5d79912aefb
fetch fig3.tar.gz dbf67323c8b26e984707a9ab512f979c
fetch fig4.tar.gz ad524c1f940d23e3280ac430a32fe286
fetch fig5.tar.gz 25bf433a60498198d5aec2e36ba287bd
fetch input.tar.gz 6f541b537c66c7de1f3f92ef30242a17
mkdir -p x
for f in fig2a-d fig2e fig3 fig4 fig5 input; do tar xzf "$f.tar.gz" -C x; done
(cd x && md5sum -c --quiet "$HERE/members.md5")
echo "unpacked and verified into $DEST/x"
