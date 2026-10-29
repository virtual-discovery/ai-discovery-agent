#!/bin/bash
# Refetch the pinned Materials Cloud archive (record qfdfa-cs517, DOI 10.24435/materialscloud:zy-qw, CC BY 4.0)
# and unpack it. strain_GaN.tar.gz is 12,934,938 bytes (md5 9a8747242f4759d235c5552aaea0ba52); README.txt is
# 973 bytes (md5 e72a76b7c622cd39d0d3a91ea75811e6). The unpacked tree is about 150 MB.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DEST="${1:-$HOME/.cache/ai-theorist/gan-hole-piezomobility}"
API="https://archive.materialscloud.org/api/records/qfdfa-cs517/files"
mkdir -p "$DEST"
cd "$DEST"
fetch() {  # name md5
  if [ ! -f "$1" ] || [ "$(md5sum "$1" | cut -d' ' -f1)" != "$2" ]; then
    curl -sSL "$API/$1/content" -o "$1"
  fi
  echo "$2  $1" | md5sum -c -
}
fetch strain_GaN.tar.gz 9a8747242f4759d235c5552aaea0ba52
fetch README.txt e72a76b7c622cd39d0d3a91ea75811e6
mkdir -p x
tar xzf strain_GaN.tar.gz -C x
(cd x && md5sum -c --quiet "$HERE/members.md5")
echo "unpacked and verified into $DEST/x"
