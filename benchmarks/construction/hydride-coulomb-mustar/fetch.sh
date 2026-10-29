#!/bin/bash
# Refetch the members of Zenodo record 23235290 used by hydride-coulomb-mustar.
# The 892 MB tarball is streamed once (never stored); only the listed members are written,
# then verified against members.md5. Peak disk ~50 MB.
#   usage: bash fetch.sh [DEST]   (default ~/.cache/ai-theorist/hydride-coulomb-mustar)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DEST=${1:-$HOME/.cache/ai-theorist/hydride-coulomb-mustar}
URL=https://zenodo.org/api/records/23235290/files/mantel_data_files_v1.tar.gz/content
MD5=31b5f98dafc279ce627f7a1473d08d3d
mkdir -p "$DEST"
cd "$DEST"
cut -c35- "$HERE/members.md5" > members.list
curl -sL "$URL" | tee >(md5sum | cut -c1-32 > stream.md5) \
  | tar -xzf - --files-from=members.list 2> >(grep -v 'unknown extended header' >&2 || true)
sleep 2
[ "$(cat stream.md5)" = "$MD5" ] || { echo "tarball md5 mismatch: $(cat stream.md5)"; exit 1; }
md5sum --quiet -c "$HERE/members.md5"
echo "verified $(wc -l < members.list) members in $DEST"
