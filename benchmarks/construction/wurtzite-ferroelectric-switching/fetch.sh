#!/bin/bash
# Refetch the pinned archive file (Zenodo record 7916565, CC BY 4.0) and unpack it.
set -euo pipefail
DEST="${1:-$HOME/.cache/ai-theorist/wurtzite-ferroelectric-switching}"
mkdir -p "$DEST"
cd "$DEST"
curl -sSL -o AlBN-system-Science.zip https://zenodo.org/api/records/7916565/files/AlBN-system-Science.zip/content
echo "814d04e76cfb237221c961154c5e8716  AlBN-system-Science.zip" | md5sum -c -
sha256sum AlBN-system-Science.zip
mkdir -p raw && cd raw && unzip -q -o ../AlBN-system-Science.zip
