#!/bin/bash
# Refetch the pinned Materials Cloud archive (record pkn2d-h8190, DOI 10.24435/materialscloud:qy-gv, CC BY 4.0)
# and extract only the members this task uses. The zip is 170,896,274 bytes (md5 e838e8716fb42200c34f06119337e909).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DEST="${1:-$HOME/.cache/ai-theorist/li3ps4-superionic-transport}"
URL="https://archive.materialscloud.org/api/records/pkn2d-h8190/files/MaterialsCloudArchive-ChemMat2024.zip/content"
MD5="e838e8716fb42200c34f06119337e909"
mkdir -p "$DEST"
cd "$DEST"
if [ ! -f archive.zip ] || [ "$(md5sum archive.zip | cut -d' ' -f1)" != "$MD5" ]; then
  curl -sSL "$URL" -o archive.zip
fi
echo "$MD5  archive.zip" | md5sum -c -
R=MaterialsCloudArchive/raw-data-4-figures
unzip -q -o archive.zip \
  "$R/Fig8-elec-conductivity/*" "$R/Fig7-transitions-Cp/*" "$R/FigS12-MSD-750K/*" \
  "$R/FigS13-sigma-breakdown/*" "$R/Fig9-SpatialCorrelation/*" "$R/Fig4-CVs/*" \
  "MaterialsCloudArchive/MD-simulations/ML-PBE0/*/input-structures/*" "MaterialsCloudArchive/README" -d x
(cd x && md5sum -c --quiet "$HERE/members.md5")
echo "extracted and verified into $DEST/x"
