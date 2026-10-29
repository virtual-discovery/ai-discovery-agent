#!/bin/bash
# Refetch the pinned source archive (Zenodo record 20591788, CC BY 4.0) and build the construction venv.
set -euo pipefail
CACHE=${CACHE:-$HOME/.cache/ai-theorist/pdh-miscibility-isotope-ace}
mkdir -p "$CACHE" && cd "$CACHE"
curl -sL -o pdh.zip "https://zenodo.org/api/records/20591788/files/Pd-H_ACE_potential_files_and_training_data.zip/content"
echo "0308b4e0ed87805975d37eaeaacf976a  pdh.zip" | md5sum -c -
unzip -oq pdh.zip 'potential_files/ace/*'
sha256sum potential_files/ace/PdH_PBE_v03a.yaml potential_files/ace/df_train.pckl.gzip potential_files/ace/df_test.pckl.gzip
# construction-only tools (host): python-ace for the C-tilde conversion and the validation reference, pandas/ase for the pickles
python3 -m venv venv && venv/bin/pip -q install python-ace==0.4.0rc1 pandas ase
