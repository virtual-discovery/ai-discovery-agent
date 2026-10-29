"""Refetch the pinned fully relativistic Ge pseudopotential for fresh-ge-soc-strained-valence and vendor it into the lab.

    python3 benchmarks/construction/fresh-ge-soc-strained-valence/fetch_pseudo.py [--install]

Ge_ONCV_PBE_FR-1.0.upf: SG15 ONCV norm-conserving, fully relativistic (D. R. Hamann, PRB 88, 085117 (2013);
M. Schlipf, F. Gygi, Comput. Phys. Commun. 196, 36 (2015); P. Scherpelz, M. Govoni, I. Hamada, G. Galli,
J. Chem. Theory Comput. 12, 3523 (2016)), from the SG15 distribution at quantum-simulation.org, and its scalar-relativistic SG15 sibling
Ge_ONCV_PBE-1.0.upf (Schlipf-Gygi 2015) for the lab's spin_orbit = false runs: QE 6.7 does not j-average the fully
relativistic file when lspinorb is off (construction test: P = -16.5 GPa near a0 and the s-like level 1.3 eV above
the p triplet, both for noncolin = .false. and for noncolin = .true., lspinorb = .false.), so that mode is not offered.
License (stated in the file's PP_INFO): Creative Commons Attribution-ShareAlike 4.0 (CC BY-SA 4.0).

The paired study's Materials Cloud record vz6hg-tvd13 (DOI 10.24435/materialscloud:yg-hd, CC BY 4.0) ships the same
potential in raw_data_inputs_outputs.zip (pseudo/Ge_ONCV_PBE_FR-1.0.upf, 167,311 bytes); it differs from the SG15
file only by one extra line holding the file name after <PP_INFO> (checked during construction). The SG15 file is the
one pinned and vendored here.

Both files carry the same license. Downloads go to ~/.cache/ai-theorist/fresh-ge-soc-strained-valence/. --install copies the files into
environment/lab/ (lab image only; the agent image never sees it).
"""
import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

BASE = 'http://www.quantum-simulation.org/potentials/sg15_oncv/upf/'
# name: (size, md5, sha256, marker that must appear in the header)
FILES = {'Ge_ONCV_PBE_FR-1.0.upf': (167286, '046fb2884f640e00dddb7bc679dd92bf',
                                    'ee520682a4103e45bff9360d87b9fb8f15a2045eb206da5ef412372984fd53af', 'relativistic="full"'),
         'Ge_ONCV_PBE-1.0.upf': (115422, '39f69ca9d3fedd83738352a2cbce5cbd',
                                 '1023970eba91ac075d8a5d85a4817afe98f98affcde8a0dd378a4365105548d7', 'relativistic="scalar"')}
CACHE = Path.home()/'.cache/ai-theorist/fresh-ge-soc-strained-valence'
LAB = Path(__file__).resolve().parents[2]/'tasks/fresh-ge-soc-strained-valence/environment/lab'


def digest(path, algo):
    return hashlib.new(algo, path.read_bytes()).hexdigest()


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    for name, (size, md5, sha, marker) in FILES.items():
        target = CACHE/name
        if not target.exists() or digest(target, 'sha256') != sha:
            print('downloading', BASE + name)
            with urllib.request.urlopen(BASE + name, timeout=120) as r, open(target, 'wb') as out:
                shutil.copyfileobj(r, out)
        assert target.stat().st_size == size, target.stat().st_size
        assert digest(target, 'md5') == md5, 'MD5 mismatch'
        assert digest(target, 'sha256') == sha, 'SHA256 mismatch'
        text = target.read_text()
        assert marker in text and 'Attribution-ShareAlike' in text
        print(name, size, 'md5', md5, 'sha256', sha)
        if '--install' in sys.argv:
            shutil.copy2(target, LAB/name)
            print('installed into', LAB)


if __name__ == '__main__':
    main()
