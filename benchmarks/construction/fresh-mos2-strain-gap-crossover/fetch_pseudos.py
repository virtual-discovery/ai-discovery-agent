"""Refetch the four pinned PSLibrary 1.0.0 pseudopotentials for fresh-mos2-strain-gap-crossover and vendor them.

    python3 benchmarks/construction/fresh-mos2-strain-gap-crossover/fetch_pseudos.py [--install]

Files (A. Dal Corso, "Pseudopotentials periodic table: From H to Pu", Comput. Mater. Sci. 95, 337 (2014);
PSLibrary 1.0.0, generated with ld1.x v6.3, 2018), downloaded from the Quantum ESPRESSO pseudopotential library
https://pseudopotentials.quantum-espresso.org/upf_files/ and checked against the SHA256 and MD5 pinned below:

- Mo.pbe-spn-rrkjus_psl.1.0.0.UPF, S.pbe-n-rrkjus_psl.1.0.0.UPF          scalar-relativistic ultrasoft (PBE)
- Mo.rel-pbe-spn-rrkjus_psl.1.0.0.UPF, S.rel-pbe-n-rrkjus_psl.1.0.0.UPF  fully relativistic ultrasoft (PBE, with
  the spin-orbit information), generated from the same ld1.x inputs with rel=2

License: GNU General Public License v2 (PSLibrary, https://github.com/dalcorso/pslibrary, file `License`, fetched
and pinned below; vendored as environment/lab/PSLIBRARY_License.txt).
Downloads go to ~/.cache/ai-theorist/fresh-mos2-strain-gap-crossover/; --install copies them into
environment/lab/ (lab image only; the agent image never sees them).
"""
import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

BASE = 'https://pseudopotentials.quantum-espresso.org/upf_files/'
FILES = {
    'Mo.pbe-spn-rrkjus_psl.1.0.0.UPF': ('03f09b12e49601d1a456a78ce8b2d514',
                                        '0d7c57996e624698242e6428a34d5f13a1e52b4cddac2ce215be4dce429073dd', 1888440),
    'S.pbe-n-rrkjus_psl.1.0.0.UPF': ('9d0fe39aa2a799ea84730884c9b010f4',
                                     '90fb585e830674f2de7b4cbace28d5dd0c05ca1a919d1f92a59b54d6b95cbfc4', 1526090),
    'Mo.rel-pbe-spn-rrkjus_psl.1.0.0.UPF': ('efb2e050edcba4ab52b616d1147efb03',
                                            'a39da9404fae8865b491eef86edf5f4dfec45fc994031b493ea7247814653e39', 3505235),
    'S.rel-pbe-n-rrkjus_psl.1.0.0.UPF': ('425a098c3228d5d3ee792f5a4250baba',
                                         '800bb92fdf3054b5e71aace4ae19a3e4b1053b1a9d74361f832b6dcc6200bf75', 3202013),
}
LICENSE = ('https://raw.githubusercontent.com/dalcorso/pslibrary/master/License', 'PSLIBRARY_License.txt',
           '204d8eff92f95aac4df6c8122bc1505f468f3a901e5a4cc08940e0ede1938994')
CACHE = Path.home()/'.cache/ai-theorist/fresh-mos2-strain-gap-crossover'
LAB = Path(__file__).resolve().parents[2]/'tasks/fresh-mos2-strain-gap-crossover/environment/lab'


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fetch(url, target, sha):
    if not target.exists() or digest(target, 'sha256') != sha:
        print('downloading', url)
        with urllib.request.urlopen(url, timeout=600) as r, open(target, 'wb') as out:
            shutil.copyfileobj(r, out)
    assert digest(target, 'sha256') == sha, f'SHA256 mismatch for {target.name}'


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    for name, (md5, sha, size) in FILES.items():
        fetch(BASE + name, CACHE/name, sha)
        assert digest(CACHE/name, 'md5') == md5 and (CACHE/name).stat().st_size == size, name
        print(name, size, 'md5', md5, 'sha256', sha)
    url, name, sha = LICENSE
    fetch(url, CACHE/name, sha)
    print(name, 'sha256', sha)
    if '--install' in sys.argv:
        for name in list(FILES) + [LICENSE[1]]:
            shutil.copy2(CACHE/name, LAB/name)
        print('installed into', LAB)


if __name__ == '__main__':
    main()
