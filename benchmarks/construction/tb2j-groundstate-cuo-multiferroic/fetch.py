"""Refetch the pinned source file for tb2j-groundstate-cuo-multiferroic and extract the members used.

    python3 benchmarks/construction/tb2j-groundstate-cuo-multiferroic/fetch.py [--dest DIR]

Source: A. Tellez-Mora, X. He, E. Bousquet, L. Wirtz, A. H. Romero, "Systematic determination of a material's
magnetic ground state from first principles", Materials Cloud Archive record rbrcp-0tc34 (2024),
DOI 10.24435/materialscloud:5m-2t, CC BY 4.0 (paper: npj Comput. Mater. 10, 20 (2024), doi 10.1038/s41524-024-01202-z).

Only CuO.tar.gz (6.5 MB) is used, so it is fetched whole (far below the 3 GB cap). The sibling task
tb2j-reference-selfconsistency-set uses FeP, FeCl2, FePS3 and NiO from the same record; this task uses none of them.
Checksums: md5 from the archive API, sha256 computed on download (2026-10-19); member sha256 below.
"""
import argparse
import hashlib
import tarfile
import urllib.request
from pathlib import Path

URL = 'https://archive.materialscloud.org/api/records/rbrcp-0tc34/files/CuO.tar.gz/content'
SIZE, MD5 = 6482521, 'b57a770ecd719fe8e9713b727fb1290b'
SHA256 = '6d2507643f20c18873e8a4de656c9bddebdaf1bcbfa83dcdd63b4dec80650a71'
MEMBERS = {  # member -> sha256 (TB2J.yaml files are not used: they duplicate exchange.out, cross-checked to 5e-5 meV)
    'README': '9740c7754bb665db04465bb0e56e57af3b23d61ed13abfb5fd8c1c415ff072b2',
    'ligand/FM/exchange.out': '62e8a6224326ab841b4dea5d864d63086261410fd01012067b8ecedffbf69389',
    'ligand/FM/dispersion_FM': 'a00cc80a1e03e4173245f84ceda2e21d1c56cded0d50343726b884550bf30a11',
    'ligand/FM/CuO_FM.xsf': '0e3e4f9a24b7716d626d7743c1399698ba311a0f5d188fadb7d5c6ad792068af',
    'ligand/AFM/exchange.out': '717f7c012e7a87582c5c97de2bba917340440577ec129789052300beda93835b',
    'ligand/AFM/dispersion_AFM': '431437fc34de7421de23f2d61ca37ef04b71080a900e102d2ab2ab81f36189d9',
    'ligand/AFM/CuO_AFM.xsf': 'baf4dcb6247183df78a07e7667e310344e12823f7a059f61676739e0780d4095',
    'noligand/FM/exchange.out': '3182c2163097b45a8c620cd00bdc5da29c440347e077f9d23d5e93719ea00fe0',
    'noligand/FM/dispersion_nol_FM': 'ee3c27496726ecf4610a40ea324187122839405e2acea469c6fac58017ab7d4e',
    'noligand/AFM/exchange.out': '1988f5136cb14cda8ede2554066cbb0bab97a3fb670a7349504d378d02978afd',
    'noligand/AFM/dispersion_nol_AFM': '537e276dc491744fa733646f6bd73f8719a0e28e9182e64093b043d6b689520b',
}
DEFAULT = Path.home()/'.cache/ai-theorist/tb2j-groundstate-cuo-multiferroic'


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dest', type=Path, default=DEFAULT)
    a = p.parse_args()
    a.dest.mkdir(parents=True, exist_ok=True)
    archive = a.dest/'CuO.tar.gz'
    if not archive.exists() or digest(archive, 'md5') != MD5:
        urllib.request.urlretrieve(URL, archive)
    assert archive.stat().st_size == SIZE
    assert digest(archive, 'md5') == MD5 and digest(archive, 'sha256') == SHA256, 'checksum mismatch: CuO.tar.gz'
    raw = a.dest/'raw'
    with tarfile.open(archive) as t:
        for member in MEMBERS:
            t.extract(member, raw)
    for member, sha in MEMBERS.items():
        assert digest(raw/member, 'sha256') == sha, f'checksum mismatch: {member}'
    print('fetched and verified; extracted', len(MEMBERS), 'members under', raw)


if __name__ == '__main__':
    main()
