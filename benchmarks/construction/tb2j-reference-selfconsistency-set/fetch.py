"""Refetch the pinned source files for tb2j-reference-selfconsistency-set and extract the members used.

    python3 benchmarks/construction/tb2j-reference-selfconsistency-set/fetch.py [--dest DIR]

Source: A. Tellez-Mora, X. He, E. Bousquet, L. Wirtz, A. H. Romero, "Systematic determination of a material's
magnetic ground state from first principles", Materials Cloud Archive record rbrcp-0tc34 (2024),
DOI 10.24435/materialscloud:5m-2t, CC BY 4.0 (paper: npj Comput. Mater. 10, 20 (2024), doi 10.1038/s41524-024-01202-z).

Four archive files are used (12.1 MB in total, fetched whole; far below the 3 GB cap): NiO, FeCl2, FePS3 and FeP.
The sibling task tb2j-groundstate-cuo-multiferroic uses CuO.tar.gz from the same record; this task does not.
Checksums: md5 from the archive API, sha256 computed on download (2026-10-19); member sha256 below.
Members are extracted to DEST/raw/<material>/<member>.
"""
import argparse
import hashlib
import tarfile
import urllib.request
from pathlib import Path

URL = 'https://archive.materialscloud.org/api/records/rbrcp-0tc34/files/{}/content'
ARCHIVES = {  # file -> (size, md5, sha256)
    'NiO.tar.gz': (327069, '02b5fac1c748ee9e612d6d803efaed71', '994cbf17a204dff1cdfdb9dd2844db2e3c94635ea2c932f49018ee360e8f8a09'),
    'FeCl2.tar.gz': (1312750, '653cd5334436fb3b69e43352704083bf', 'ab4f4bb5dd88558325a408e384943a8739841a82972024e23f4889a2c118716f'),
    'FePS3.tar.gz': (5750718, '54ebb117ece8326050225bf6fc141d78', 'dca986c02a035c4d213638fb9a72315e09c4ac54fc9dffe1c277b7f9ab2a5cdc'),
    'FeP.tar.gz': (4700382, '9ca0793b00444416ef5a7ea596c741e0', '27ca2d78df076f0a6a6d0e876bda44f2ec97528b36cbfcdfd8d394c3287405b4'),
}
MEMBERS = {  # (archive, member) -> sha256.  TB2J.yaml files duplicate exchange.out and are not used.
    ('NiO.tar.gz', 'README'): '9740c7754bb665db04465bb0e56e57af3b23d61ed13abfb5fd8c1c415ff072b2',
    ('NiO.tar.gz', 'FM/exchange.out'): 'e7b3998367543a476a511607a236805a6175a17ec3acac1054f58cbcbbd04a03',
    ('NiO.tar.gz', 'AFM/exchange.out'): '39bda9368a59f5532ae6e13b0cbecf3c4108c0c54fa0de335e33e941a3e21021',
    ('FeCl2.tar.gz', 'README'): '9740c7754bb665db04465bb0e56e57af3b23d61ed13abfb5fd8c1c415ff072b2',
    ('FeCl2.tar.gz', 'FM/exchange.out'): 'a42b30acf81d2ccc672be10b380ab701345be5515af8f79662a36ed567eb054c',
    ('FeCl2.tar.gz', 'AFM/exchange.out'): 'c4eefc608958a309c7304ced94d0d8863ce2d0dd0288fbc3f85ff42d2ef1ee2b',
    ('FeCl2.tar.gz', 'FM/dispersion_FM'): '010b4d555c8a803a245c4fd01c842a5c012365a72ac266ad051c9568f1378053',
    ('FeCl2.tar.gz', 'AFM/dispersion_AFM'): '440255ae8a38191d16b50d2e3332289110f367f16246fa38d972cd36f3ad92d6',
    ('FePS3.tar.gz', 'README'): '9740c7754bb665db04465bb0e56e57af3b23d61ed13abfb5fd8c1c415ff072b2',
    ('FePS3.tar.gz', 'DFT/FM/exchange.out'): 'c229e2e49079b9626b2ca59dc2915b0deb637f9fbb1dbce619378806adbf08ad',
    ('FePS3.tar.gz', 'DFT/AFM/exchange.out'): '15f3cea9c384cfade7bc4a3ea1d3ca68c780a38d145427842a98eb0d2deb4c16',
    ('FePS3.tar.gz', 'DFT+U/FM/exchange.out'): '10901ec776160d4b9c79ee4279407fb2fa59a180c3955f64870c5684c232bd86',
    ('FePS3.tar.gz', 'DFT+U/AFM/exchange.out'): '384cc4150fc74cefc48f0322add2c5715b8d67202236be349eb78a97495fe918',
    ('FePS3.tar.gz', 'DFT/FM/dispersion_FM'): '734a6a131e8bca905883c758f314a04000e8b6e47378caf52c4c3780aaab18b0',
    ('FeP.tar.gz', 'README'): '453787392000154bad82e536dfa25803dc0297f0cc724eed71487cd5c01afe60',
    ('FeP.tar.gz', 'FM/exchange.out'): 'a8f4adec5f612e3bff6ceec5d6e47535686b00a5b5987e12c3d919d7639b5707',
    ('FeP.tar.gz', 'AFM/exchange.out'): 'f1144a28638c1228892031d3c3f879ab224c1b9a3ce3eee8543fb213eb40b4eb',
    ('FeP.tar.gz', 'AFM/dispersion_AFM'): '5e278550b1bb7944102c29356758d99ad0fe8f10705ed356e1830b8447405742',
}
DEFAULT = Path.home()/'.cache/ai-theorist/tb2j-reference-selfconsistency-set'


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
    for name, (size, md5, sha) in ARCHIVES.items():
        archive = a.dest/name
        if not archive.exists() or digest(archive, 'md5') != md5:
            urllib.request.urlretrieve(URL.format(name), archive)
        assert archive.stat().st_size == size, name
        assert digest(archive, 'md5') == md5 and digest(archive, 'sha256') == sha, 'checksum mismatch: ' + name
        material = name.split('.')[0]
        with tarfile.open(archive) as t:
            for (arc, member) in MEMBERS:
                if arc == name:
                    t.extract(member, a.dest/'raw'/material)
    for (arc, member), sha in MEMBERS.items():
        path = a.dest/'raw'/arc.split('.')[0]/member
        assert digest(path, 'sha256') == sha, 'checksum mismatch: ' + str(path)
    print('fetched and verified; extracted', len(MEMBERS), 'members under', a.dest/'raw')


if __name__ == '__main__':
    main()
