"""Refetch the pinned Materials Cloud record p8297-g9n39 (DOI 10.24435/materialscloud:s5-a8, CC BY 4.0).

    python3 -B fetch.py [CACHE]     (default ~/.cache/ai-theorist/mnwo4-halfmoon-frustration/raw)

Downloads the four zip archives the task uses through the record's per-file API links, checks each md5
against the record's own checksums (listed below), and unpacks only the members extract.py reads.
Total download 113 MB; unpacked members about 200 MB.
"""
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/p8297-g9n39/files/%s/content'
FILES = {  # name: (bytes, md5) from https://archive.materialscloud.org/api/records/p8297-g9n39/files
    'Fig_2.zip': (15325551, 'aa3841d1f4b7be0718f4cf502ff12a9f'),             # paper Fig. 1b,c: 1.5 K S(q,w), model inputs
    'Fig_4.zip': (27103293, '9520e07c20ddde03d70f42770ec7890c'),             # paper Fig. 3: S(q) maps at 14/20/30/60 K
    'Fig_7.zip': (51334042, 'c2080da562d727f4651f863e0c7a040e'),             # paper Fig. 5b-f: 20 K maps of reduced models
    'Fig_5_and_Fig_10.zip': (19989965, 'eb1da6949f380e604eb8610030a2d924'),  # paper SI Figs. 7/9: 11 K full-zone S(q)
}
NEEDED = ('jfile', 'kfile', 'posfile', 'momfile', 'inpsd.dat', 'qfile', 'sq.loopTemp.out', 'sqw.loopTemp.out',
          'swdos.loopTemp.out', 'qpointsdir.out', 'totenergy.loopTemp.out', 'README.txt')


def main():
    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/mnwo4-halfmoon-frustration/raw'
    dest.mkdir(parents=True, exist_ok=True)
    for name, (size, md5) in FILES.items():
        path = dest/name
        if not path.exists() or hashlib.md5(path.read_bytes()).hexdigest() != md5:
            with urllib.request.urlopen(RECORD % name, timeout=600) as r:
                path.write_bytes(r.read())
        data = path.read_bytes()
        assert len(data) == size and hashlib.md5(data).hexdigest() == md5, name
        with zipfile.ZipFile(path) as z:
            for member in z.namelist():
                if member.rsplit('/', 1)[-1] in NEEDED and '/qfiles_1d/' not in member:
                    z.extract(member, dest/'x')
    print('verified and unpacked', len(FILES), 'archives into', dest)


if __name__ == '__main__':
    main()
