"""Download and verify the MnTe files of Zenodo record 18863304 (CC BY 4.0).

    python3 -B benchmarks/construction/mnte-altermagnet-ahc-doping/fetch.py [DEST]

Source: K. P. Kluczyk et al., "Coexistence of anomalous Hall effect and weak net magnetization in collinear
antiferromagnet MnTe", Phys. Rev. B 110, 155201 (2024), arXiv:2310.09134; data DOI 10.5281/zenodo.18863304.

Downloads dft_data.zip (4,006,544 bytes) and fig6_script.zip (93,091 bytes) through the record's per-file links,
checks their md5, extracts every member and checks it against members.md5. fig.zip (7.2 MB, figures only) is
not needed: its figS10.dat repeats the withheld sigma_ahc.txt curve on a wider energy window.
DEST defaults to ~/.cache/ai-theorist/mnte-altermagnet-ahc-doping; members land in DEST/<member path>.
"""
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

RECORD = 'https://zenodo.org/api/records/18863304/files'
FILES = {'dft_data.zip': 'cef0196e73cec5deefad7e707137b940',
         'fig6_script.zip': '15e2f717c3d2c28f11c19025122769ce'}
HERE = Path(__file__).resolve().parent


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    dest = Path(sys.argv[1] if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/mnte-altermagnet-ahc-doping')
    dest.mkdir(parents=True, exist_ok=True)
    for name, checksum in FILES.items():
        target = dest/name
        if not target.exists() or md5(target) != checksum:
            urllib.request.urlretrieve(f'{RECORD}/{name}/content', target)
        if md5(target) != checksum:
            raise SystemExit(f'{name}: md5 mismatch')
    members = dict(reversed(line.split(None, 1)) for line in (HERE/'members.md5').read_text().splitlines() if line.strip())
    for name in FILES:
        with zipfile.ZipFile(dest/name) as archive:
            for member in archive.namelist():
                if member.endswith('/'):
                    continue
                archive.extract(member, dest)
                if md5(dest/member) != members[member]:
                    raise SystemExit(f'{member}: md5 mismatch')
    print(f'verified {len(FILES)} files and {len(members)} members in {dest}')


if __name__ == '__main__':
    main()
