"""Refetch the members of the two Zenodo archives used by hgte-gated-subband-densities.

    python3 -B fetch.py [DEST]          (default ~/.cache/ai-theorist/hgte-gated-subband-densities)

Only single members are read, through HTTP range requests on the remote zip (central
directory plus each member's local header and data). The 0.5 GB data_repository.zip is
never downloaded whole; the 807 MB / 528 MB dispersion CSVs and the Origin projects are
skipped. Each extracted member is checked against members.md5 when that file exists;
with --write-md5 the checksums are (re)written instead.

Archives (both CC BY 4.0):
  https://zenodo.org/records/18709991  data_repository.zip  498,134,382 B  md5 a948d21c70f3f911d3bbc8554acccbce
  https://zenodo.org/records/18429111  data_and_code.zip     57,037,973 B  md5 e777270ff41037b2072b819a8d4b9961
"""
import hashlib
import io
import sys
import zipfile
from pathlib import Path

import requests

HERE = Path(__file__).parent
ARCHIVES = {
    'calc': 'https://zenodo.org/api/records/18709991/files/data_repository.zip/content',
    'meas': 'https://zenodo.org/api/records/18429111/files/data_and_code.zip/content',
}
SKIP = ('dispersion', '.opju')


class RangeFile(io.RawIOBase):
    def __init__(self, url):
        head = requests.head(url, allow_redirects=True, timeout=60)
        self.url, self.size, self.pos = head.url, int(head.headers['Content-Length']), 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def readinto(self, buf):
        if self.pos >= self.size or len(buf) == 0:
            return 0
        end = min(self.pos + len(buf), self.size) - 1
        for _ in range(5):
            try:
                data = requests.get(self.url, headers={'Range': f'bytes={self.pos}-{end}'}, timeout=120).content
                break
            except requests.RequestException:
                continue
        buf[:len(data)] = data
        self.pos += len(data)
        return len(data)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    dest = Path(args[0]) if args else Path.home()/'.cache/ai-theorist/hgte-gated-subband-densities'
    write = '--write-md5' in sys.argv
    known = {}
    if (HERE/'members.md5').exists() and not write:
        for line in (HERE/'members.md5').read_text().splitlines():
            md5, name = line.split('  ', 1)
            known[name] = md5
    out = []
    for tag, url in ARCHIVES.items():
        z = zipfile.ZipFile(io.BufferedReader(RangeFile(url), buffer_size=1 << 20))
        for info in z.infolist():
            if info.is_dir() or any(s in info.filename for s in SKIP):
                continue
            target = dest/tag/info.filename
            if not target.exists() or target.stat().st_size != info.file_size:
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(z.read(info))
            md5 = hashlib.md5(target.read_bytes()).hexdigest()
            name = f'{tag}/{info.filename}'
            if known and known.get(name) != md5:
                raise SystemExit(f'checksum mismatch: {name}')
            out.append(f'{md5}  {name}')
    if write:
        (HERE/'members.md5').write_text('\n'.join(out) + '\n')
    print(f'{len(out)} members in {dest}')


if __name__ == '__main__':
    main()
