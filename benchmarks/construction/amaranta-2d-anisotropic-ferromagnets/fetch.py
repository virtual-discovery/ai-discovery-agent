"""Refetch the pinned AMaRaNTA source files (Materials Cloud record wygeh-jrc57, CC BY 4.0).

    python3 fetch.py [CACHE_DIR]

Downloads (about 57 MB in total, never the multi-GB archives in full):
  * data_readable.zip (53,905 bytes, md5 2020bf7c7813bb99f8af19c2030ef0e1), unpacked;
  * the member `db.sqlite3` (AiiDA provenance database) of group1.aiida, group4.aiida,
    group5.aiida, group7.aiida and group10.aiida, read with HTTP range requests from the remote zip (central directory,
    local header, deflated member only). Each member is checked against the CRC32 in the
    archive's central directory and against the sha256 recorded below.
"""
import hashlib
import json
import struct
import subprocess
import sys
import zipfile
import zlib
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/wygeh-jrc57/files/%s/content'
CACHE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/amaranta-2d-anisotropic-ferromagnets'
ZIP = {'name': 'data_readable.zip', 'size': 53905, 'md5': '2020bf7c7813bb99f8af19c2030ef0e1'}
ARCHIVES = {  # archive: (size in bytes, archive md5 from the record, sha256 of the extracted db.sqlite3)
    'group1.aiida': (2228772135, '12d6c201e55cf12df2279b7d4874870d', '1bf3450a671b696c889a79dcef80c1260f155b9e4fa42724659e30f70867d1a5'),
    'group4.aiida': (3455702018, '63429af96ec9978edb12a769fa52c5d5', 'aa4540abf9304336b0343220e2f8bf5e0f3fc81a1ff054b7625487bd1d23b94f'),
    'group5.aiida': (1160855479, 'd977634a0b59b0b5b08ed47cadeb5829', '07898d7f148a70d252600e4c8d7b9e4160b5a69a18926808c5387e5545d8ea94'),
    'group7.aiida': (989415385, '6d7a6de21ba04f0263af09fea094c685', 'a88112265cc4228f3450a303101c84e116d25952149e0b93170866713a4aabf5'),
    'group10.aiida': (862462950, 'f33b7f9621f32aa4e199d15078794d9b', '80e067e68f2127bd08529222122781a3b4f3f660ccac238793fbae39a4d1e9ef'),
}


def signed_url(name):
    return subprocess.run(['curl', '-s', '-o', '/dev/null', '-w', '%{redirect_url}', RECORD % name],
                          capture_output=True, text=True, check=True).stdout


def ranged(url, a, b):
    return subprocess.run(['curl', '-s', '-f', '-r', f'{a}-{b}', url], capture_output=True, check=True).stdout


def zip_member(name, size, member):
    url = signed_url(name)
    tail = ranged(url, size - 65536, size - 1)
    i = tail.rfind(b'PK\x05\x06')
    count, cd_size, cd_off = struct.unpack('<HII', tail[i + 10:i + 20])
    if cd_off == 0xFFFFFFFF:  # zip64 end of central directory record
        j = tail.rfind(b'PK\x06\x06')
        cd_size, cd_off = struct.unpack('<QQ', tail[j + 40:j + 56])
    cd = ranged(url, cd_off, cd_off + cd_size - 1)
    p = 0
    while cd[p:p + 4] == b'PK\x01\x02':
        _, _, _, method, _, _, crc, csize, usize, fl, el, cl, _, _, _, off = struct.unpack('<HHHHHHIIIHHHHHII', cd[p + 4:p + 46])
        fname = cd[p + 46:p + 46 + fl].decode()
        extra, q = cd[p + 46 + fl:p + 46 + fl + el], 0
        while q < len(extra):  # zip64 extended information
            hid, hsz = struct.unpack('<HH', extra[q:q + 4])
            if hid == 1:
                vals, k = extra[q + 4:q + 4 + hsz], 0
                if usize == 0xFFFFFFFF:
                    usize = struct.unpack('<Q', vals[k:k + 8])[0]; k += 8
                if csize == 0xFFFFFFFF:
                    csize = struct.unpack('<Q', vals[k:k + 8])[0]; k += 8
                if off == 0xFFFFFFFF:
                    off = struct.unpack('<Q', vals[k:k + 8])[0]; k += 8
            q += 4 + hsz
        if fname == member:
            head = ranged(url, off, off + 29)
            hfl, hel = struct.unpack('<HH', head[26:30])
            start = off + 30 + hfl + hel
            raw = ranged(url, start, start + csize - 1)
            data = raw if method == 0 else zlib.decompressobj(-15).decompress(raw)
            assert len(data) == usize and zlib.crc32(data) == crc, f'{name}:{member} failed CRC check'
            return data
        p += 46 + fl + el + cl
    raise KeyError(member)


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    zpath = CACHE/ZIP['name']
    if not zpath.exists() or hashlib.md5(zpath.read_bytes()).hexdigest() != ZIP['md5']:
        subprocess.run(['curl', '-sL', '-f', RECORD % ZIP['name'], '-o', str(zpath)], check=True)
    assert hashlib.md5(zpath.read_bytes()).hexdigest() == ZIP['md5']
    zipfile.ZipFile(zpath).extractall(CACHE)
    for name, (size, _, sha) in ARCHIVES.items():
        out = CACHE/(name.split('.')[0] + '_db.sqlite3')
        if not out.exists() or hashlib.sha256(out.read_bytes()).hexdigest() != sha:
            out.write_bytes(zip_member(name, size, 'db.sqlite3'))
        assert hashlib.sha256(out.read_bytes()).hexdigest() == sha, out
        print('ok', out)
    print(json.dumps({'cache': str(CACHE)}))


if __name__ == '__main__':
    main()
