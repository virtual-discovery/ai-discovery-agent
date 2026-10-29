"""Range-request access to members of the remote Materials Cloud zips (no full download).

    python3 remotezip.py list CaF2.zip > CaF2.list.json
    python3 remotezip.py get CaF2.zip <member> <outfile>

Used by fetch.py. Works with zip and zip64 central directories; members are stored or deflated.
"""
import json
import struct
import sys
import time
import zlib
import requests

RECORD = 'https://archive.materialscloud.org/api/records/yzhmw-vje16/files'


def url(name):
    return f'{RECORD}/{name}/content'


def fetch(name, start, end, tries=6):
    for k in range(tries):
        try:
            r = requests.get(url(name), headers={'Range': f'bytes={start}-{end}'}, timeout=300)
            if r.status_code == 206:
                return r.content
            raise RuntimeError(f'HTTP {r.status_code}')
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(5 * (k + 1))


def total_size(name):
    r = requests.get(url(name), headers={'Range': 'bytes=0-0'}, timeout=60)
    return int(r.headers['Content-Range'].split('/')[-1])


def central_directory(name):
    size = total_size(name)
    tail = fetch(name, max(0, size - 70000), size - 1)
    i = tail.rfind(b'PK\x05\x06')
    _, _, _, _, n, cd_size, cd_off, _ = struct.unpack('<IHHHHIIH', tail[i:i + 22])
    if cd_off == 0xFFFFFFFF or n == 0xFFFF:
        j = tail.rfind(b'PK\x06\x07')
        z64_off = struct.unpack('<Q', tail[j + 8:j + 16])[0]
        z = fetch(name, z64_off, z64_off + 55)
        n, cd_size, cd_off = struct.unpack('<QQQ', z[32:56])
    cd = fetch(name, cd_off, cd_off + cd_size - 1)
    entries, p = [], 0
    while p < len(cd) and cd[p:p + 4] == b'PK\x01\x02':
        (method, crc, csize, usize, nlen, xlen, clen, lho) = (
            struct.unpack('<H', cd[p + 10:p + 12])[0], *struct.unpack('<III', cd[p + 16:p + 28]),
            *struct.unpack('<HHH', cd[p + 28:p + 34]), struct.unpack('<I', cd[p + 42:p + 46])[0])
        fname = cd[p + 46:p + 46 + nlen].decode('utf-8', 'replace')
        extra = cd[p + 46 + nlen:p + 46 + nlen + xlen]
        q = 0
        while q + 4 <= len(extra):
            hid, hlen = struct.unpack('<HH', extra[q:q + 4])
            if hid == 1:
                vals = list(struct.unpack('<' + 'Q' * (hlen // 8), extra[q + 4:q + 4 + 8 * (hlen // 8)]))
                if usize == 0xFFFFFFFF:
                    usize = vals.pop(0)
                if csize == 0xFFFFFFFF:
                    csize = vals.pop(0)
                if lho == 0xFFFFFFFF:
                    lho = vals.pop(0)
            q += 4 + hlen
        entries.append({'name': fname, 'method': method, 'crc': crc, 'csize': csize, 'usize': usize, 'offset': lho})
        p += 46 + nlen + xlen + clen
    return entries


def extract(name, entry):
    head = fetch(name, entry['offset'], entry['offset'] + 29)
    nlen, xlen = struct.unpack('<HH', head[26:30])
    start = entry['offset'] + 30 + nlen + xlen
    data = fetch(name, start, start + entry['csize'] - 1) if entry['csize'] else b''
    if entry['method'] == 8:
        data = zlib.decompress(data, -15)
    elif entry['method'] != 0:
        raise RuntimeError(f"unsupported method {entry['method']}")
    if zlib.crc32(data) & 0xFFFFFFFF != entry['crc']:
        raise RuntimeError(f"CRC mismatch for {entry['name']}")
    return data


if __name__ == '__main__':
    cmd, name = sys.argv[1], sys.argv[2]
    if cmd == 'list':
        print(json.dumps(central_directory(name)))
    else:
        entry = [e for e in central_directory(name) if e['name'] == sys.argv[3]][0]
        open(sys.argv[4], 'wb').write(extract(name, entry))
