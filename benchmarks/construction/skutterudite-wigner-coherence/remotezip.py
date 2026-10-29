"""Minimal HTTP range reader for members of a remote zip (no full download).

Used by fetch.py: reads the central directory of the Materials Cloud archive zip and
extracts single members through HTTP Range requests.
"""
import struct
import time
import zlib

URL = ('https://archive.materialscloud.org/api/records/5jsfj-cdg70/files/'
       'materials_cloud_archive_submission.zip/content')
SIZE = 1213559335


def get_range(start, end, url=URL, tries=5):
    import requests
    for attempt in range(tries):
        try:
            r = requests.get(url, headers={'Range': f'bytes={start}-{end}'}, timeout=600)
            if r.status_code == 206 and len(r.content) == end - start + 1:
                return r.content
            raise IOError(f'HTTP {r.status_code}, {len(r.content)} bytes')
        except (IOError, requests.RequestException):
            if attempt == tries - 1:
                raise
            time.sleep(5 * (attempt + 1))


def central_directory():
    tail = get_range(SIZE - 65536, SIZE - 1)
    i = tail.rfind(b'PK\x05\x06')
    _, _, _, _, n, cd_size, cd_off, _ = struct.unpack('<4s4H2LH', tail[i:i + 22])
    cd = get_range(cd_off, cd_off + cd_size - 1)
    members, p = {}, 0
    for _ in range(n):
        (sig, _, _, flag, method, _, _, crc, csize, usize, nlen, xlen, clen, _, _, _,
         off) = struct.unpack('<4s6H3L5H2L', cd[p:p + 46])
        assert sig == b'PK\x01\x02'
        name = cd[p + 46:p + 46 + nlen].decode('utf-8')
        extra = cd[p + 46 + nlen:p + 46 + nlen + xlen]
        if 0xFFFFFFFF in (csize, usize, off):          # zip64 extra field
            q = 0
            while q < len(extra):
                hid, hlen = struct.unpack('<2H', extra[q:q + 4])
                if hid == 1:
                    vals, r = [], q + 4
                    for field in (usize, csize, off):
                        if field == 0xFFFFFFFF:
                            vals.append(struct.unpack('<Q', extra[r:r + 8])[0]); r += 8
                        else:
                            vals.append(field)
                    usize, csize, off = vals
                q += 4 + hlen
        members[name] = {'method': method, 'crc': crc, 'csize': csize, 'usize': usize, 'offset': off}
        p += 46 + nlen + xlen + clen
    return members


def read_member(info, chunk=64 * 2**20):
    head = get_range(info['offset'], info['offset'] + 29)
    nlen, xlen = struct.unpack('<2H', head[26:30])
    start = info['offset'] + 30 + nlen + xlen
    end = start + info['csize'] - 1
    out = bytearray()
    d = zlib.decompressobj(-15) if info['method'] == 8 else None
    pos = start
    while pos <= end:
        block = get_range(pos, min(end, pos + chunk - 1))
        out += d.decompress(block) if d else block
        pos += len(block)
    if d:
        out += d.flush()
    data = bytes(out)
    assert len(data) == info['usize'] and zlib.crc32(data) == info['crc'], 'size/CRC mismatch'
    return data
