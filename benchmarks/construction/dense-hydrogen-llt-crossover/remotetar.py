"""Range-request access to members of the remote Zenodo tar files (no full download).

    python3 -B remotetar.py list N_128_T1300K.tar > listing.json
    python3 -B remotetar.py get N_128_T1300K.tar <member> <outfile>

Walks the 512-byte ustar/GNU headers with HTTP range reads; handles GNU long names ('L') and pax ('x') headers.
"""
import json
import sys
import time
import requests

RECORD = 'https://zenodo.org/api/records/14790037/files'


def url(name):
    return f'{RECORD}/{name}/content'


def fetch(name, start, end, tries=12):
    for k in range(tries):
        try:
            r = requests.get(url(name), headers={'Range': f'bytes={start}-{end}'}, timeout=300)
            if r.status_code == 206:
                return r.content
            if r.status_code == 429:   # Zenodo rate limit: wait as instructed
                time.sleep(float(r.headers.get('Retry-After', 60)) + 1)
                continue
            raise RuntimeError(f'HTTP {r.status_code}')
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(5 * (k + 1))
    raise RuntimeError('rate limited')


def octal(field):
    field = field.rstrip(b'\0 ').strip()
    if not field:
        return 0
    if field[0] & 0x80:   # GNU base-256
        return int.from_bytes(field[1:], 'big')
    return int(field, 8)


def listing(name, base=0, end=None, prefetch=262144):
    """List a tar stored in `name` starting at byte `base` (a nested tar member when base > 0)."""
    entries, off, longname = [], base, None
    buf, buf_off = b'', 0
    while True:
        if not (buf_off <= off and off + 512 <= buf_off + len(buf)):
            buf, buf_off = fetch(name, off, off + prefetch - 1), off
        h = buf[off - buf_off:off - buf_off + 512]
        if len(h) < 512 or h == b'\0' * 512 or (end is not None and off >= end):
            break
        size = octal(h[124:136])
        typ = h[156:157]
        fname = h[:100].rstrip(b'\0').decode('utf-8', 'replace')
        prefix = h[345:500].rstrip(b'\0').decode('utf-8', 'replace') if h[257:262] == b'ustar' and h[257:265] != b'ustar  \0' else ''
        if prefix:
            fname = prefix + '/' + fname
        data_off = off + 512
        nxt = data_off + ((size + 511) // 512) * 512
        if typ in (b'L', b'x'):
            body = fetch(name, data_off, data_off + size - 1)
            if typ == b'L':
                longname = body.rstrip(b'\0').decode('utf-8', 'replace')
            else:
                for rec in body.decode('utf-8', 'replace').split('\n'):
                    if ' path=' in rec:
                        longname = rec.split(' path=', 1)[1]
            off = nxt
            continue
        if longname:
            fname, longname = longname, None
        entries.append({'name': fname, 'type': typ.decode() or '0', 'size': size, 'offset': data_off})
        off = nxt
    return entries


def extract(name, entry, start=0, length=None):
    length = entry['size'] - start if length is None else min(length, entry['size'] - start)
    if length <= 0:
        return b''
    out, pos, chunk = [], entry['offset'] + start, 64 * 2**20
    end = entry['offset'] + start + length
    while pos < end:
        stop = min(pos + chunk, end)
        out.append(fetch(name, pos, stop - 1))
        pos = stop
    return b''.join(out)


if __name__ == '__main__':
    cmd, name = sys.argv[1], sys.argv[2]
    if cmd == 'list':
        print(json.dumps(listing(name), indent=0))
    else:
        entry = [e for e in listing(name) if e['name'] == sys.argv[3]][0]
        open(sys.argv[4], 'wb').write(extract(name, entry))
