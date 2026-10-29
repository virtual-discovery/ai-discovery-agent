"""Refetch the pinned archive file and unpack the members used.

    python3 -B fetch.py [CACHE]     (default ~/.cache/ai-theorist/mnte-altermagnet-chiral-magnons)
"""
import hashlib
import os
import sys
import urllib.request
import zipfile

URL = 'https://zenodo.org/api/records/20049187/files/Figures.zip/content'
MD5 = 'fb87bc34319008f92adee9d7dbe3ada8'
SIZE = 5319363
MEMBERS = ['Fig1c.mat', 'Figure2/Fig2b_gray.dat', 'Figure2/Fig2b_green.dat', 'Figure2/Fig2c.dat', 'Figure2/Fig2e.dat',
           'Figure2/Fig2g.dat', 'Figure2/Fig2h.dat', 'Figure2/Fig2i.dat', 'Figure2/Fig2j.dat', 'Figure2/Fig2df.mat',
           'Figure2/Fig2ghij.mat', 'Figure3/Fig3b_red', 'Figure3/Fig3b_blue']


def main():
    cache = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/.cache/ai-theorist/mnte-altermagnet-chiral-magnons')
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, 'Figures.zip')
    if not os.path.exists(path) or os.path.getsize(path) != SIZE:
        urllib.request.urlretrieve(URL, path)
    md5 = hashlib.md5(open(path, 'rb').read()).hexdigest()
    assert md5 == MD5, md5
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        wanted = [n for n in names if n in MEMBERS or n.startswith('Figure3/Fig3c/') or n.startswith('Figure3/Fig3d/')]
        z.extractall(os.path.join(cache, 'x'), members=wanted)
    print('ok', md5, len(wanted), 'members')


if __name__ == '__main__':
    main()
