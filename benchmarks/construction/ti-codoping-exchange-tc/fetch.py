"""Download the pinned archive files of Materials Cloud record 6txnj-70460 (DOI 10.24435/materialscloud:b7-6k).

    python3 fetch.py [CACHE]     (default ~/.cache/ai-theorist/ti-codoping-exchange-tc)
Only small files are fetched through the per-file API links; the AiiDA exports (0.4-46 GB) are not needed.
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

RECORD = 'https://archive.materialscloud.org/api/records/6txnj-70460/files/{}/content'
FILES = {  # name: (bytes, md5) as listed by the record API
    'README.md': (6985, '7b0d99083467b6cc0a2215d8b07873d7'),
    'Jij_table.csv': (422637, 'e73c49c9da7b0653f2de3e0e4466eaa4'),
    'Bi2Te3.cif': (1539, 'f714130613532a4063571c361bc6dc46'),
    'Tcs_Bi2Te3.txt': (10047, '1c8a734b01a5b6c2aa18e95a0e2ddc49'),
    'data_analysis_Bi2Te3.ipynb': (1505163, '690b41fab7da399fd6828b77b36e327a'),
    'files_description.md': (1326, '9ce97c9529b32ed4c75cdf4a6f49ee08'),
}


def main():
    cache = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home()/'.cache/ai-theorist/ti-codoping-exchange-tc'
    cache.mkdir(parents=True, exist_ok=True)
    for name, (size, md5) in FILES.items():
        path = cache/name
        if not path.exists() or hashlib.md5(path.read_bytes()).hexdigest() != md5:
            with urllib.request.urlopen(RECORD.format(name), timeout=300) as r:
                path.write_bytes(r.read())
        data = path.read_bytes()
        assert len(data) == size and hashlib.md5(data).hexdigest() == md5, f'checksum mismatch for {name}'
        print(name, size, md5, 'sha256', hashlib.sha256(data).hexdigest())


if __name__ == '__main__':
    main()
