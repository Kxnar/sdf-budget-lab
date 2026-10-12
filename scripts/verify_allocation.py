"""Verify a frozen allocation run directory or evidence ZIP without regenerating it."""
import argparse
import hashlib
import json
import zipfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('evidence', type=Path)
    args = parser.parse_args()
    archive = zipfile.ZipFile(args.evidence) if args.evidence.is_file() else None
    def read(name):
        return archive.read(name) if archive else (args.evidence/name).read_bytes()
    manifest = json.loads(read('manifest.json'))
    for name, expected in manifest['artifact_sha256'].items():
        if hashlib.sha256(read(name)).hexdigest() != expected:
            raise ValueError(f'Artifact hash mismatch: {name}')
    for name, expected in manifest['source_sha256'].items():
        if hashlib.sha256(read('source/'+name)).hexdigest() != expected:
            raise ValueError(f'Source hash mismatch: {name}')
    if archive:
        archive.close()
    print(f'Verified {len(manifest["artifact_sha256"])} artifacts and {len(manifest["source_sha256"])} source snapshots')


if __name__ == '__main__':
    main()
