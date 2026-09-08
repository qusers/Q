"""Export a clean, committed source snapshot for the bounded cluster smoke job."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def export(repository, destination):
    dirty = subprocess.run(['git', 'status', '--porcelain'], cwd=repository,
                           capture_output=True, text=True, check=True, timeout=10).stdout
    if dirty:
        raise ValueError('Commit the release implementation before export')
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=repository,
                            capture_output=True, text=True, check=True, timeout=10).stdout.strip()
    destination.mkdir(parents=True, exist_ok=False)
    for name, paths in [
            ('source.tar', ['src/q6', 'src/QligFEP', 'experiments/analytical-charge-validation']),
            ('native-source.tar', ['src/q6'])]:
        subprocess.run(['git', 'archive', '--format=tar', '--output='+str((destination/name).resolve()),
                        commit, *paths], cwd=repository, check=True, timeout=30)
    manifest = {'commit': commit,
                'source_archive_sha256': hashlib.sha256((destination/'source.tar').read_bytes()).hexdigest(),
                'native_archive_sha256': hashlib.sha256((destination/'native-source.tar').read_bytes()).hexdigest(),
                'purpose': 'bounded_snellius_target_smoke', 'production_ready': False}
    with (destination/'release.json').open('x') as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    export(Path(__file__).resolve().parents[2], args.destination.resolve())
