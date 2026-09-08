"""Build serial Qdyn/Qprep from a pinned Git archive in a fresh directory.

Retain source, compiler identity, build command/log and binary hashes. This is
reproducibility evidence, not a signed attestation or numerical qualification.
"""
from __future__ import annotations

import argparse
import io
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import signal
import subprocess
import tarfile

from .charge_protocol import fingerprint


def _local(root, name):
    path = PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts:
        raise ValueError('Unsafe build artifact path')
    resolved = (root/name).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError('Build artifact escaped its directory')
    return resolved


def build(repository, directory, compiler, commit='HEAD'):
    resolved = subprocess.run(['git', 'rev-parse', '--verify', '--end-of-options', commit+'^{commit}'],
                              cwd=repository, capture_output=True, text=True, check=True, timeout=10).stdout.strip()
    archive = subprocess.run(['git', 'archive', '--format=tar', resolved, 'src/q6'], cwd=repository,
                             capture_output=True, check=True, timeout=10).stdout
    return _build_archive(archive, directory, compiler, resolved)


def build_archive(path, directory, compiler, *, expected_commit, expected_sha256):
    """Build a transferred Git archive without requiring a remote Git checkout."""
    archive = path.read_bytes()
    if hashlib.sha256(archive).hexdigest() != expected_sha256:
        raise ValueError('Transferred native archive hash mismatch')
    return _build_archive(archive, directory, compiler, expected_commit)


def _build_archive(archive, directory, compiler, resolved):
    if not re.fullmatch('[0-9a-f]{40}', resolved):
        raise ValueError('Require a full source commit')
    with tarfile.open(fileobj=io.BytesIO(archive), mode='r:') as source:
        if source.pax_headers.get('comment') != resolved:
            raise ValueError('Source commit differs from Git archive metadata')
    compiler_path = shutil.which(str(compiler))
    make_path = shutil.which('make')
    if compiler_path is None or make_path is None:
        raise ValueError('Compiler and make must already be installed')
    compiler_path = str(Path(compiler_path).resolve())
    if not re.fullmatch(r'[A-Za-z0-9_./+@\-]+', compiler_path):
        raise ValueError('Compiler path contains unsupported make/shell characters')
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    (directory/'native-source.tar').write_bytes(archive)
    source_files = {}
    with tarfile.open(fileobj=io.BytesIO(archive), mode='r:') as source:
        for member in source.getmembers():
            target = _local(directory, member.name)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            if not member.isfile() or not member.name.startswith('src/q6/'):
                raise ValueError('Archive contains unsupported native source member')
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open('xb') as stream:
                stream.write(source.extractfile(member).read())
            source_files[member.name] = fingerprint(target)
    if any(name.endswith(('.o', '.mod', '/qdyn', '/qprep')) for name in source_files):
        raise ValueError('Source archive contains prebuilt native artifacts')
    version = subprocess.run([compiler_path, '--version'], capture_output=True, text=True,
                             check=True, timeout=10).stdout
    compiler_hash, make_hash = fingerprint(Path(compiler_path)), fingerprint(Path(make_path))
    command = [make_path, '-j1', '-C', str(directory/'src/q6'), 'qdyn', 'qprep', 'FC='+compiler_path]
    # Do not inherit MAKEFLAGS, FC/FFLAGS or library injection environment settings.
    environment = {'PATH': os.environ.get('PATH', '/usr/bin:/bin'), 'LC_ALL': 'C'}
    if 'TMPDIR' in os.environ:
        environment['TMPDIR'] = os.environ['TMPDIR']
    with (directory/'build.log').open('x') as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   env=environment, start_new_session=True)
        try:
            returncode = process.wait(timeout=60)
        except subprocess.TimeoutExpired:
            # The build owns this new process group, including compiler children.
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            raise
    if returncode:
        raise ValueError(f'Native build failed ({returncode}); retain {directory}/build.log')
    if fingerprint(Path(compiler_path)) != compiler_hash or fingerprint(Path(make_path)) != make_hash:
        raise ValueError('Compiler or make executable changed during build')
    if any(fingerprint(_local(directory, name)) != checksum for name, checksum in source_files.items()):
        raise ValueError('Native source changed during build')
    report = {'schema_version': 1, 'gate': 'isolated_native_build_recorded', 'production_ready': False,
              'source_commit': resolved, 'source_archive_sha256': fingerprint(directory/'native-source.tar'),
              'source_files_sha256': source_files, 'compiler': {'binary': compiler_path,
              'sha256': compiler_hash, 'version': version},
              'command': command, 'make_sha256': make_hash, 'environment': environment,
              'platform': platform.platform(),
              'build_log_sha256': fingerprint(directory/'build.log'),
              'binaries': {name: {'path': f'src/q6/{name}', 'sha256': fingerprint(directory/f'src/q6/{name}')}
                           for name in ('qdyn', 'qprep')}}
    with (directory/'build.json').open('x') as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    return report


def validate(path):
    report = json.loads(path.read_text())
    root = path.parent
    if report['schema_version'] != 1 or report['gate'] != 'isolated_native_build_recorded':
        raise ValueError('Unsupported native build report')
    if fingerprint(root/'native-source.tar') != report['source_archive_sha256']:
        raise ValueError('Native source archive changed')
    with tarfile.open(root/'native-source.tar', mode='r:') as source:
        if source.pax_headers.get('comment') != report['source_commit']:
            raise ValueError('Source commit differs from Git archive metadata')
        archived = {}
        for member in source.getmembers():
            _local(root, member.name)
            if member.isdir():
                continue
            if not member.isfile() or not member.name.startswith('src/q6/') or member.name in archived:
                raise ValueError('Unsupported or duplicate archive source member')
            archived[member.name] = hashlib.sha256(source.extractfile(member).read()).hexdigest()
        if archived != report['source_files_sha256']:
            raise ValueError('Recorded source inventory differs from archive contents')
    if fingerprint(root/'build.log') != report['build_log_sha256']:
        raise ValueError('Native build log changed')
    for name, checksum in report['source_files_sha256'].items():
        if fingerprint(_local(root, name)) != checksum:
            raise ValueError('Native source snapshot changed')
    for binary in report['binaries'].values():
        if fingerprint(_local(root, binary['path'])) != binary['sha256']:
            raise ValueError('Built native executable changed')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repository', type=Path)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--commit', default='HEAD')
    args = parser.parse_args()
    try:
        report = build(**vars(args))
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        parser.exit(2, f'Charge build failed: {error}\n')
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
