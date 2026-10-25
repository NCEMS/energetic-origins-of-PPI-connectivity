#!/usr/bin/env python3
"""Stage selected PDBs and generate one portable FastRelax command per replicate.

Generation requires pandas and numpy. Execution on the cluster needs only
Python 3 and Rosetta. Relative source PDB paths are resolved from --source-root
(default: the current working directory), as in the original runner.
"""
import argparse
import csv
import hashlib
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def run_pose(args):
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    stem = f'{args.label}_{args.pose:04d}'
    final_pdb = output / f'{stem}.pdb'
    final_sc = output / f'{stem}.sc'
    if final_pdb.is_file() and final_sc.is_file():
        print(f'Skipping {stem}: both outputs exist')
        return
    scratch_root = output / '.scratch'
    scratch_root.mkdir(exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix=f'{stem}_', dir=scratch_root))
    command = [
        args.relax_executable, '-s', str(Path(args.pdb).resolve()),
        '-relax:fast', '-relax:constrain_relax_to_start_coords',
        '-nstruct', '1', '-score:weights', 'ref2015',
        '-out:path:all', str(scratch.resolve()),
        '-out:file:scorefile', str((scratch / f'{stem}.sc').resolve()),
        '-out:pdb', 'true',
    ]
    try:
        with (output / f'{stem}.log').open('w') as log:
            log.write('=== CMD ===\n' + shlex.join(command) + '\n\n')
            log.flush()
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        pdbs = sorted(scratch.glob('*.pdb'))
        scores = sorted(scratch.glob('*.sc'))
        if len(pdbs) != 1 or len(scores) != 1:
            raise RuntimeError(f'{stem}: expected one PDB and one scorefile; found {len(pdbs)} and {len(scores)}')
        os.replace(pdbs[0], final_pdb)
        os.replace(scores[0], final_sc)
    except Exception:
        print(f'{stem} failed; see {output / (stem + ".log")}; scratch retained at {scratch}', file=sys.stderr)
        raise
    else:
        shutil.rmtree(scratch)
        print(f'{stem}: OK')


def generate(args):
    import numpy as np
    import pandas as pd

    if args.nstruct < 1:
        raise ValueError('--nstruct must be at least 1')
    remote = Path(args.remote_bundle)
    if not remote.is_absolute():
        raise ValueError('--remote-bundle must be an absolute cluster path')
    bundle = Path(args.bundle).resolve()
    # Avoid leaving stale command lists or mixing different preparation runs.
    if bundle.exists() and any(bundle.iterdir()):
        raise ValueError(f'Bundle directory must be new or empty: {bundle}')
    source_root = Path(args.source_root).resolve()
    required = {'node', 'has_verified_sequence', 'DeepTMHMM_class',
                'final_structure_path', 'final_structure_source'}
    selected = {}
    report = []
    for filename in args.nodes:
        filename = Path(filename)
        frame = (pd.read_pickle(filename) if filename.suffix.lower() in {'.pkl', '.pickle'}
                 else pd.read_csv(filename))
        frame.columns = frame.columns.str.strip()
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f'{filename}: missing columns {sorted(missing)}. Use the structure-annotated node file, not a node-only list.')
        counts = {'rows': len(frame), 'filtered': 0, 'eligible': 0, 'duplicate': 0}
        for _, row in frame.iterrows():
            verified = row['has_verified_sequence']
            verified = isinstance(verified, (bool, np.bool_)) and bool(verified)
            source = row['final_structure_source']
            path_value = row['final_structure_path']
            if not (verified and row['DeepTMHMM_class'] in ['GLOB', 'SP']
                    and isinstance(path_value, str) and path_value.strip()
                    and pd.notna(source) and str(source).strip() not in {'', 'None'}):
                counts['filtered'] += 1
                continue
            if pd.isna(row['node']):
                raise ValueError(f'{filename}: eligible row has no node identifier')
            label = str(row['node']).strip()
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', label):
                raise ValueError(f'{filename}: unsafe or grouped node label: {label!r}')
            path = Path(path_value).expanduser()
            path = (path if path.is_absolute() else source_root / path).resolve()
            if not path.is_file():
                raise FileNotFoundError(f'{label}: selected PDB not found: {path}')
            sha = digest(path)
            counts['eligible'] += 1
            if label in selected:
                if selected[label]['sha256'] != sha:
                    raise ValueError(f'Conflicting structures for node {label}: {selected[label]["source_pdb"]} and {path}')
                selected[label]['node_files'].append(str(filename))
                counts['duplicate'] += 1
            else:
                selected[label] = {'node': label, 'source_pdb': str(path),
                                   'sha256': sha, 'node_files': [str(filename)]}
        report.append((str(filename), counts))
    if not selected:
        raise ValueError('No eligible structures selected; no bundle created.')
    inputs = bundle / 'inputs'
    inputs.mkdir(parents=True, exist_ok=True)
    script_name = 'make-rosetta-commands.py'
    shutil.copy2(Path(__file__).resolve(), bundle / script_name)
    commands = []
    with (bundle / 'selected_structures.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['node', 'source_pdb', 'sha256', 'node_files'], delimiter='\t')
        writer.writeheader()
        for label, record in sorted(selected.items()):
            staged = inputs / f'{label}.pdb'
            shutil.copy2(record['source_pdb'], staged)
            if digest(staged) != record['sha256']:
                raise RuntimeError(f'Source changed during staging: {label}')
            writer.writerow({**record, 'node_files': ';'.join(record['node_files'])})
            for pose in range(1, args.nstruct + 1):
                commands.append(shlex.join([
                    args.python, str(remote / script_name), 'run-pose',
                    '--pdb', str(remote / 'inputs' / f'{label}.pdb'),
                    '--label', label, '--pose', str(pose),
                    '--output-dir', str(remote / 'scores'),
                    '--relax-executable', args.relax_executable,
                ]))
    (bundle / 'commands.txt').write_text('\n'.join(commands) + '\n')
    for filename, counts in report:
        print(f'{filename}: ' + ', '.join(f'{key}={value:,}' for key, value in counts.items()))
    print(f'Unique selected proteins: {len(selected):,}')
    print(f'Replicates per protein: {args.nstruct:,}')
    print(f'Commands: {len(commands):,}')
    print(f'Bundle to transfer: {bundle}')
    print(f'Expected cluster location: {remote}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    gen = sub.add_parser('generate', help='Read node files, copy PDBs, and write commands')
    gen.add_argument('--nodes', nargs='+', required=True)
    gen.add_argument('--bundle', required=True, help='New or empty local bundle directory')
    gen.add_argument('--remote-bundle', required=True, help='Absolute cluster destination of the bundle')
    gen.add_argument('--source-root', default='.', help='Base for relative final_structure_path values')
    gen.add_argument('--relax-executable', required=True, help='Rosetta executable path on the cluster')
    gen.add_argument('--nstruct', type=int, default=10)
    gen.add_argument('--python', default='python3', help='Python executable on the cluster')
    run = sub.add_parser('run-pose', help='Execute one staged replicate on the cluster')
    run.add_argument('--pdb', required=True)
    run.add_argument('--label', required=True)
    run.add_argument('--pose', type=int, required=True)
    run.add_argument('--output-dir', required=True)
    run.add_argument('--relax-executable', required=True)
    args = parser.parse_args()
    if args.mode == 'generate':
        generate(args)
    else:
        run_pose(args)


if __name__ == '__main__':
    main()
