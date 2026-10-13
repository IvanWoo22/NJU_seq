#!/usr/bin/env python3
"""Run a command with Linux disk, memory and time limits.

Monitor regular files under --output without following symlinks, and stop the
process group if a limit is exceeded.
"""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import time


def size(path):
    total, largest = 0, 0
    if not path.exists():
        return total, largest
    for folder, dirs, files in os.walk(path, followlinks=False):
        for name in files:
            p = Path(folder)/name
            try:
                if p.is_symlink() or not p.is_file():
                    continue
                n = p.stat().st_size
            except FileNotFoundError:
                continue
            total += n
            largest = max(n, largest)
    return total, largest


def memory(group):
    """Linux process-group RSS/swap; RSS sums may double-count shared pages."""
    rss = swap = processes = 0
    for folder in Path('/proc').iterdir():
        if not folder.name.isdigit():
            continue
        try:
            fields = (folder/'stat').read_text().rsplit(')', 1)[1].split()
            if int(fields[2]) != group:
                continue
            status = dict(line.split(':', 1) for line in (folder/'status').read_text().splitlines()
                          if ':' in line)
            rss += int(status.get('VmRSS', '0 kB').split()[0]) * 1024
            swap += int(status.get('VmSwap', '0 kB').split()[0]) * 1024
            processes += 1
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    info = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    available = int(info['MemAvailable'].split()[0]) * 1024
    return rss, swap, available, processes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-gib', type=float, default=4)
    parser.add_argument('--max-file-gib', type=float, default=1)
    parser.add_argument('--min-free-gib', type=float, default=20)
    parser.add_argument('--max-address-gib', type=float, default=0,
                        help='Per-process virtual address-space hard limit; 0 disables')
    parser.add_argument('--max-rss-gib', type=float, default=0,
                        help='Sampled process-group RSS limit; 0 disables')
    parser.add_argument('--max-swap-gib', type=float, default=0,
                        help='Sampled process-group swap limit; 0 disables')
    parser.add_argument('--min-available-gib', type=float, default=0,
                        help='Minimum host MemAvailable; 0 disables')
    parser.add_argument('--timeout', type=float, default=1800)
    parser.add_argument('--poll', type=float, default=2)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command[:1] == ['--']:
        args.command = args.command[1:]
    for key in ('max_gib', 'max_file_gib', 'timeout', 'poll'):
        value = getattr(args, key)
        if not math.isfinite(value) or value <= 0:
            parser.error(f'{key} must be finite and positive')
    for key in ('max_address_gib', 'max_rss_gib', 'max_swap_gib', 'min_available_gib'):
        if not math.isfinite(getattr(args, key)) or getattr(args, key) < 0:
            parser.error(f'{key} must be finite and nonnegative')
    if not math.isfinite(args.min_free_gib) or args.min_free_gib < 0 or not args.command:
        parser.error('Require command and finite nonnegative min-free-gib')
    out = args.output.absolute()
    report = out.with_name(out.name + '.resources.json')
    if out.exists() or out.is_symlink() or report.exists() or report.is_symlink():
        parser.error('Use a new output directory and resource report name')
    if shutil.disk_usage(out.parent).free < args.min_free_gib * 1024**3:
        parser.error('Insufficient free space before launch')
    started = time.monotonic()
    def interrupted(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGHUP, interrupted)
    data = dict(argv=args.command, output=str(out), limits={k:v for k,v in vars(args).items()
                if k not in ('command', 'output')}, peak_output_bytes=0, largest_file_bytes=0,
                peak_rss_bytes=0, peak_swap_bytes=0,
                status='running', monitor_scope='regular files under output only; symlinks excluded; '
                'Linux process-group RSS/swap sampled; summed RSS may count shared pages twice')
    def save():
        temporary = report.with_suffix(report.suffix+'.partial')
        temporary.write_text(json.dumps(data, indent=2)+'\n')
        temporary.replace(report)
    def limits():
        n = int(args.max_file_gib * 1024**3)
        resource.setrlimit(resource.RLIMIT_FSIZE, (n, n))
        if args.max_address_gib:
            n = int(args.max_address_gib * 1024**3)
            resource.setrlimit(resource.RLIMIT_AS, (n, n))
    save()
    child = None
    reason = None
    try:
        child = subprocess.Popen(args.command, start_new_session=True, preexec_fn=limits)
        while True:
            rss, swap, available, processes = memory(child.pid)
            data.update(rss_bytes=rss, swap_bytes=swap, available_memory_bytes=available,
                        process_count=processes, process_group=child.pid)
            data['peak_rss_bytes'] = max(rss, data['peak_rss_bytes'])
            data['peak_swap_bytes'] = max(swap, data['peak_swap_bytes'])
            n, largest = size(out)
            data['peak_output_bytes'] = max(n, data['peak_output_bytes'])
            data['largest_file_bytes'] = max(largest, data['largest_file_bytes'])
            data['elapsed_seconds'] = time.monotonic() - started
            data['free_bytes'] = shutil.disk_usage(out.parent).free
            if n > args.max_gib * 1024**3:
                reason = 'output_size_limit'
            elif data['free_bytes'] < args.min_free_gib * 1024**3:
                reason = 'free_space_limit'
            elif data['elapsed_seconds'] > args.timeout:
                reason = 'wall_time_limit'
            elif args.max_rss_gib and rss > args.max_rss_gib * 1024**3:
                reason = 'rss_limit'
            elif args.max_swap_gib and swap > args.max_swap_gib * 1024**3:
                reason = 'swap_limit'
            elif args.min_available_gib and available < args.min_available_gib * 1024**3:
                reason = 'available_memory_limit'
            code = child.poll()
            if reason or code is not None:
                break
            save()
            time.sleep(args.poll)
    except (KeyboardInterrupt, SystemExit):
        reason = 'interrupted'
    except Exception as e:
        reason = f'launch_or_monitor_error: {e}'
    finally:
        if child is not None:
            if reason or child.poll() not in (None, 0):
                try:
                    os.killpg(child.pid, signal.SIGTERM)
                    child.wait(timeout=5)
                except ProcessLookupError:
                    pass
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                # Descendants may survive a parent that exits promptly on TERM.
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            data['returncode'] = child.wait()
        else:
            data['returncode'] = None
        data['status'] = 'success' if not reason and data['returncode'] == 0 else 'failed'
        data['reason'] = reason or ('command_failed' if data['returncode'] else None)
        data['elapsed_seconds'] = time.monotonic() - started
        save()
        if data['status'] == 'failed' and out.is_dir():
            (out/'GUARD_FAILED').write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps(data), flush=True)
    return 0 if data['status'] == 'success' else 1


if __name__ == '__main__':
    raise SystemExit(main())
