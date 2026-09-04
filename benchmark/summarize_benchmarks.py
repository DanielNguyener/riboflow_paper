#!/usr/bin/env python3
"""Regenerate the performance table (RiboFlow_v2 computational performance) from the sanitized benchmark evidence shipped in this directory."""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

SCENARIOS = ('one_library', 'three_libraries')
REPLICATES = (1, 2, 3)

EXPECTED_TASKS = {'one_library': 50, 'three_libraries': 125}

CPU_DIVISOR_CANDIDATES = (100.0, 10000.0)

BYTES_PER_GIB = 1024 ** 3

DP_WALL = 0
DP_CPU_H = 3
DP_UTIL = 1

PUBLISHED = {
    'one_library': {'input_reads': 14187683, 'wall_time_s': 516,
                    'cpu_time_h': 0.942, 'cpu_utilization_pct': 41.10,
                    'max_task_rss_gib': 30.02},
    'three_libraries': {'input_reads': 62536429, 'wall_time_s': 1646,
                        'cpu_time_h': 3.878, 'cpu_utilization_pct': 53.00,
                        'max_task_rss_gib': 30.08},
    'change': {'input_reads': 4.41, 'wall_time_s': 3.19, 'cpu_time_h': 4.12,
               'cpu_utilization_pct_pp': 11.90, 'max_task_rss_pct': 0.22},
}

def run_id(scenario, replicate):
    return '%s.rep%d' % (scenario, replicate)

def read_trace(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def load_input_reads(reads_csv):
    """Input fragments per run = ribo + RNA-seq input_fragments.

    reads.csv reports input_fragments once per route (genome and transcriptome)
    for the same FASTQ, so only the genome rows are summed -- adding the
    transcriptome rows would double-count every read.
    """
    total = {}
    with open(reads_csv) as fh:
        for r in csv.DictReader(fh):
            if r['metric'] != 'input_fragments':
                continue
            if r['path'] not in ('ribo_genome', 'rnaseq_genome'):
                continue
            key = (r['scenario'], int(r['replicate']))
            total[key] = total.get(key, 0) + int(round(float(r['value'])))
    return total

def calibrate_cpu_divisor(tasks, wall_s, executor_cpus):
    """Pick the %cpu divisor from the trace alone. See the module docstring."""
    capacity = wall_s * executor_cpus
    viable = []
    for div in CPU_DIVISOR_CANDIDATES:
        cpu_s = sum(float(t['realtime']) / 1000.0 * float(t['%cpu']) / div
                    for t in tasks)
        util = cpu_s / capacity * 100.0 if capacity else 0.0
        if 1.0 < util <= 100.0:
            viable.append((div, cpu_s, util))
    if len(viable) != 1:
        raise SystemExit(
            'ERROR: %d candidate %%cpu divisors are physically plausible '
            '(expected exactly 1): %s. Refusing to guess.'
            % (len(viable), [v[0] for v in viable]))
    return viable[0]

def summarize_run(bench_dir, scenario, replicate, executor_cpus, input_reads):
    rid = run_id(scenario, replicate)
    rdir = os.path.join(bench_dir, 'runs', rid)

    exit_code = int(open(os.path.join(rdir, 'exit_code')).read().strip())
    if exit_code != 0:
        raise SystemExit('ERROR: %s exited %d; refusing to summarize a failed run'
                         % (rid, exit_code))

    wall_s = int(open(os.path.join(rdir, 'wall_seconds')).read().strip())
    manifest = json.load(open(os.path.join(rdir, 'manifest.json')))
    tasks = read_trace(os.path.join(rdir, 'trace.txt'))

    n_completed = sum(1 for t in tasks if t['status'] == 'COMPLETED')
    n_cached = sum(1 for t in tasks if t['status'] == 'CACHED')
    n_retried = sum(1 for t in tasks if int(t['attempt']) > 1)
    n_task_failed = sum(1 for t in tasks if t['exit'] not in ('0', '-'))
    expected = EXPECTED_TASKS[scenario]
    for label, got, want in (('n_tasks', len(tasks), expected),
                             ('n_completed', n_completed, expected),
                             ('n_cached', n_cached, 0),
                             ('n_retried', n_retried, 0),
                             ('n_task_nonzero_exit', n_task_failed, 0)):
        if got != want:
            raise SystemExit('ERROR: %s %s = %d, expected %d' % (rid, label, got, want))

    divisor, cpu_s, util_pct = calibrate_cpu_divisor(tasks, wall_s, executor_cpus)

    max_task_peak_rss_b = max(int(t['peak_rss']) for t in tasks)
    sum_realtime_s = sum(float(t['realtime']) / 1000.0 for t in tasks)

    return {
        'run_id': rid,
        'scenario': scenario,
        'replicate': replicate,
        'description': manifest.get('description', ''),
        'exit_code': exit_code,
        'n_tasks': len(tasks),
        'n_completed': n_completed,
        'n_cached': n_cached,
        'n_retried': n_retried,
        'input_reads': input_reads,
        'wall_time_s': wall_s,
        'cpu_time_s': cpu_s,
        'cpu_time_h': cpu_s / 3600.0,
        'cpu_utilization_pct': util_pct,
        'max_task_peak_rss_bytes': max_task_peak_rss_b,
        'max_task_peak_rss_gib': max_task_peak_rss_b / BYTES_PER_GIB,
        'sum_task_realtime_s': sum_realtime_s,
        'executor_cpus': executor_cpus,
        'cpu_pct_divisor': int(divisor),
    }

def mean(xs):
    return sum(xs) / float(len(xs))

def mean_of_rounded(xs, dp):
    """Convention 3: average the per-run values AFTER rounding to display precision."""
    return mean([round(x, dp) for x in xs])

def summarize_scenario(runs):
    """One performance-table column. Every averaging rule here is deliberate."""
    n = len(runs)
    reads = {r['input_reads'] for r in runs}
    if len(reads) != 1:
        raise SystemExit('ERROR: replicates disagree on input_reads: %s' % sorted(reads))

    rss_bytes = [r['max_task_peak_rss_bytes'] for r in runs]
    mean_rss_b = mean(rss_bytes)

    return {
        'scenario': runs[0]['scenario'],
        'n_repeats': n,
        'input_reads': reads.pop(),
        'wall_time_s': mean_of_rounded([r['wall_time_s'] for r in runs], DP_WALL),
        'cpu_time_h': mean_of_rounded([r['cpu_time_h'] for r in runs], DP_CPU_H),
        'cpu_utilization_pct': mean_of_rounded(
            [r['cpu_utilization_pct'] for r in runs], DP_UTIL),
        'wall_time_s_unrounded': mean([r['wall_time_s'] for r in runs]),
        'cpu_time_h_unrounded': mean([r['cpu_time_h'] for r in runs]),
        'cpu_utilization_pct_unrounded': mean([r['cpu_utilization_pct'] for r in runs]),
        'max_task_rss_bytes_mean': mean_rss_b,
        'max_task_rss_gib': mean_rss_b / BYTES_PER_GIB,
    }

def build_change_row(one, three):
    """The change row. Utilization is PERCENTAGE POINTS, never a percentage."""
    return {
        'scenario': 'change',
        'n_repeats': '',
        'input_reads': 'x%.2f' % (three['input_reads'] / float(one['input_reads'])),
        'wall_time_s': 'x%.2f' % (three['wall_time_s'] / one['wall_time_s']),
        'cpu_time_h': 'x%.2f' % (three['cpu_time_h'] / one['cpu_time_h']),
        'cpu_utilization_pct': '%+.2f pp' % (three['cpu_utilization_pct']
                                             - one['cpu_utilization_pct']),
        'wall_time_s_unrounded': '',
        'cpu_time_h_unrounded': '',
        'cpu_utilization_pct_unrounded': '%+.4f pp' % (
            three['cpu_utilization_pct_unrounded'] - one['cpu_utilization_pct_unrounded']),
        'max_task_rss_bytes_mean': '',
        'max_task_rss_gib': '%+.2f%%' % (
            (three['max_task_rss_bytes_mean'] / one['max_task_rss_bytes_mean'] - 1.0)
            * 100.0),
    }

INDIVIDUAL_COLUMNS = [
    'run_id', 'scenario', 'replicate', 'exit_code',
    'n_tasks', 'n_completed', 'n_cached', 'n_retried',
    'input_reads', 'wall_time_s', 'cpu_time_s', 'cpu_time_h', 'cpu_utilization_pct',
    'max_task_peak_rss_bytes', 'max_task_peak_rss_gib', 'sum_task_realtime_s',
    'executor_cpus', 'cpu_pct_divisor', 'description',
]

SUMMARY_COLUMNS = [
    'scenario', 'n_repeats', 'input_reads', 'wall_time_s', 'cpu_time_h',
    'cpu_utilization_pct', 'max_task_rss_gib',
    'wall_time_s_unrounded', 'cpu_time_h_unrounded', 'cpu_utilization_pct_unrounded',
    'max_task_rss_bytes_mean',
]

FMT = {
    'cpu_time_s': '%.1f', 'cpu_time_h': '%.3f', 'cpu_utilization_pct': '%.2f',
    'max_task_peak_rss_gib': '%.4f', 'sum_task_realtime_s': '%.1f',
    'wall_time_s': '%.0f', 'max_task_rss_gib': '%.2f',
    'wall_time_s_unrounded': '%.4f', 'cpu_time_h_unrounded': '%.6f',
    'cpu_utilization_pct_unrounded': '%.4f', 'max_task_rss_bytes_mean': '%.3f',
}

def render(col, val):
    if isinstance(val, str) or val is None:
        return val
    if col in FMT and isinstance(val, float):
        return FMT[col] % val
    return val

def write_csv(path, columns, rows):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(columns)
        for r in rows:
            w.writerow([render(c, r.get(c, '')) for c in columns])

def check_published(scen_rows, change):
    """Assert every published performance-table value. Returns a list of failures."""
    fails = []

    def cmp(label, got, want, tol):
        if abs(float(got) - float(want)) > tol:
            fails.append('%s: got %s, published %s' % (label, got, want))

    for name, row in scen_rows.items():
        p = PUBLISHED[name]
        cmp('%s input_reads' % name, row['input_reads'], p['input_reads'], 0)
        cmp('%s wall_time_s' % name, row['wall_time_s'], p['wall_time_s'], 0.5)
        cmp('%s cpu_time_h' % name, row['cpu_time_h'], p['cpu_time_h'], 5e-4)
        cmp('%s cpu_utilization_pct' % name, row['cpu_utilization_pct'],
            p['cpu_utilization_pct'], 5e-3)
        cmp('%s max_task_rss_gib' % name, row['max_task_rss_gib'],
            p['max_task_rss_gib'], 5e-3)

    p = PUBLISHED['change']
    cmp('change input_reads', change['input_reads'].lstrip('x'), p['input_reads'], 5e-3)
    cmp('change wall_time_s', change['wall_time_s'].lstrip('x'), p['wall_time_s'], 5e-3)
    cmp('change cpu_time_h', change['cpu_time_h'].lstrip('x'), p['cpu_time_h'], 5e-3)
    cmp('change cpu_utilization_pp',
        change['cpu_utilization_pct'].replace(' pp', ''), p['cpu_utilization_pct_pp'], 5e-3)
    cmp('change max_task_rss_pct',
        change['max_task_rss_gib'].rstrip('%'), p['max_task_rss_pct'], 5e-3)
    return fails

def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bench-dir', default=HERE,
                    help='directory holding runs/, host.json, reads.csv (default: this one)')
    ap.add_argument('--outdir', default=None,
                    help='where to write the two CSVs (default: --bench-dir)')
    ap.add_argument('--check', action='store_true',
                    help='additionally assert every published performance-table value')
    args = ap.parse_args(argv)

    bench = os.path.abspath(args.bench_dir)
    outdir = os.path.abspath(args.outdir or bench)
    os.makedirs(outdir, exist_ok=True)

    host = json.load(open(os.path.join(bench, 'host.json')))
    executor_cpus = int(host['executor_cpus'])
    input_reads = load_input_reads(os.path.join(bench, 'reads.csv'))

    runs = []
    for scenario in SCENARIOS:
        for rep in REPLICATES:
            key = (scenario, rep)
            if key not in input_reads:
                raise SystemExit('ERROR: reads.csv has no input_fragments for %s'
                                 % run_id(scenario, rep))
            runs.append(summarize_run(bench, scenario, rep, executor_cpus,
                                      input_reads[key]))

    divisors = {r['cpu_pct_divisor'] for r in runs}
    if len(divisors) != 1:
        raise SystemExit('ERROR: runs disagree on the %%cpu divisor: %s' % sorted(divisors))

    scen_rows = {s: summarize_scenario([r for r in runs if r['scenario'] == s])
                 for s in SCENARIOS}
    change = build_change_row(scen_rows['one_library'], scen_rows['three_libraries'])

    ind_path = os.path.join(outdir, 'individual_runs.csv')
    sum_path = os.path.join(outdir, 'benchmark_summary.csv')
    write_csv(ind_path, INDIVIDUAL_COLUMNS, runs)
    write_csv(sum_path, SUMMARY_COLUMNS,
              [scen_rows['one_library'], scen_rows['three_libraries'], change])

    print('wrote %s  (%d rows)' % (ind_path, len(runs)))
    print('wrote %s  (3 rows)' % sum_path)
    print()
    print('Performance table  (executor: %d CPU, %s; %%cpu divisor %d)'
          % (executor_cpus, host.get('executor_memory', '?'), divisors.pop()))
    print('%-24s %16s %16s %14s' % ('', 'one_library', 'three_libraries', 'change'))
    o, t = scen_rows['one_library'], scen_rows['three_libraries']
    print('%-24s %16s %16s %14s' % ('input reads', '{:,}'.format(o['input_reads']),
                                    '{:,}'.format(t['input_reads']), change['input_reads']))
    print('%-24s %16.0f %16.0f %14s' % ('wall time (s)', o['wall_time_s'],
                                        t['wall_time_s'], change['wall_time_s']))
    print('%-24s %16.3f %16.3f %14s' % ('CPU time (h)', o['cpu_time_h'],
                                        t['cpu_time_h'], change['cpu_time_h']))
    print('%-24s %16.2f %16.2f %14s' % ('CPU utilization (%)', o['cpu_utilization_pct'],
                                        t['cpu_utilization_pct'],
                                        change['cpu_utilization_pct']))
    print('%-24s %16.2f %16.2f %14s' % ('max task RSS (GiB)', o['max_task_rss_gib'],
                                        t['max_task_rss_gib'], change['max_task_rss_gib']))

    if args.check:
        fails = check_published(scen_rows, change)
        print()
        if fails:
            print('CHECK FAILED (%d):' % len(fails))
            for f in fails:
                print('  ' + f)
            return 1
        print('CHECK OK: all 14 published performance-table values reproduce exactly.')
    return 0

if __name__ == '__main__':
    sys.exit(main())
