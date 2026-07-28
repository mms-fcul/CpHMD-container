#!/usr/bin/env python3
"""Audit GROMACS performance records from completed CpHMD run directories."""

from __future__ import annotations

import argparse
import gzip
import math
import re
import statistics
import sys
from pathlib import Path

PERF_RE = re.compile(
    r"^\s*Performance:\s+([0-9.+\-Ee]+)\s+([0-9.+\-Ee]+)\s*$"
)
TIME_RE = re.compile(
    r"^\s*Time:\s+([0-9.+\-Ee]+)\s+([0-9.+\-Ee]+)\s+([0-9.+\-Ee]+)\s*$"
)
HOST_RE = re.compile(r"^bio\d+$")
GPU_USE_PATTERNS = (
    "Using 1 GPU",
    "Using 1 compatible GPU",
    "PP tasks will do",
    "PME tasks will do",
    "GPU activities",
    "GPU timings",
    "CUDA",
)


def read_text(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", errors="replace") as handle:
            return handle.read()
    return path.read_text(errors="replace")


def find_one(run: Path, patterns: list[str]) -> Path | None:
    found: list[Path] = []
    for pattern in patterns:
        found.extend(run.glob(pattern))
    found = sorted(set(found))
    if not found:
        return None
    if len(found) > 1:
        print(
            f"WARNING: multiple files matched in {run}; using {found[-1].name}",
            file=sys.stderr,
        )
    return found[-1]


def parse_records(text: str) -> list[dict[str, float]]:
    records: list[dict[str, float]] = []
    pending_time: tuple[float, float, float] | None = None

    for line in text.splitlines():
        match = TIME_RE.match(line)
        if match:
            pending_time = tuple(float(value) for value in match.groups())
            continue

        match = PERF_RE.match(line)
        if not match:
            continue

        ns_day = float(match.group(1))
        hour_ns = float(match.group(2))
        record: dict[str, float] = {
            "ns_day": ns_day,
            "hour_ns": hour_ns,
        }

        if pending_time is not None:
            core_s, wall_s, percent = pending_time
            record.update(
                core_s=core_s,
                wall_s=wall_s,
                parallel_percent=percent,
            )
            # Recover the amount of simulated time represented by this record.
            record["simulated_ns"] = ns_day * wall_s / 86400.0

        records.append(record)
        pending_time = None

    return records


def detect_host(run: Path) -> str:
    serial = find_one(run, ["serial_test_*.log", "serial_*.log", "slurm-*.out"])
    if serial is None:
        return "unknown"
    for line in read_text(serial).splitlines():
        stripped = line.strip()
        if HOST_RE.match(stripped):
            return stripped
    return "unknown"


def summarize(label: str, run: Path) -> dict[str, object]:
    log = find_one(run, ["*_CpHrun_cycles*.log.gz", "*_CpHrun_cycles*.log"])
    if log is None:
        raise RuntimeError(f"No finalized GROMACS segment log found in {run}")

    text = read_text(log)
    records = parse_records(text)
    if not records:
        raise RuntimeError(f"No GROMACS Performance records found in {log}")

    ns_days = [float(record["ns_day"]) for record in records]
    wall_records = [record for record in records if "wall_s" in record]

    aggregate_ns_day = math.nan
    total_wall_s = math.nan
    total_simulated_ns = math.nan
    if wall_records:
        total_wall_s = sum(float(record["wall_s"]) for record in wall_records)
        total_simulated_ns = sum(
            float(record["simulated_ns"]) for record in wall_records
        )
        if total_wall_s > 0:
            aggregate_ns_day = total_simulated_ns / total_wall_s * 86400.0

    gpu_markers = [
        pattern for pattern in GPU_USE_PATTERNS if pattern.lower() in text.lower()
    ]

    return {
        "label": label,
        "run": run,
        "host": detect_host(run),
        "log": log,
        "records": records,
        "mean_ns_day": statistics.fmean(ns_days),
        "median_ns_day": statistics.median(ns_days),
        "min_ns_day": min(ns_days),
        "max_ns_day": max(ns_days),
        "aggregate_ns_day": aggregate_ns_day,
        "total_wall_s": total_wall_s,
        "total_simulated_ns": total_simulated_ns,
        "gpu_markers": gpu_markers,
    }


def fmt(value: object, digits: int = 3) -> str:
    if not isinstance(value, (int, float)) or math.isnan(float(value)):
        return "n/a"
    return f"{float(value):.{digits}f}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cpu_run", type=Path)
    parser.add_argument("gpu_run", type=Path)
    args = parser.parse_args()

    for run in (args.cpu_run, args.gpu_run):
        if not run.is_dir():
            parser.error(f"run directory not found: {run}")

    summaries = [
        summarize("CPU", args.cpu_run.resolve()),
        summarize("GPU", args.gpu_run.resolve()),
    ]

    print(
        "mode\thost\trecords\taggregate_ns/day\tmedian_ns/day\t"
        "min_ns/day\tmax_ns/day\tMD_wall_s\tlog"
    )
    for item in summaries:
        print(
            f"{item['label']}\t{item['host']}\t{len(item['records'])}\t"
            f"{fmt(item['aggregate_ns_day'])}\t{fmt(item['median_ns_day'])}\t"
            f"{fmt(item['min_ns_day'])}\t{fmt(item['max_ns_day'])}\t"
            f"{fmt(item['total_wall_s'])}\t{item['log']}"
        )

    cpu, gpu = summaries
    if (
        isinstance(cpu["aggregate_ns_day"], float)
        and isinstance(gpu["aggregate_ns_day"], float)
        and not math.isnan(cpu["aggregate_ns_day"])
        and not math.isnan(gpu["aggregate_ns_day"])
        and cpu["aggregate_ns_day"] > 0
    ):
        speedup = gpu["aggregate_ns_day"] / cpu["aggregate_ns_day"]
        print(f"\nObserved GPU/CPU raw-MD throughput ratio: {speedup:.3f}x")
    else:
        print("\nObserved GPU/CPU raw-MD throughput ratio: unavailable")

    if cpu["host"] != gpu["host"]:
        print(
            "CAUTION: CPU and GPU jobs ran on different hosts; the ratio is not "
            "a controlled hardware comparison."
        )

    print(
        "CAUTION: these CpHMD segments contain only three 20-ps effective-MD "
        "cycles. Startup/autotuning and PB/MC overhead make them a smoke test, "
        "not a final benchmark."
    )

    print("\nPer-record GROMACS performance:")
    for item in summaries:
        for index, record in enumerate(item["records"], start=1):
            print(
                f"{item['label']} record {index}: "
                f"{record['ns_day']:.3f} ns/day, "
                f"{record['hour_ns']:.3f} hour/ns, "
                f"wall={fmt(record.get('wall_s'))} s"
            )

    print("\nGPU-offload markers found:")
    if gpu["gpu_markers"]:
        for marker in gpu["gpu_markers"]:
            print(f"  {marker}")
    else:
        print("  none; inspect the GPU log before accepting the benchmark")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
