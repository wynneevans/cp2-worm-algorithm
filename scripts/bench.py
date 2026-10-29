#!/usr/bin/env python3
"""Time the simulator on the fixed-seed runs in tests/bench/bench.toml.

    scripts/bench.py --save before-change              # baseline, saved in tests/bench/results/
    scripts/bench.py --compare tests/bench/results/<file>.json --save after-change
    scripts/bench.py --only long --repeats 5           # just the 60 s run
    scripts/bench.py --list                            # the saved results

Each repeat runs every selected run once (rounds are interleaved, so slow drifts in the machine
spread over all runs), into a temporary folder that is deleted afterwards. The time is wall_seconds
from the run's run_<time>.toml, which covers the simulation and its file output but not start-up.
Prints the median and the spread (max - min) per run. Run on an otherwise idle machine; a change
counts only if it is larger than the spread.

--save LABEL writes <UTC time>_<LABEL>.json to the results folder with the times, every run's full
parameters (from its run record, including the 48 weights), the git version, CPU, OS and build.
--compare warns when the parameters, CPU or compiler differ from the baseline's, since the times
are then not comparable.
"""

import argparse
import ast
import json
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "tests/bench/results"
NAME = re.compile(r"^name\s*=\s*\"([^\"]+)\"", re.MULTILINE)


def run_names(config):
    return NAME.findall(config.read_text())


def read_record(text):
    """The simple TOML of a run record as {section: {key: value}} ('run', 'run.probabilities', ...)."""
    sections, current = {}, None
    for line in text.splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        if line.startswith("["):
            current = sections.setdefault(line.strip("[]"), {})
            continue
        key, _, value = line.partition("=")
        value = value.strip().replace("true", "True").replace("false", "False")
        current[key.strip()] = ast.literal_eval(value)
    return sections


def parameters(record):
    """Everything that determines the simulation: the [run] settings and the probability row."""
    p = {k: v for k, v in record["run"].items() if k not in ("name", "data_dir")}
    p["probabilities"] = record["run.probabilities"]
    return p


def time_run(binary, config, name, work):
    """(wall_seconds, record) of one run of `name`, its output written to (and removed from) `work`."""
    out = work / name
    proc = subprocess.run(
        [str(binary), str(config), "--only", name, "--data-dir", str(out)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        sys.exit(f"run '{name}' failed:\n{proc.stdout}{proc.stderr}")
    records = list(out.glob("*/run_*.toml"))
    if len(records) != 1:
        sys.exit(f"run '{name}': expected one run record in {out}, found {len(records)}")
    record = read_record(records[0].read_text())
    shutil.rmtree(out)
    return float(record["provenance"]["wall_seconds"]), record


def command(*args):
    try:
        return subprocess.run(args, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def environment(binary):
    """git version, machine and build of this benchmark."""
    if sys.platform == "darwin":
        cpu = command("sysctl", "-n", "machdep.cpu.brand_string")
        os_name = f"macOS {platform.mac_ver()[0]} (Darwin {platform.release()})"
    else:
        m = re.search(r"^model name\s*:\s*(.+)$", Path("/proc/cpuinfo").read_text(), re.MULTILINE)
        cpu = m.group(1) if m else platform.processor() or "unknown"
        os_name = f"{platform.system()} {platform.release()}"
    cache = {}
    cache_file = binary.resolve().parent.parent / "CMakeCache.txt"
    if cache_file.exists():
        for line in cache_file.read_text().splitlines():
            key, sep, value = line.partition("=")
            if sep and ":" in key:
                cache[key.split(":")[0]] = value
    cxx = cache.get("CMAKE_CXX_COMPILER", "c++")
    build_type = cache.get("CMAKE_BUILD_TYPE", "unknown")
    return {
        "git": command("git", "-C", str(REPO), "describe", "--always", "--dirty"),
        "machine": {"cpu": cpu, "os": os_name},
        "build": {
            "compiler": f"{command(cxx, '--version').splitlines()[0]} ({cxx})",
            "type": build_type,
            "flags": " ".join(
                f for f in (cache.get("CMAKE_CXX_FLAGS", ""), cache.get(f"CMAKE_CXX_FLAGS_{build_type.upper()}", "")) if f
            ),
        },
    }


def summarise(times):
    return {"median": statistics.median(times), "spread": max(times) - min(times)}


def save(result, label, folder):
    folder.mkdir(parents=True, exist_ok=True)
    stamp = datetime.fromisoformat(result["date"].replace("Z", "+00:00")).strftime("%Y%m%d-%H%M%S")
    path = folder / f"{stamp}_{label}.json"
    if path.exists():
        sys.exit(f"{path} exists already")
    path.write_text(json.dumps(result, indent=1) + "\n")
    return path


def warn_if_different(result, baseline):
    for key in ("cpu",):
        if result["machine"][key] != baseline["machine"][key]:
            print(f"warning: {key} differs from the baseline ({baseline['machine'][key]})")
    if result["build"] != baseline["build"]:
        print(f"warning: build differs from the baseline ({baseline['build']})")
    for name, run in result["runs"].items():
        if name in baseline["runs"] and run["parameters"] != baseline["runs"][name]["parameters"]:
            print(f"warning: run '{name}' has different parameters from the baseline; its times are not comparable")


def list_results(folder):
    files = sorted(folder.glob("*.json"))
    if not files:
        print(f"no results in {folder}")
        return
    for f in files:
        r = json.loads(f.read_text())
        medians = "  ".join(f"{n} {run['median']:.2f}" for n, run in r["runs"].items())
        print(f"{f.name}\n    {r['date']}  git {r['git']}  {r['machine']['cpu']}\n    {medians}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--binary", type=Path, default=REPO / "build/release/src/wormCP2Mu2dAFMplaq")
    p.add_argument("--config", type=Path, default=REPO / "tests/bench/bench.toml")
    p.add_argument("--only", action="append", help="time only this run (can be repeated)")
    p.add_argument("--repeats", type=int, default=3, help="times each run is timed (default 3)")
    p.add_argument("--save", metavar="LABEL", help="save the result as <UTC time>_<LABEL>.json")
    p.add_argument("--compare", type=Path, help="a saved result to compare with")
    p.add_argument("--results-dir", type=Path, default=RESULTS, help=f"where results are saved (default {RESULTS})")
    p.add_argument("--list", action="store_true", help="list the saved results and exit")
    args = p.parse_args()

    if args.list:
        list_results(args.results_dir)
        return
    if args.save and not re.fullmatch(r"[A-Za-z0-9._-]+", args.save):
        sys.exit("--save LABEL may contain only letters, digits, '.', '_' and '-'")
    names = run_names(args.config)
    if args.only:
        unknown = set(args.only) - set(names)
        if unknown:
            sys.exit(f"no run named {', '.join(sorted(unknown))} in {args.config}")
        names = [n for n in names if n in args.only]
    baseline = json.loads(args.compare.read_text()) if args.compare else None

    result = {
        "label": args.save,
        "date": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        **environment(args.binary),
        "repeats": args.repeats,
        "runs": {},
    }
    times = {n: [] for n in names}
    params = {}
    with tempfile.TemporaryDirectory(prefix="cp2bench_") as tmp:
        for r in range(args.repeats):
            for n in names:
                t, record = time_run(args.binary, args.config, n, Path(tmp))
                times[n].append(t)
                params.setdefault(n, parameters(record))
                print(f"round {r + 1}/{args.repeats}  {n:<12} {t:8.3f} s", flush=True)
    for n in names:
        result["runs"][n] = {"parameters": params[n], "seconds": times[n], **summarise(times[n])}

    header = f"\n{'run':<12} {'median s':>9} {'spread s':>9}"
    print(header + (f" {'baseline s':>11} {'change':>8}" if baseline else ""))
    for n, run in result["runs"].items():
        line = f"{n:<12} {run['median']:9.3f} {run['spread']:9.3f}"
        if baseline and n in baseline["runs"]:
            b = baseline["runs"][n]["median"]
            line += f" {b:11.3f} {100 * (run['median'] - b) / b:+7.2f}%"
        print(line)
    if baseline:
        warn_if_different(result, baseline)

    if args.save:
        print(f"saved to {save(result, args.save, args.results_dir)}")


if __name__ == "__main__":
    main()
