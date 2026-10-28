"""Finished runs in the data folder and reading of their M files.

Each finished run leaves M_<time>.dat, C_<time>.dat and run_<time>.toml (its record) in its parameter
folder, data/<bc>_<measurements>_<beta>_<Lx>_<Ly>_<Lt>_mu3_<mu3>_mu8_<mu8>/. All parameters are read
from the record.
"""

import tomllib
from dataclasses import dataclass
from pathlib import Path

import polars as pl

# M_*.dat columns: spatial/temporal windings, temporal charges (doubled units), worm lengths, worm tail.
COLUMNS = ["w3x", "w8x", "w3y", "w8y", "q3", "q8", "l", "lp", "lm", "tail"]
SERIES = COLUMNS[:9]  # numeric columns; `tail` is a string label


@dataclass(frozen=True)
class RunInfo:
    bc: str
    beta: float
    lx: int
    ly: int
    lt: int
    mu3: float
    mu8: float
    measurements: int = 0
    name: str = ""  # the run's name in its config
    stamp: str = ""  # UTC finish time from the file names
    seed: int | None = None
    path: Path | None = None  # the M_<time>.dat file

    @property
    def lattice(self) -> str:
        """Spatial lattice, as in the exact results: "LxxLy"."""
        return f"{self.lx}x{self.ly}"

    @property
    def size(self) -> str:
        return f"{self.lx}x{self.ly}x{self.lt}"

    @property
    def n_sites(self) -> int:
        return self.lx * self.ly

    @property
    def folder(self) -> str:
        """Name of the parameter folder holding this run."""
        return self.path.parent.name if self.path else ""

    @property
    def label(self) -> str:
        """Unique per run: <name>_<finish time>."""
        return f"{self.name}_{self.stamp}"


def record_for(m_file: str | Path) -> Path:
    """The run record belonging to M_<time>.dat: run_<time>.toml in the same folder."""
    p = Path(m_file)
    return p.with_name("run_" + p.stem.removeprefix("M_") + ".toml")


def read_run(m_file: str | Path) -> RunInfo:
    """The run that wrote this M_<time>.dat, from its record."""
    m = Path(m_file)
    record = record_for(m)
    if not m.name.startswith("M_") or not record.exists():
        raise SystemExit(f"error: '{m}' is not a finished run (no {record.name} next to it)")
    with open(record, "rb") as f:
        r = tomllib.load(f)["run"][0]
    lx, ly, lt = r["lattice"]
    return RunInfo(
        bc=r["bc"],
        beta=float(r["beta"]),
        lx=lx,
        ly=ly,
        lt=lt,
        mu3=float(r["mu3"]),
        mu8=float(r["mu8"]),
        measurements=int(r["measurements"]),
        name=r["name"],
        stamp=m.stem.removeprefix("M_"),
        seed=r.get("seed"),
        path=m,
    )


def find_runs(inputs: list[str]) -> list[RunInfo]:
    """Finished runs in the given M files and folders, sorted by parameter folder and finish time.

    A folder (the data folder or a parameter folder) is searched recursively for M_<time>.dat files
    with a run_<time>.toml; unfinished runs (M_running_*, no record) are skipped. An M file given
    explicitly must have its record.
    """
    paths: set[Path] = set()
    for i in inputs:
        p = Path(i)
        if p.is_dir():
            paths.update(f.resolve() for f in p.rglob("M_*.dat") if record_for(f).exists())
        elif p.exists():
            read_run(p)  # raises if it has no record
            paths.add(p.resolve())
        else:
            raise SystemExit(f"error: '{p}' does not exist")
    runs = [read_run(p) for p in paths]
    return sorted(runs, key=lambda r: (r.folder, r.stamp))


def read_m(path: str | Path) -> pl.DataFrame:
    """Read an M_*.dat file (no header, space-padded columns); numeric columns become Float64."""
    df = pl.read_csv(path, has_header=False, new_columns=COLUMNS, infer_schema=False)
    return df.with_columns([pl.col(c).str.strip_chars().cast(pl.Float64) for c in SERIES])
