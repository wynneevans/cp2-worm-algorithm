"""Run parameters and reading of M files.

Each finished run leaves M_<time>.dat, C_<time>.dat and run_<time>.toml (its parameters) in
data/<bc>_<measurements>_<beta>_<Lx>_<Ly>_<Lt>_mu3_<mu3>_mu8_<mu8>/. Old flat files
(M_wormCP2AFMmu2dplaq_BC_..._nSim_1.dat) are still understood through their long names.
"""

import re
import tomllib
from dataclasses import dataclass
from pathlib import Path

import polars as pl

# M_*.dat columns: spatial/temporal windings, temporal charges (doubled units), worm lengths, worm tail.
COLUMNS = ["w3x", "w8x", "w3y", "w8y", "q3", "q8", "l", "lp", "lm", "tail"]
SERIES = COLUMNS[:9]  # numeric columns; `tail` is a string label

_NAME_RE = re.compile(
    r"BC_(?P<bc>[^_]+)_nMeas_(?P<nmeas>\d+)_beta_(?P<beta>[\d.]+)_XYT_(?P<lx>\d+)_(?P<ly>\d+)_(?P<lt>\d+)"
    r"_mu3_(?P<mu3>-?[\d.]+)_mu8_(?P<mu8>-?[\d.]+)_nSim_(?P<sim>\d+)"
)


@dataclass(frozen=True)
class RunInfo:
    bc: str
    beta: float
    lx: int
    ly: int
    lt: int
    mu3: float
    mu8: float
    name: str = ""
    sim: int = 1
    seed: int | None = None

    @property
    def lattice(self) -> str:
        return f"{self.lx}x{self.ly}"

    @property
    def n_sites(self) -> int:
        return self.lx * self.ly


def _legacy_name(run: RunInfo) -> str:
    name = f"{run.lattice}_{run.bc}_b{run.beta:g}"
    return name + (f"_mu3_{run.mu3:g}" if run.mu3 else "")


def parse_name(path: str | Path) -> RunInfo:
    """Run parameters encoded in an old flat M_/C_ file name."""
    m = _NAME_RE.search(Path(path).name)
    if not m:
        raise ValueError(f"cannot parse run parameters from '{path}'")
    d = m.groupdict()
    run = RunInfo(
        bc=d["bc"],
        beta=float(d["beta"]),
        lx=int(d["lx"]),
        ly=int(d["ly"]),
        lt=int(d["lt"]),
        mu3=float(d["mu3"]),
        mu8=float(d["mu8"]),
        sim=int(d["sim"]),
    )
    return RunInfo(**{**run.__dict__, "name": _legacy_name(run)})


def read_record(path: str | Path) -> RunInfo:
    """Run parameters from a run record (run_<time>.toml)."""
    with open(path, "rb") as f:
        r = tomllib.load(f)["run"][0]
    lx, ly, lt = r["lattice"]
    stamp = Path(path).stem.removeprefix("run_")
    return RunInfo(
        bc=r["bc"],
        beta=float(r["beta"]),
        lx=lx,
        ly=ly,
        lt=lt,
        mu3=float(r["mu3"]),
        mu8=float(r["mu8"]),
        name=f"{r['name']}_{stamp}",
        seed=r.get("seed"),
    )


def record_for(m_file: str | Path) -> Path:
    """The run record belonging to M_<time>.dat: run_<time>.toml in the same folder."""
    p = Path(m_file)
    return p.with_name("run_" + p.stem.removeprefix("M_") + ".toml")


def run_info(m_file: str | Path) -> RunInfo:
    """Parameters of the run that wrote this M file (from its record, or an old flat file's name)."""
    record = record_for(m_file)
    return read_record(record) if record.exists() else parse_name(m_file)


def find_m_files(inputs: list[str]) -> list[Path]:
    """Expand files and directories into a sorted list of M files.

    A directory is searched recursively for finished runs (M_<time>.dat with its run_<time>.toml).
    Old flat M_*.dat files are only used when given as files.
    """
    files: list[Path] = []
    for i in inputs:
        p = Path(i)
        if p.is_dir():
            files += sorted(f for f in p.rglob("M_*.dat") if record_for(f).exists())
        else:
            files.append(p)
    return files


def read_m(path: str | Path) -> pl.DataFrame:
    """Read an M_*.dat file (no header, space-padded columns); numeric columns become Float64."""
    df = pl.read_csv(path, has_header=False, new_columns=COLUMNS, infer_schema=False)
    return df.with_columns([pl.col(c).str.strip_chars().cast(pl.Float64) for c in SERIES])
