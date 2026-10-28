"""Each subcommand runs end to end on a small synthetic M_ file."""

import numpy as np
import pytest

from cp2analysis.cli import main

NAME = "M_wormCP2AFMmu2dplaq_BC_closed_nMeas_20000_beta_0.05_XYT_2_2_4_mu3_0.000_mu8_0.000_nSim_1.dat"


RECORD = """[[run]]
name = "mini"
data_dir = ".."
bc = "closed"
measurements = 20000
beta = 0.05
lattice = [2, 2, 4]
mu3 = 0.0
mu8 = 0.0
seed = 1
"""

EXACT = (
    "lattice,boundary,beta_label,mu3,mu8,quantity,value,status,source,note\n"
    "2x2,closed,0.050,0.000,0.000,susceptibility,0.008,labelled file,x,\n"
    "2x2,closed,0.050,0.000,0.000,magnetization,0.,labelled file,x,\n"
)


def write_m(path):
    rng = np.random.default_rng(0)
    q = rng.choice([-1, 1], size=(20_000, 2))
    with open(path, "w") as f:
        for q3, q8 in q:
            vals = [0, 0, 0, 0, q3, q8, 4, 4, 1]
            f.write(",".join(f"{v:<15f}" for v in vals) + ",vm             \n")


@pytest.fixture
def data(tmp_path):
    """A data folder holding one finished run (M_<time>.dat + run_<time>.toml), and an exact-results CSV."""
    folder = tmp_path / "data" / "closed_20000_0.05_2_2_4_mu3_0.000_mu8_0.000"
    folder.mkdir(parents=True)
    write_m(folder / "M_20261028-120000.dat")
    (folder / "run_20261028-120000.toml").write_text(RECORD)
    exact = tmp_path / "exact.csv"
    exact.write_text(EXACT)
    return tmp_path / "data", exact


def test_summary(data):
    d, exact = data
    out = d / "out"
    main(["summary", str(d), "--exact", str(exact), "--out", str(out)])
    assert {p.name for p in out.iterdir()} == {"m_summary.csv", "susceptibility.png", "magnetisation.png"}
    assert "mini_20261028-120000" in (out / "m_summary.csv").read_text()


def test_default_output_is_plots_next_to_the_input(data):
    d, exact = data
    main(["summary", str(d), "--exact", str(exact)])
    assert (d / "plots" / "m_summary.csv").exists()


def test_raw(data):
    d, _ = data
    out = d / "out"
    main(["raw", str(d), "--out", str(out), "--block", "100"])
    assert (out / "raw_mini_20261028-120000.png").exists()


def test_acf(data):
    d, _ = data
    out = d / "out"
    main(["acf", str(d), "--out", str(out), "--maxlag", "10"])
    assert (out / "acf_mini_20261028-120000.png").exists() and (out / "tau_int.csv").exists()


def test_thermalisation(data):
    d, exact = data
    out = d / "out"
    main(["thermalisation", str(d), "--exact", str(exact), "--out", str(out)])
    assert (out / "skipscan_mini_20261028-120000.png").exists() and (out / "mser.csv").exists()


def test_old_flat_file_still_works(tmp_path):
    write_m(tmp_path / NAME)
    exact = tmp_path / "exact.csv"
    exact.write_text(EXACT)
    out = tmp_path / "out"
    main(["summary", str(tmp_path / NAME), "--exact", str(exact), "--out", str(out)])
    assert "2x2_closed_b0.05" in (out / "m_summary.csv").read_text()


def test_no_files_is_an_error(tmp_path):
    with pytest.raises(SystemExit):
        main(["raw", str(tmp_path), "--out", str(tmp_path / "o")])
