import pytest

from cp2analysis.runs import RunInfo, find_m_files, parse_name, run_info

NAME = "M_wormCP2AFMmu2dplaq_BC_y-open_nMeas_1000000_beta_0.05_XYT_2_3_4_mu3_0.000_mu8_0.000_nSim_1.dat"


def test_parse_name():
    run = parse_name(NAME)
    assert (run.bc, run.beta, run.lx, run.ly, run.lt, run.mu3, run.mu8) == ("y-open", 0.05, 2, 3, 4, 0.0, 0.0)
    assert run.name == "2x3_y-open_b0.05"
    assert run.lattice == "2x3" and run.n_sites == 6


def test_parse_name_rejects_other_files():
    with pytest.raises(ValueError):
        parse_name("notes.txt")


def test_find_m_files(tmp_path):
    (tmp_path / NAME).write_text("")
    (tmp_path / "C_other.dat").write_text("")
    assert find_m_files([str(tmp_path)]) == []  # no run record next to it
    assert find_m_files([str(tmp_path / NAME)]) == [tmp_path / NAME]


RECORD = """# Record of what was run.
[[run]]
name = "2x3_yopen"
data_dir = ".."
bc = "y-open"
measurements = 1000000
beta = 0.05
lattice = [2, 3, 4]
mu3 = 0.0
mu8 = 0.0
seed = 20250103

[run.probabilities]
dt = 0.05
mu3 = 0.0
mu8 = 0.0
p01 = 1.0

[provenance]
started = '2026-10-28T10:41:15Z'
git = 'abc1234'
"""


def test_run_folder(tmp_path):
    folder = tmp_path / "y-open_1000000_0.05_2_3_4_mu3_0.000_mu8_0.000"
    folder.mkdir()
    for stamp in ("20261028-114631", "20261028-114631-2"):
        (folder / f"M_{stamp}.dat").write_text("")
        (folder / f"run_{stamp}.toml").write_text(RECORD)
    (folder / "M_running_123.dat").write_text("")  # an unfinished run has no record and is skipped
    (tmp_path / NAME).write_text("")  # old flat files in a searched folder are skipped too

    found = find_m_files([str(tmp_path)])
    assert found == [folder / "M_20261028-114631-2.dat", folder / "M_20261028-114631.dat"]
    run = run_info(folder / "M_20261028-114631.dat")
    assert run == RunInfo(
        bc="y-open",
        beta=0.05,
        lx=2,
        ly=3,
        lt=4,
        mu3=0.0,
        mu8=0.0,
        name="2x3_yopen_20261028-114631",
        seed=20250103,
    )
    assert run_info(tmp_path / NAME).name == "2x3_y-open_b0.05"  # given explicitly, an old file still works
