import pytest

from conftest import add_run
from cp2analysis.runs import find_runs, read_run


def test_find_runs_reads_records(data):
    d, _ = data
    runs = find_runs([str(d)])
    assert len(runs) == 5
    mini = [r for r in runs if r.name == "mini"]
    assert [r.stamp for r in mini] == ["20261028-120000", "20261028-130000"]
    r = mini[0]
    assert (r.bc, r.beta, r.lx, r.ly, r.lt, r.mu3, r.mu8, r.measurements) == (
        "closed",
        0.05,
        2,
        2,
        4,
        0.0,
        0.0,
        6000,
    )
    assert r.folder == "closed_6000_0.05_2_2_4_mu3_0.000_mu8_0.000"
    assert r.label == "mini_20261028-120000" and r.lattice == "2x2" and r.size == "2x2x4"


def test_parameter_folder_and_file_inputs_without_duplicates(data):
    d, _ = data
    folder = d / "closed_6000_0.05_2_2_4_mu3_0.000_mu8_0.000"
    assert len(find_runs([str(folder)])) == 2
    assert len(find_runs([str(folder), str(folder / "M_20261028-120000.dat"), str(d)])) == 5


def test_unfinished_and_stray_files_are_skipped(tmp_path):
    folder = add_run(tmp_path, "20261028-120000", "a")
    (folder / "M_running_123.dat").write_text("")
    (tmp_path / "M_wormCP2AFMmu2dplaq_BC_closed_nMeas_10_nSim_1.dat").write_text("")  # old flat file
    assert [r.label for r in find_runs([str(tmp_path)])] == ["a_20261028-120000"]


def test_file_without_record_is_an_error(tmp_path):
    flat = tmp_path / "M_wormCP2AFMmu2dplaq_BC_closed_nMeas_10_nSim_1.dat"
    flat.write_text("")
    with pytest.raises(SystemExit, match="not a finished run"):
        read_run(flat)
    with pytest.raises(SystemExit, match="not a finished run"):
        find_runs([str(flat)])
    with pytest.raises(SystemExit, match="does not exist"):
        find_runs([str(tmp_path / "nope")])
