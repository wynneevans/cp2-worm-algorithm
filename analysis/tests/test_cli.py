"""Each subcommand runs end to end on the synthetic data folder (see conftest.py)."""

import polars as pl
import pytest

from conftest import add_run
from cp2analysis.cli import main
from cp2analysis.commands.summary import choose_x, series


def test_summary_pools_each_parameter_folder(data):
    d, exact = data
    out = d / "out"
    main(["summary", str(d), "--exact", str(exact), "--out", str(out)])
    assert {p.name for p in out.iterdir()} == {"m_summary.csv", "susceptibility.png", "magnetisation.png"}
    s = pl.read_csv(out / "m_summary.csv")
    assert s.height == 4
    mini = s.filter(pl.col("bc") == "closed").row(0, named=True)
    assert mini["n_runs"] == 2 and mini["n_used"] == 2 * (6000 - 60)  # 1% skipped from each run
    assert mini["exact_chi"] == 0.008


def test_summary_per_run(data):
    d, exact = data
    out = d / "out"
    main(["summary", str(d), "--per-run", "--exact", str(exact), "--out", str(out)])
    s = pl.read_csv(out / "m_summary.csv")
    assert s.height == 5 and set(s["n_runs"]) == {1}
    assert "mini_20261028-130000" in s["label"].to_list()


def test_summary_of_a_scan_plots_against_mu3(data):
    d, exact = data
    main(["summary", "--beta", "1.5", str(d), "--exact", str(exact), "--out", str(d / "out")])
    s = pl.read_csv(d / "out" / "m_summary.csv")
    assert s["mu3"].to_list() == [0.0, 0.01, 0.02]
    assert choose_x(s, None) == "mu3"
    assert choose_x(pl.read_csv(d / "out" / "m_summary.csv").head(1), None) == "run"


def test_scan_over_mu3_and_ly_plots_mu3_with_a_series_per_ly(data):
    d, exact = data
    for i, mu3 in enumerate((0.0, 0.01)):
        add_run(
            d,
            f"20261030-10000{i}",
            f"ly3_{mu3}",
            bc="y-open",
            beta=1.5,
            lattice=(4, 3, 8),
            mu3=mu3,
            seed=9 + i,
        )
    main(
        [
            "summary",
            "--beta",
            "1.5",
            "--lx",
            "4",
            "--lt",
            "8",
            str(d),
            "--exact",
            str(exact),
            "--out",
            str(d / "o"),
        ]
    )
    s = pl.read_csv(d / "o" / "m_summary.csv")
    assert s.height == 5 and sorted(set(s["ly"])) == [2, 3]
    assert choose_x(s, None) == "mu3"
    assert [name for name, _ in series(s, "mu3")] == ["ly=2", "ly=3"]


def test_default_output_folder(data):
    d, exact = data
    folder = d / "closed_6000_0.05_2_2_4_mu3_0.000_mu8_0.000"
    main(["summary", str(folder), "--exact", str(exact)])
    assert (folder / "plots" / "m_summary.csv").exists()  # one parameter set: its own folder
    main(["summary", str(d), "--beta", "1.5", "--exact", str(exact)])
    assert (d / "plots" / "m_summary.csv").exists()  # several: the folder holding them


def test_raw(data):
    d, _ = data
    main(["raw", str(d), "--mu3", "0.02", "--out", str(d / "out"), "--block", "100"])
    assert [p.name for p in (d / "out").iterdir()] == ["raw_scan_0.02_20261029-100002.png"]


def test_acf(data):
    d, _ = data
    out = d / "out"
    main(["acf", str(d), "--bc", "closed", "--out", str(out), "--maxlag", "10"])
    assert (out / "acf_mini_20261028-120000.png").exists() and (out / "tau_int.csv").exists()


def test_thermalisation(data):
    d, exact = data
    out = d / "out"
    main(
        [
            "thermalisation",
            str(d),
            "--bc",
            "closed",
            "--exact",
            str(exact),
            "--out",
            str(out),
            "--nblocks",
            "20",
        ]
    )
    assert (out / "skipscan_mini_20261028-130000.png").exists() and (out / "mser.csv").exists()
    chunks = pl.read_csv(out / "chunks.csv")
    assert chunks.height == 2 * 8 and chunks["chi2dof_M3"].is_not_null().all()  # chi3: q3^2 is constant here


def test_nothing_selected_is_an_error(data, tmp_path):
    d, _ = data
    with pytest.raises(SystemExit, match="no finished runs"):
        main(["raw", str(d), "--mu3", "5"])
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(SystemExit, match="no finished runs"):
        main(["raw", str(empty)])


def test_check_exact_fails_on_a_wrong_value(data):
    d, exact = data
    mini = next(d.glob("closed_*"))
    # Synthetic q = +-1 gives chi3 = beta / (4 * 2 * 2) = 0.003125 exactly, not the CSV's 0.008.
    with pytest.raises(SystemExit, match="values disagree"):
        main(["check-exact", str(mini), "--exact", str(exact)])


def test_check_exact_fails_without_an_exact_value(data):
    d, exact = data
    with pytest.raises(SystemExit, match="no exact value"):
        main(["check-exact", str(d), "--beta", "1.5", "--mu3", "0.01", "--exact", str(exact)])


def test_check_exact_compares_the_charge_with_the_chemical_potential(data, tmp_path):
    d, _ = data
    csv = tmp_path / "typed.csv"
    csv.write_text(
        "lattice,boundary,beta_label,mu3,mu8,quantity,value,status,source,note\n"
        f"4x2,y-open,1.5,0.01,0.000,susceptibility,{1.5 / 32},typed in notebook; x,x,\n"
        "4x2,y-open,1.5,0.01,0.000,magnetization,0.0,typed in notebook; x,x,\n"
    )
    args = ["check-exact", str(d), "--beta", "1.5", "--mu3", "0.01", "--exact", str(csv)]
    with pytest.raises(SystemExit, match="no exact value"):  # typed rows are only used with --typed
        main(args)
    main([*args, "--typed"])
