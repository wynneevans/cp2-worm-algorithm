import pytest

from cp2analysis.filters import apply_filters, parse_filter
from cp2analysis.runs import find_runs


def test_number_terms():
    p = parse_filter("0, 0.02", "number")
    assert p(0.0) and p(0.02) and not p(0.01)
    r = parse_filter("0.01:0.03", "number")
    assert r(0.01) and r(0.03) and not r(0.04)
    assert parse_filter("56.3", "number")(56.3) and parse_filter(":1", "number")(-5)


def test_lattice_terms():
    assert parse_filter("4x2", "lattice")("4x2x8")
    assert parse_filter("4x2x8", "lattice")("4x2x8") and not parse_filter("4x2x16", "lattice")("4x2x8")
    with pytest.raises(SystemExit):
        parse_filter("4by2", "lattice")


def test_apply_filters(data):
    d, _ = data
    runs = find_runs([str(d)])
    assert len(apply_filters(runs, {"beta": "1.5"})) == 3
    assert [r.mu3 for r in apply_filters(runs, {"beta": "1.5", "mu3": "0.01:"})] == [0.01, 0.02]
    assert len(apply_filters(runs, {"bc": "closed", "lattice": "2x2"})) == 2
    assert apply_filters(runs, {"mu8": "1"}) == []
    assert len(apply_filters(runs, {"mu8": "0.00"})) == 5
    assert len(apply_filters(runs, {"lx": "4", "lt": "8"})) == 3  # any Ly
    assert len(apply_filters(runs, {"ly": "1:2"})) == 5
