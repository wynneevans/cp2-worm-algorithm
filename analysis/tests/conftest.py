"""A synthetic data folder: finished runs in parameter folders, as the simulator writes them."""

import numpy as np
import pytest

EXACT = (
    "lattice,boundary,beta_label,mu3,mu8,quantity,value,status,source,note\n"
    "2x2,closed,0.050,0.000,0.000,susceptibility,0.008,labelled file,x,\n"
    "2x2,closed,0.050,0.000,0.000,magnetization,0.,labelled file,x,\n"
)


def folder_name(bc, n, beta, lattice, mu3, mu8=0.0):
    lx, ly, lt = lattice
    return f"{bc}_{n}_{beta:.2f}_{lx}_{ly}_{lt}_mu3_{mu3:.3f}_mu8_{mu8:.3f}"


def record(name, bc, n, beta, lattice, mu3, mu8=0.0, seed=1):
    return (
        f'[[run]]\nname = "{name}"\ndata_dir = ".."\nbc = "{bc}"\nmeasurements = {n}\nbeta = {beta}\n'
        f"lattice = {list(lattice)}\nmu3 = {mu3}\nmu8 = {mu8}\nseed = {seed}\n"
    )


def write_m(path, n, seed=0):
    rng = np.random.default_rng(seed)
    q = rng.choice([-1, 1], size=(n, 2))
    with open(path, "w") as f:
        for q3, q8 in q:
            vals = [0, 0, 0, 0, q3, q8, 4, 4, 1]
            f.write(",".join(f"{v:<15f}" for v in vals) + ",vm             \n")


def add_run(data, stamp, name, bc="closed", n=6000, beta=0.05, lattice=(2, 2, 4), mu3=0.0, seed=0):
    folder = data / folder_name(bc, n, beta, lattice, mu3)
    folder.mkdir(parents=True, exist_ok=True)
    write_m(folder / f"M_{stamp}.dat", n, seed)
    (folder / f"run_{stamp}.toml").write_text(record(name, bc, n, beta, lattice, mu3, seed=seed))
    return folder


@pytest.fixture
def data(tmp_path):
    """data/ with an exact-size case (two runs) and a small mu3 scan at beta=1.5; plus an exact CSV."""
    d = tmp_path / "data"
    add_run(d, "20261028-120000", "mini", seed=1)
    add_run(d, "20261028-130000", "mini", seed=2)  # a repeat: same parameter folder
    for i, mu3 in enumerate((0.0, 0.01, 0.02)):
        add_run(
            d,
            f"20261029-10000{i}",
            f"scan_{mu3}",
            bc="y-open",
            beta=1.5,
            lattice=(4, 2, 8),
            mu3=mu3,
            seed=3 + i,
        )
    exact = tmp_path / "exact.csv"
    exact.write_text(EXACT)
    return d, exact
