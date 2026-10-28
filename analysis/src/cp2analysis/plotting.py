"""Plot helpers. The only module that touches matplotlib directly (seaborn draws through it)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402

BLUE = "#0072B2"  # Okabe-Ito, colour-blind safe
ORANGE = "#D55E00"


def apply_theme(context: str = "notebook") -> None:
    sns.set_theme(style="whitegrid", context=context)


def new_figure(nrows: int = 1, ncols: int = 1, **kwargs):
    return plt.subplots(nrows, ncols, **kwargs)


def errorbar(ax, x, y, yerr, **kwargs) -> None:
    """Points with error bars from precomputed errors (seaborn only draws errors it computes itself)."""
    ax.errorbar(x, y, yerr=yerr, color=kwargs.pop("color", BLUE), **kwargs)


def save(fig, path: str | Path, dpi: int = 130) -> None:
    sns.despine(fig)
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
