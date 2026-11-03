"""Small plotting helpers shared across metric modules."""
import numpy as np


def sig_stars(p) -> str:
    if np.isnan(p):
        return ""
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
