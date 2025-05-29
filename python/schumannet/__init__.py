"""
scHumanNet: Construction of cell-type-specific functional gene network

This package provides tools for constructing cell-type-specific functional gene networks
using SCINET and HumanNetv3, with functions for downstream analysis including hub gene
identification and differential network analysis.
"""

__version__ = "0.1.0"
__author__ = "Junha Cha"
__email__ = "junhacha@yonsei.ac.kr"

from .core import (
    sort_add_lls,
    get_centrality,
    combine_perc_rank,
    top_hub,
    diff_pr,
    find_diff_hub,
    find_all_hub,
)

from .utils import (
    load_humannet_v3,
    connectivity,
    deconvolute_net,
)

__all__ = [
    "sort_add_lls",
    "get_centrality", 
    "combine_perc_rank",
    "top_hub",
    "diff_pr",
    "find_diff_hub",
    "find_all_hub",
    "load_humannet_v3",
    "connectivity",
    "deconvolute_net",
]