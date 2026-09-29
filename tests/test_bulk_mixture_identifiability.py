import importlib.util
import sys
from pathlib import Path


MODULE = Path(__file__).parents[1] / "scripts" / "exploratory" / "bulk_mixture_identifiability.py"
SPEC = importlib.util.spec_from_file_location("bulk_mixtures", MODULE)
bulk_mixtures = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = bulk_mixtures
SPEC.loader.exec_module(bulk_mixtures)


def test_equal_gain_loss_haplotype_mixture_has_diploid_bulk_marginals():
    a, b, rd, baf = bulk_mixtures.totals(
        {"AB": 0.60, "AAB": 0.10, "ABB": 0.10, "A": 0.10, "B": 0.10}
    )
    assert (a, b, rd, baf) == (1.0, 1.0, 1.0, 0.5)


def test_fixed_cell_enumeration_matches_balanced_near_cancellation_formula():
    a, b, rd, baf = bulk_mixtures.enumerate_cells(
        {"AB": 160, "AAB": 11, "ABB": 11, "A": 9, "B": 9}
    )
    assert (a, b, rd, baf) == (1.01, 1.01, 1.01, 0.5)


def load_tests(loader, tests, pattern):
    import unittest
    tests.addTests(unittest.FunctionTestCase(function) for name, function in globals().items()
                   if name.startswith('test_') and callable(function))
    return tests
