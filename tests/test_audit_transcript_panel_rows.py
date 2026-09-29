#!/usr/bin/env python3
"""Focused geometry check for primary versus secondary transcript reporting."""
import importlib.util
import sys
from pathlib import Path

MODULE = Path(__file__).parents[1] / "scripts" / "gene_audits" / "audit_transcript_panel_rows.py"
sys.path.insert(0, str(MODULE.parent))
SPEC = importlib.util.spec_from_file_location("audit", MODULE)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_secondary_coding_and_splice_flank_are_distinguished():
    primary = AUDIT.transcript_geometry({"strand": 1, "Translation": {"start": 100, "end": 150}, "Exon": [{"start": 100, "end": 150}]}, 8)
    secondary = AUDIT.transcript_geometry({"strand": 1, "Translation": {"start": 100, "end": 250}, "Exon": [{"start": 100, "end": 150}, {"start": 200, "end": 250}]}, 8)
    assert AUDIT.region([(196, 196)], primary) == "other"
    assert AUDIT.region([(196, 196)], secondary) == "splice_acceptor"
    assert AUDIT.region([(220, 220)], secondary) == "coding"


def test_gt_value_preserves_phasing_separator():
    class Sample(dict):
        phased = True
    class Record:
        samples = {"sample": Sample(GT=(0, 1))}
    assert AUDIT.gt_value(Record(), "sample") == "0|1"


def test_filter_value_does_not_turn_dot_into_pass():
    class Filter:
        def keys(self):
            return []
    class Record:
        filter = Filter()
    assert AUDIT.filter_value(Record()) == "."


def load_tests(loader, tests, pattern):
    import unittest
    tests.addTests(unittest.FunctionTestCase(function) for name, function in globals().items()
                   if name.startswith('test_') and callable(function))
    return tests
