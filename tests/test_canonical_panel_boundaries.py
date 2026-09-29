#!/usr/bin/env python3
"""Focused allele-span test for internal canonical splice boundaries."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE = Path(__file__).parents[1] / "scripts" / "gene_audits" / "classify_canonical_panel_rows.py"
SPEC = importlib.util.spec_from_file_location("classifier", MODULE)
CLASSIFIER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CLASSIFIER)


class CanonicalBoundaryTest(unittest.TestCase):
    def test_minimal_allele_spans_and_strand_aware_internal_splice_roles(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            coordinates = directory / "coordinates.tsv"
            coordinates.write_text("PLUS\tENSGPLUS\t1\t100\t350\t1\tplus\nMINUS\tENSGMINUS\t1\t100\t350\t-1\tminus\n")
            for gene, strand in (("PLUS", 1), ("MINUS", -1)):
                (directory / f"ensembl_{gene}_canonical_structure.json").write_text(json.dumps({
                    "strand": strand,
                    "Translation": {"start": 100, "end": 350},
                    "Exon": [{"start": 100, "end": 150}, {"start": 200, "end": 250}],
                }))
            geometry = CLASSIFIER.load_geometry(coordinates, directory)
            plus_coding, plus_splice = geometry["PLUS"]
            minus_coding, minus_splice = geometry["MINUS"]

            # REF AT>ALT A changes base 151, not only its POS 150 anchor.
            deletion = CLASSIFIER.changed_span(150, "AT", "A")
            # REF A>ALT AG is inserted between 197/198; both flanks are retained.
            insertion = CLASSIFIER.changed_span(197, "A", "AG")
            self.assertEqual(deletion, (151, 151))
            self.assertEqual(insertion, (197, 198))
            self.assertEqual(plus_splice, [("donor", (151, 152)), ("acceptor", (198, 199))])
            self.assertEqual(minus_splice, [("acceptor", (151, 152)), ("donor", (198, 199))])
            self.assertEqual(CLASSIFIER.classify([deletion], plus_coding, plus_splice), "canonical_splice_intronic_donor")
            self.assertEqual(CLASSIFIER.classify([deletion], minus_coding, minus_splice), "canonical_splice_intronic_acceptor")
            self.assertEqual(CLASSIFIER.classify([insertion], plus_coding, plus_splice), "canonical_splice_intronic_acceptor")
            # A deletion reaching coding base 150 remains coding despite adjacency to donor flank.
            self.assertEqual(CLASSIFIER.classify([CLASSIFIER.changed_span(149, "AT", "A")], plus_coding, plus_splice), "coding")


if __name__ == "__main__":
    unittest.main()
