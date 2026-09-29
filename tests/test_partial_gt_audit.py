"""Focused representation/inheritance regression: never invent the second allele."""
import importlib.util
import io
from pathlib import Path
import unittest

path = Path(__file__).resolve().parents[1] / 'scripts/qc/audit_partial_gt.py'
spec = importlib.util.spec_from_file_location('partial_gt_audit', path)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class PartialGTTest(unittest.TestCase):
    def test_complete_versus_presence_and_original_likelihood_semantics(self):
        import pysam
        header = pysam.VariantHeader()
        header.contigs.add('chr1', length=1000)
        for name, number, kind in [('GT', 1, 'String'), ('AD', 'R', 'Integer'), ('PL', 'G', 'Integer')]:
            header.formats.add(name, number, kind, 'synthetic')
        header.add_sample('synthetic')
        original = header.new_record(contig='chr1', start=100, alleles=('A', 'C', 'G'))
        original.samples[0]['GT'] = (1, 2)
        original.samples[0]['AD'] = (3, 12, 15)
        original.samples[0]['PL'] = (90, 80, 70, 60, 0, 50)
        self.assertEqual(audit.called_alts(original), [1, 2])
        for index, gt in [(1, (1, None)), (2, (None, 1))]:
            split = header.new_record(contig='chr1', start=100, alleles=('A', original.alts[index-1]))
            split.samples[0]['GT'] = gt
            self.assertEqual(audit.state(split), 'missing_or_partial')
            self.assertEqual(audit.called_alts(split), [1])
            self.assertEqual(audit.gt_text(original.samples[0]['GT']), '1/2')
        self.assertEqual(original.samples[0]['AD'], (3, 12, 15))
        self.assertEqual(original.samples[0]['PL'][4], 0)
        self.assertEqual(list(audit.proto_fields(io.BytesIO(b'\x08\x96\x01\x12\x02ok'))), [(1, 150), (2, b'ok')])

    def test_saved_transcript_strand_coordinates_and_insertion_flanks(self):
        # chr id 1, begin 100, end 110; no embedded Contig (the retained format).
        forward = b'\x18\x01\x20\x64\x28\x6e'
        reverse = b'\x10\x01' + forward
        self.assertEqual(audit.interval(forward, {1: '1'}, {1: 1000}), ('chr1', 101, 110))
        self.assertEqual(audit.interval(reverse, {1: '1'}, {1: 1000}), ('chr1', 891, 900))
        bins = {('chr1', 0): {(101, 110, 'synthetic', 'tx')}}
        self.assertEqual(audit.coding_hits(('chr1', 100, 'A', 'AT'), bins), [('synthetic', 'tx')])
        self.assertEqual(audit.coding_hits(('chr1', 100, 'A', 'T'), bins), [])


if __name__ == '__main__':
    unittest.main()
