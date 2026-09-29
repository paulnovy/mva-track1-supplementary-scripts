"""Exercise spanning SVs and fresh-output safety without genomic inputs."""
import csv
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('structure', ROOT / 'scripts/gene_audits/analyze_trim37_cep192_structure.py')
STRUCTURE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(STRUCTURE)


class GeneAuditReplayTest(unittest.TestCase):
    def test_whole_span_called_alt_and_normalized_bins(self):
        # Both deletion breakpoints are outside CEP192; an INV encloses TRIM37.
        # Reference-genotype DEL is retained as evidence, but not dosage-reviewed.
        fields = 'id chrom pos end svtype filter qual pe sr mapq srmapq gt gq ft dr dv rr rv alt consensus'.split()
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sv = tmp / 'sv.tsv'
            with sv.open('w') as handle:
                writer = csv.writer(handle, delimiter='\t')
                for ident, chrom, start, end, kind, gt, flt in [
                    ('D', 'chr18', 12000000, 14000000, 'DEL', '0/1', '.'),
                    ('R', 'chr18', 12000000, 14000000, 'DEL', '0/0', 'LowQual'),
                    ('I', 'chr17', 58000000, 60000000, 'INV', '0|1', 'PASS'),
                ]:
                    row = dict.fromkeys(fields, '.')
                    row.update(id=ident, chrom=chrom, pos=start, end=end, svtype=kind, gt=gt, filter=flt, alt=f'<{kind}>')
                    writer.writerow([row[key] for key in fields])
            bins = tmp / 'bins.tsv'
            bins.write_text('chrom\tstart\tend\tclean\tcorrected_rd\tmedian_abs_baf_minus_half\n'
                            'chr18\t12900000\t13000000\t1\t1.02\t0.06\n'
                            'chr18\t13100000\t13200000\t0\t.\t.\n'
                            'chr18\t14100000\t14200000\t1\t1.01\t0.05\n'
                            'chr17\t58900000\t59000000\t1\t1.00\t.\n')
            hits = STRUCTURE.load_sv(sv)
            self.assertEqual(len(hits), 3)
            reviewed = {row['id']: row for row in STRUCTURE.span_bin_audit(bins, hits)}
            self.assertEqual(set(reviewed), {'D', 'I'})
            self.assertEqual(reviewed['D']['filter'], '.')
            self.assertEqual(reviewed['D']['inside']['mean_corrected_rd'], 1.02)
            self.assertEqual(reviewed['D']['masked_or_uncallable_overlapping_bins'], 1)
            self.assertEqual(reviewed['I']['breakpoints_inside_gene']['TRIM37'], [False, False])
            self.assertIn('unvalidated', reviewed['I']['interpretation'])
            self.assertIsNone(reviewed['I']['inside']['median_abs_baf_minus_half'])

    def test_runner_help_and_existing_output_are_nonmutating(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = dict(os.environ, AUDIT_OUTPUT=tmp, AUDIT_SHARED='missing', CORE_BAM='missing',
                       DELLY_BCF='missing', MASKED_BINS='missing', REFERENCE_FASTA='missing',
                       VENDOR_VCF='missing', NORMALIZED_VCF='missing', CORE_VCF='missing')
            marker = Path(tmp) / 'preserve.txt'
            marker.write_text('untouched')
            for name in ('variants', 'structure'):
                script = ROOT / f'scripts/gene_audits/run_trim37_cep192_{name}.sh'
                self.assertEqual(subprocess.run(['bash', str(script), '--help'], env=env, capture_output=True).returncode, 0)
                run = subprocess.run(['bash', str(script)], env=env, capture_output=True, text=True)
                self.assertEqual(run.returncode, 2)
                self.assertIn('existing output', run.stderr)
            self.assertEqual(marker.read_text(), 'untouched')


if __name__ == '__main__':
    unittest.main()
