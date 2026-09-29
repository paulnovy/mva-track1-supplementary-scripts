"""Called ALT survives review selection without inventing the missing allele."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import pysam

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts/clinvar' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PartialReview(unittest.TestCase):
    def test_selection_counts_and_original_review_boundary(self):
        exact, triage = load('clinvar_exact'), load('clinvar_triage')
        header = pysam.VariantHeader()
        header.contigs.add('chr1', length=1000)
        header.add_sample('synthetic')
        header.formats.add('GT', 1, 'String', 'Genotype')
        header.formats.add('AD', 'R', 'Integer', 'Allelic depths')
        header.formats.add('PL', 'G', 'Integer', 'Likelihoods')
        header.formats.add('DP', 1, 'Integer', 'Depth')
        header.formats.add('GQ', 1, 'Integer', 'Genotype quality')
        header.info.add('ALLELEID', 1, 'Integer', 'Exact allele ID')
        header.info.add('CLNSIG', '.', 'String', 'Classification')
        for field in ('CLNREVSTAT', 'CLNDN', 'GENEINFO', 'MC'):
            header.info.add(field, '.', 'String', field)
        rows = [((0, 1), 'Pathogenic'), ((None, 1), 'Pathogenic'),
                ((1, None), 'Uncertain_significance'), ((0, None), 'Pathogenic'),
                ((None, None), 'Pathogenic'), ((0, 0), 'Pathogenic')]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'synthetic.vcf'
            with pysam.VariantFile(path, 'w', header=header) as out:
                for ident, (gt, classification) in enumerate(rows, 1):
                    record = out.new_record(contig='chr1', start=ident, alleles=('A', 'C'))
                    record.info['ALLELEID'] = ident
                    record.info['CLNSIG'] = classification
                    record.samples[0]['GT'] = gt
                    record.samples[0]['AD'] = (0, 12)
                    record.samples[0]['PL'] = (80, 60, 20)
                    out.write(record)
            counts = exact.aggregate(path)['counts']
            self.assertEqual(counts['matched_complete_nonreference_pathogenic_or_likely_pathogenic'], 1)
            self.assertEqual(counts['matched_called_alt'], 3)
            self.assertEqual(counts['matched_partial_called_alt'], 2)
            self.assertEqual(counts['matched_partial_called_alt_vus_or_uncertain'], 1)
            self.assertEqual(counts['matched_partial_called_alt_pathogenic_or_likely_pathogenic'], 1)
            self.assertEqual(counts['matched_partial_without_called_alt'], 2)
            _, selected = triage.collect(path)
            self.assertEqual(len(selected), 2)
            partial = next(x['call'] for key, x in selected.items() if key[-1] == 2)
            self.assertEqual(partial['gt'], [None, 1])
            self.assertTrue(partial['requires_original_genotype_review'])
            self.assertFalse(partial['checks']['complete_nonreference'])
            self.assertIsNone(partial['called_pl_is_min'])
            self.assertIsNone(partial['pl_margin'])


if __name__ == '__main__':
    unittest.main()
