"""A listed but uncalled ALT must not become a clinical match."""
import importlib.util
from pathlib import Path
import unittest

import pysam

path = Path(__file__).resolve().parents[1] / 'scripts/gene_audits/audit_panel_clinvar_snapshot.py'
spec = importlib.util.spec_from_file_location('snapshot',path)
snapshot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(snapshot)


class CalledAlleleTest(unittest.TestCase):
    def test_only_genotyped_alt(self):
        header = pysam.VariantHeader()
        header.contigs.add('chr18',length=80373285)
        header.formats.add('GT',1,'String','Genotype')
        header.add_sample('synthetic')
        record = header.new_record(contig='chr18',start=100,alleles=('A','G','T'))
        for genotype, expected in [((0,1),[1]),((2,2),[2]),((0,0),[]),((None,None),[])]:
            record.samples['synthetic']['GT'] = genotype
            self.assertEqual(snapshot.called_alts(record),expected)


if __name__ == '__main__':
    unittest.main()
