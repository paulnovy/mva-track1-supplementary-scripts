"""Synthetic checks that candidate reports use inputs, not historical results."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam


SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts' / 'candidates'


def load_annotation():
    spec = importlib.util.spec_from_file_location('candidate_annotation', SCRIPTS/'focused_candidate_annotation.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CandidateReportTest(unittest.TestCase):
    def test_targeted_counts_follow_synthetic_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sam = root/'reads.sam'
            sam.write_text('')
            depth = root/'depth.tsv'
            depth.write_text('')
            dfam = root/'dfam.json'
            dfam.write_text(json.dumps({'consensus_sequence':'ACGT' * 20}))
            auxiliary = root/'empty.txt'
            auxiliary.write_text('')
            out = root/'output'
            command = [sys.executable, str(SCRIPTS/'analyze_targeted_sv_mei.py'),
                '--bub-sam', str(sam), '--bub-depth', str(depth), '--clcn-sam', str(sam),
                '--clcn-depth', str(depth), '--clcn-depth-mapq20', str(depth),
                '--clcn-depth-mapq60', str(depth), '--dfam', str(dfam),
                '--delly-records', str(auxiliary), '--clcn-realignment', str(auxiliary),
                '--output', str(out)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            empty = json.loads((out/'targeted-evidence.json').read_text())
            self.assertEqual(empty['bub1b_published_alu']['clip_support']['unique_fragments'], 0)
            self.assertIsNone(empty['clcnkb_deletion']['target_to_mean_control_depth_ratio_mapq20'])
            sam.write_text('\t'.join(['synthetic-read', '0', 'chr15', '40196543', '60',
                '10S20M', '*', '0', '0', 'A'*30, 'I'*30]) + '\n')
            subprocess.run(command, check=True, capture_output=True, text=True)
            populated = json.loads((out/'targeted-evidence.json').read_text())
            self.assertEqual(populated['bub1b_published_alu']['clip_support']['unique_fragments'], 1)
            self.assertEqual(populated['bub1b_published_alu']['status'], empty['bub1b_published_alu']['status'])
            self.assertIn('interpretation_required', populated['bub1b_published_alu']['status'])

    def test_annotation_select_assemble_verify_uses_manifest_and_actual_counts(self):
        module = load_annotation()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root/'reference.fa'
            fasta.write_text('>chr1\n' + 'A'*200 + '\n')
            pysam.faidx(str(fasta))
            header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            for name, record in [('source', 'chr1\t101\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\n'), ('clinvar', '')]:
                path = root/(name+'.vcf')
                path.write_text(header + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsynthetic\n' + record)
                pysam.tabix_compress(str(path), str(path)+'.gz', force=True)
                pysam.tabix_index(str(path)+'.gz', preset='vcf', force=True)
            exroot = root/'exomiser'
            for profile in ('full', 'axis'):
                directory = exroot/profile/'work/rare-disease/native'
                directory.mkdir(parents=True)
                (directory/f'MVA1-{profile}-exomiser.variants.tsv').write_text(
                    'GENE_SYMBOL\tCONTIG\tSTART\tREF\tALT\tALL_FREQ\tHGVS\n'
                    'SYNTHETIC\t1\t101\tA\tT\t\t\n')
            manifest = root/'manifest.json'
            manifest.write_text(json.dumps({'genes':['SYNTHETIC'], 'audit_genes':['SYNTHETIC'],
                'input_paths':{'reference_fasta':fasta.name, 'normalized_vcf':'source.vcf.gz',
                    'clinvar_vcf':'clinvar.vcf.gz', 'exomiser_root':'exomiser'},
                'inheritance_annotations':[{'gene':'SYNTHETIC', 'phase':'cis_supported'}]}))
            module.configure_manifest(manifest)
            module.RES, module.REP, module.OPS = (root/name for name in ('results', 'reports', 'operations'))
            for path in (module.RES, module.REP, module.OPS):
                path.mkdir()
            module.select()
            selected = json.loads((module.RES/'selection.json').read_text())
            self.assertEqual(len(selected), 1)
            self.assertIsNone(selected[0]['clinvar_snapshot_date'])
            self.assertIn('present_in_ranked_small_variant_outputs', (module.RES/'shortlist_exclusions.tsv').read_text())
            archive = module.RES/'vep_release115'
            archive.mkdir()
            candidate_id = selected[0]['id']
            (archive/(candidate_id+'.json')).write_text(json.dumps([{
                'assembly_name':'GRCh38', 'seq_region_name':'1', 'start':101, 'allele_string':'A/T',
                'transcript_consequences':[{'gene_symbol':'SYNTHETIC', 'canonical':1, 'transcript_id':'SYNTHETIC_TX'}]}]))
            (archive/(candidate_id+'.stdout.txt')).write_text('')
            (archive/'ensembl_data_release.json').write_text(json.dumps({'releases':[115]}))
            gnomad = module.RES/'gnomad_qc'
            gnomad.mkdir()
            (gnomad/(candidate_id+'.json')).write_text(json.dumps({'data':{'variant':{
                'variant_id':candidate_id, 'flags':[], 'exome':None, 'genome':None, 'joint':None}}}))
            module.assemble()
            with contextlib.redirect_stdout(io.StringIO()):
                module.verify()
            verification = json.loads((module.OPS/'verification.json').read_text())
            self.assertEqual(verification['checks']['selected_variants'], 1)
            self.assertEqual(verification['checks']['clinvar_exact_snapshot_matches'], 0)
            self.assertNotIn('focused_unit_check', verification['checks'])
            inheritance = json.loads((module.RES/'inheritance_audit.json').read_text())
            self.assertEqual(len(inheritance), 1)
            self.assertIn('not inferred', inheritance[0]['annotation_source'])
            self.assertTrue(inheritance[0]['phase_rule_interpretation'].startswith('nominated_pair_not_biallelic'))
            completion = json.loads((module.OPS/'completion.json').read_text())
            self.assertTrue(Path(completion['report']).is_file())


if __name__ == '__main__':
    unittest.main()
