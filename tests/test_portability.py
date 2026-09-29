import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]

class PortabilityTest(unittest.TestCase):
    def test_offline_cooccurrence_extraction(self):
        response = {"data": {"variant_cooccurrence": {
            "variant_ids": ["v1", "v2"], "genotype_counts": {"double_het": 3},
            "haplotype_counts": {"trans": 2, "cis": 1},
            "p_compound_heterozygous": 0.5, "populations": []}}}
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / "response.json", Path(tmp) / "result.json"
            source.write_text(json.dumps(response))
            subprocess.run([sys.executable, str(ROOT / "scripts/phase/query_gnomad_schema.py"),
                "--variants", "v1", "v2", "--saved-response", str(source),
                "--output", str(output)], check=True)
            result = json.loads(output.read_text())
            self.assertEqual(result["summary"]["haplotype_counts"]["trans"], 2)
            self.assertEqual(result["summary"]["requested_variant_ids"], ["v1", "v2"])

    def test_network_modes_are_opt_in_and_help_is_safe(self):
        env = os.environ.copy()
        env.pop("MVA_ALLOW_NETWORK", None)
        gated = subprocess.run([sys.executable, str(ROOT / "scripts/phase/query_gnomad_schema.py"),
            "--variants", "v1", "v2", "--wrapper", str(ROOT / "scripts/public_api/get_variant_frequency.py"),
            "--output", "/dev/null"], env=env, capture_output=True, text=True)
        self.assertEqual(gated.returncode, 2)
        self.assertIn("MVA_ALLOW_NETWORK=1", gated.stderr)
        for script in ("scripts/candidates/focused_candidate_annotation.py", "scripts/clcnkb/audit_clcnkb_ins79.py"):
            help_run = subprocess.run([sys.executable, str(ROOT / script), "--help"], env=env, capture_output=True)
            self.assertEqual(help_run.returncode, 0)
        direct = subprocess.run([sys.executable, str(ROOT / "scripts/public_api/get_variant_frequency.py"),
            "--variant_id", "1-1-A-C", "--output", "/dev/null"], env=env, capture_output=True, text=True)
        self.assertEqual(direct.returncode, 2)
        self.assertIn("MVA_ALLOW_NETWORK=1", direct.stderr)

if __name__ == "__main__":
    unittest.main()
