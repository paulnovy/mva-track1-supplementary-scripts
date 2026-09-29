"""Offline replay using invented ontology entries and no case data."""

import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/exomiser/prepare_mva1_phenotype_sensitivity.py"


class PhenotypeSensitivityTest(unittest.TestCase):
    def test_explicit_private_config_drives_profiles_and_audit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            docx = root / "synthetic.docx"
            with zipfile.ZipFile(docx, "w") as archive:
                archive.writestr(
                    "word/document.xml",
                    '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
                    '<w:body><w:p><w:r><w:t>Synthetic source A.</w:t></w:r></w:p>'
                    '<w:p><w:r><w:t>Synthetic source B.</w:t></w:r></w:p>'
                    '<w:p><w:r><w:t>Synthetic relative context.</w:t></w:r></w:p>'
                    '</w:body></w:document>',
                )
            hpo = root / "synthetic.obo"
            hpo.write_text(
                "format-version: 1.2\n\n[Term]\nid: HP:9999997\nname: Synthetic A\n"
                "\n[Term]\nid: HP:9999998\nname: Synthetic B\n"
                "\n[Term]\nid: HP:9999999\nname: Synthetic relative feature\n"
            )
            vcf = root / "synthetic.vcf.gz"
            with gzip.open(vcf, "wt") as handle:
                handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n")
            curation = []
            for number, (term, label) in enumerate(
                (("HP:9999997", "Synthetic A"), ("HP:9999998", "Synthetic B"),
                 ("HP:9999999", "Synthetic relative feature")), 1
            ):
                curation.append({
                    "key": f"synthetic_{number}", "hpo_id": term, "label": label,
                    "feature_paragraphs": [number], "onset": "not stated",
                    "subject": "proband" if number < 3 else "family history",
                    "justification": "Synthetic test input only.",
                })
            payload = {
                "curation": curation,
                "profiles": {
                    "baseline": ["HP:9999997"],
                    "full": ["HP:9999998", "HP:9999997"],
                    "axis": ["HP:9999998"],
                },
                "family_history": {"synthetic_context": True},
            }
            config = root / "private.json"
            config.write_text(json.dumps(payload))
            output = root / "outputs"
            command = [sys.executable, str(SCRIPT), "--docx", str(docx), "--hpo", str(hpo),
                       "--vcf", str(vcf), "--output", str(output)]
            missing = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(missing.returncode, 0)
            self.assertIn("--config", missing.stderr)
            self.assertFalse(output.exists())
            subprocess.run(command + ["--config", str(config)], check=True, capture_output=True, text=True)

            for profile, expected_terms in payload["profiles"].items():
                packet_path = output / profile / "phenopacket.json"
                packet = json.loads(packet_path.read_text())
                self.assertEqual([item["type"]["id"] for item in packet["phenotypicFeatures"]], expected_terms)
                self.assertEqual(packet["subject"]["id"], "SYNTHETIC")
                self.assertEqual(packet_path.stat().st_mode & 0o777, 0o600)
                self.assertEqual((output / profile / "family.ped").read_text(), "FAMILY\tSYNTHETIC\t0\t0\t0\t2\n")
            self.assertEqual(json.loads((output / "family-history.json").read_text()), payload["family_history"])
            provenance = json.loads((output / "provenance.json").read_text())
            self.assertEqual(provenance["curation_config"]["sha256"], hashlib.sha256(config.read_bytes()).hexdigest())
            with (output / "curation.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[1]["original_wording"], "Synthetic source B.")
            self.assertEqual(rows[1]["scored_profiles"], "full,axis")
            self.assertEqual(rows[2]["scored_profiles"], "none_interpretive_input")


if __name__ == "__main__":
    unittest.main()
