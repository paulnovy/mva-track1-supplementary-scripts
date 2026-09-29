import json
import os
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/copy_number/run_mocha_wgs.sh"


class RunMochaWgsArgvTest(unittest.TestCase):
    def test_excluded_sites_are_passed_as_mocha_exclusion(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            input_bcf = work / "input.bcf"
            sites_bcf = work / "excluded.bcf"
            cnps_bed = work / "cnps.bed"
            for path in (input_bcf, sites_bcf, cnps_bed):
                path.write_text("synthetic\n")
            (work / "input.bcf.csi").write_text("index\n")
            (work / "excluded.bcf.csi").write_text("index\n")
            output = work / "out"
            captured = work / "docker-argv.json"
            fake_bin = work / "bin"
            fake_bin.mkdir()
            fake_docker = fake_bin / "docker"
            fake_docker.write_text(
                "#!/bin/sh\n"
                "python3 -c 'import json,os,sys; json.dump(sys.argv[1:], open(os.environ[\"CAPTURE\"], \"w\"))' \"$@\"\n"
            )
            fake_docker.chmod(fake_docker.stat().st_mode | stat.S_IXUSR)

            env = os.environ.copy()
            env["PATH"] = str(fake_bin) + os.pathsep + env["PATH"]
            env["CAPTURE"] = str(captured)
            subprocess.run(
                [str(SCRIPT), str(input_bcf), str(sites_bcf), str(cnps_bed), str(output), "SYNTH"],
                check=True,
                env=env,
            )

            argv = json.loads(captured.read_text())
            self.assertIn("^/excluded/excluded.bcf", argv)
            self.assertNotIn("/excluded/excluded.bcf", argv)
            self.assertIn("mva-mocha:bcftools-1.20-95686b7", argv)
            self.assertFalse((output / "annotated.bcf").exists())


if __name__ == "__main__":
    unittest.main()
