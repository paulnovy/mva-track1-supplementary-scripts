import importlib.util
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch


class NativeClientTest(unittest.TestCase):
    def test_explicit_transport_and_key_boundary_without_network(self):
        path = Path(__file__).resolve().parents[1] / "scripts/alphagenome/alphagenome_native_client.py"
        spec = importlib.util.spec_from_file_location("native_client_test", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(RuntimeError):
                module.create_client()
            os.environ["MVA_ALLOW_NETWORK"] = "1"
            with self.assertRaises(RuntimeError):
                module.create_client()
            os.environ["MVA_ALPHAGENOME_TRANSPORT"] = "native"
            for credential in ("", "oc-sent-v2.synthetic-test-sentinel"):
                os.environ["ALPHAGENOME_API_KEY"] = credential
                with self.assertRaises(RuntimeError):
                    module.create_client()
            os.environ["ALPHAGENOME_API_KEY"] = "synthetic-test-key-not-a-credential"
            fake_models = SimpleNamespace(dna_client=SimpleNamespace(create=lambda **kwargs: kwargs))
            with patch.dict(sys.modules, {"alphagenome.models": fake_models}):
                result = module.create_client(timeout_seconds=123)
            self.assertEqual(result["timeout"], 123)
            self.assertEqual(result["address"], "dns:///gdmscience.googleapis.com:443")
            self.assertEqual(result["api_key"], "synthetic-test-key-not-a-credential")


if __name__ == "__main__":
    unittest.main()
