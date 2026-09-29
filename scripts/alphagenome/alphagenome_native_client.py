"""Explicit native SDK transport for standalone reproduction, never a fallback."""
import os


def create_client(*, target="dns:///gdmscience.googleapis.com:443", timeout_seconds=900):
    if os.environ.get("MVA_ALLOW_NETWORK") != "1":
        raise RuntimeError("Set MVA_ALLOW_NETWORK=1 for online prediction")
    if os.environ.get("MVA_ALPHAGENOME_TRANSPORT") != "native":
        raise RuntimeError("Explicitly select MVA_ALPHAGENOME_TRANSPORT=native")
    key = os.environ.get("ALPHAGENOME_API_KEY", "")
    if not key or key.startswith("oc-sent-v2."):
        raise RuntimeError("Native transport requires an externally provisioned API key, not a gateway sentinel")
    from alphagenome.models import dna_client

    # SDK timeout is channel-readiness timeout, not a per-request deadline.
    return dna_client.create(api_key=key, timeout=timeout_seconds, address=target)
