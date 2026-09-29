#!/usr/bin/env python3
"""Query bounded gnomAD co-occurrence data, or replay a saved JSON response."""
import argparse
import importlib.util
import json
from pathlib import Path

QUERY = """
query($variants: [String!]!, $dataset: DatasetId!) {
  variant_cooccurrence(variants: $variants, dataset: $dataset) {
    variant_ids genotype_counts haplotype_counts p_compound_heterozygous
    populations { id genotype_counts haplotype_counts p_compound_heterozygous }
  }
}
"""

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", nargs="+", required=True, help="gnomAD variant IDs")
    ap.add_argument("--dataset", default="gnomad_r2_1")
    ap.add_argument("--output", type=Path, required=True)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--saved-response", type=Path, help="offline replay input")
    mode.add_argument("--wrapper", type=Path, help="opt-in local gnomAD client wrapper")
    args = ap.parse_args()
    if args.wrapper and __import__('os').environ.get('MVA_ALLOW_NETWORK') != '1':
        ap.error('set MVA_ALLOW_NETWORK=1 for opt-in gnomAD query')
    if args.saved_response:
        response = json.loads(args.saved_response.read_text())
    else:
        spec = importlib.util.spec_from_file_location("gnomad_wrapper", args.wrapper)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        response = module.CLIENT.fetch_json(
            "https://gnomad.broadinstitute.org/api", method="POST",
            json_body={"query": QUERY, "variables": {"variants": args.variants, "dataset": args.dataset}},
        )
    if response.get("errors"):
        raise SystemExit("gnomAD response contains GraphQL errors; raw response was not interpreted")
    cooccurrence = response.get("data", {}).get("variant_cooccurrence")
    summary = {
        "dataset": args.dataset,
        "requested_variant_ids": args.variants,
        "available": cooccurrence is not None,
        "variant_ids": None if cooccurrence is None else cooccurrence.get("variant_ids"),
        "genotype_counts": None if cooccurrence is None else cooccurrence.get("genotype_counts"),
        "haplotype_counts": None if cooccurrence is None else cooccurrence.get("haplotype_counts"),
        "p_compound_heterozygous": None if cooccurrence is None else cooccurrence.get("p_compound_heterozygous"),
        "populations": [] if cooccurrence is None else cooccurrence.get("populations", []),
        "limitations": ["population co-occurrence is not sample phase", "missing data are not zero"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"summary": summary, "raw_response": response}, indent=2, sort_keys=True) + "\n")

if __name__ == "__main__":
    main()
