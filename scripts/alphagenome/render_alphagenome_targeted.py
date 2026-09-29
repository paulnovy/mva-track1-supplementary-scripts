#!/usr/bin/env python3
"""Render protected exact-allele AlphaGenome Atlas results."""

import argparse
import json
from pathlib import Path


FEATURE_NAMES = [
    "MERGED_SPLICING", "MAX_ABS_ATAC", "MAX_ABS_CONTACT_MAPS",
    "MAX_ABS_DNASE", "MAX_ABS_CHIP_TF", "MAX_ABS_CHIP_HISTONE",
    "MAX_ABS_CAGE", "MAX_ABS_PROCAP", "MAX_ABS_RNA_SEQ",
    "MAX_ABS_POLYADENYLATION", "ALPHAMISSENSE", "CACTUS_241_WAY",
    "PROTEIN_TERMINATION", "START_LOST", "STOP_LOST",
    "PHASTCONS_470_WAY", "IS_INSERTION", "IS_DELETION",
]


def rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text().splitlines() if line]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--evidence", type=Path, required=True)
    ap.add_argument("--lookup", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.out.chmod(0o700)
    ev = json.loads(args.evidence.read_text())
    wanted = []
    for i, item in enumerate(ev["candidates"], 1):
        chrom, pos = item["coordinate"].rsplit(":", 1)
        wanted.append({"id": f"BUB1B-{i}", "chrom": chrom.removeprefix("chr"),
                       "pos": int(pos), "ref": item["ref"], "alt": item["alt"][0]})

    table_rows = {n: rows(args.lookup / f"{n}.exact.tsv")
                  for n in ("splicing", "avi", "feature_importance")}
    results = []
    for item in wanted:
        key = (item["chrom"], item["pos"], item["ref"], item["alt"])
        matched = {}
        for name, data in table_rows.items():
            hits = [r for r in data if (r[0].removeprefix("chr"), int(r[1]), r[2], r[3]) == key]
            if len(hits) != 1:
                raise RuntimeError(f"{name}: expected one exact allele match, got {len(hits)}")
            matched[name] = hits[0]
        features = {name: float(value) for name, value in zip(FEATURE_NAMES, matched["feature_importance"][4:])}
        ranked = sorted(features.items(), key=lambda kv: abs(kv[1]), reverse=True)
        results.append({
            **item,
            "coordinate": f"chr{item['chrom']}:{item['pos']}",
            "exact_match": {"splicing": True, "avi": True, "feature_importance": True},
            "alphagenome_splicing": float(matched["splicing"][4]),
            "splicing_gt_1_large_effect_heuristic": float(matched["splicing"][4]) > 1.0,
            "avi_raw_score": float(matched["avi"][4]),
            "avi_phred": float(matched["avi"][5]),
            "avi_feature_importance": features,
            "top_absolute_feature_contributions": [
                {"feature": name, "value": value} for name, value in ranked[:5]
            ],
        })

    manifest = json.loads((args.lookup / "lookup-manifest.json").read_text())
    summary = {
        "schema_version": 1,
        "status": "complete",
        "scope": "local exact-allele AlphaGenome Atlas annotation of two BUB1B SNVs",
        "reference_build": "GRCh38",
        "privacy": "local-only; no patient allele or phenotype sent to an external service",
        "candidate_count": len(results),
        "exact_match_counts": {"splicing": len(results), "avi": len(results), "feature_importance": len(results)},
        "results": results,
        "interpretation": {
            "any_splicing_gt_1": any(x["splicing_gt_1_large_effect_heuristic"] for x in results),
            "statement": "Both exact SNVs are present. Neither exceeds the provisional >1 merged-splicing large-effect heuristic. Elevated AVI prioritization is dominated by protein-termination/AlphaMissense and conservation features rather than a dominant AlphaGenome regulatory or merged-splicing contribution; all scores remain predictive and non-diagnostic.",
            "lztr1_queried": False,
            "lztr1_reason": "Not needed to answer the exact BUB1B second-allele mechanism question; avoiding scope expansion.",
            "shap_wholesale_extracted": False,
        },
        "provenance": manifest,
        "references": [
            "https://deepmind.google/science/alphagenome/",
            "https://www.alphagenomedocs.com/colabs/splicing_variant_scoring.html",
            "https://github.com/google-deepmind/alphagenome/blob/main/docs/source/variant_scoring.md",
        ],
        "limitations": [
            "Atlas scores are in-silico molecular-effect predictions and are not clinical classifications, diagnoses, phase evidence, or functional validation.",
            "Absence or a low score would not establish benignity; this lookup only supports SNVs represented by the supplied Atlas tables.",
            "The provisional splicing >1 heuristic is not a validated clinical cutoff.",
            "The exact Atlas build/REF/ALT were preserved and required for matching; classification was not transferred by protein consequence.",
            "Tabix validated the supplied indices and BGZF blocks accessed by all three exact regional queries. Whole ZIP-member CRCs were not recomputed because no member was extracted; archive-declared CRC32 values are recorded in provenance.",
        ],
    }
    summary_path = args.out / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    summary_path.chmod(0o600)

    lines = [
        "# Ukierunkowana lokalna adnotacja AlphaGenome Atlas — BUB1B",
        "",
        "**Status:** zakończono. **Zakres:** dwa dokładne allele SNV BUB1B, GRCh38. Analiza wyłącznie lokalna; żadnych danych pacjenta nie wysłano do API ani serwisu zewnętrznego.",
        "",
        "## Wynik dokładnego dopasowania",
        "",
        "| Kandydat | Dokładny allel Atlas | AlphaGenome splicing | AVI raw | AVI PHRED | Największe bezwzględne wkłady AVI |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for x in results:
        top = ", ".join(f"{v['feature']}={v['value']:.6g}" for v in x["top_absolute_feature_contributions"][:3])
        lines.append(f"| {x['id']} ({x['coordinate']} {x['ref']}>{x['alt']}) | tak | {x['alphagenome_splicing']:.6g} | {x['avi_raw_score']:.6g} | {x['avi_phred']:.6g} | {top} |")
    lines += [
        "",
        "Oba dokładne allele występują we wszystkich trzech lokalnych tabelach. Żaden wynik `alphagenome_splicing` nie przekracza roboczej, nieklinicznej heurystyki >1, którą dokumentacja opisuje jako poziom zwykle związany z dużym efektem. To nie jest dowód braku wpływu na splicing ani łagodności allelu.",
        "",
        "AVI nadaje obu allelom podwyższone pozycje priorytetyzacyjne, dlatego wykonano wyłącznie dokładne zapytanie do tabeli wkładów cech. Wkłady wskazują, które cechy modelu napędzają wynik; nie dowodzą mechanizmu biologicznego. Szczegółowe wartości i komplet wkładów są w `summary.json`.",
        "",
        "Dominujące wkłady AVI są zgodne przede wszystkim z konsekwencją kodującą (odpowiednio `PROTEIN_TERMINATION` oraz `ALPHAMISSENSE`) i konserwacją, a nie z dominującym sygnałem regulatorowym lub merged-splicing AlphaGenome. Atlas nie dostarcza więc tutaj nowego, silnego wsparcia dla hipotezy regulatorowego/splicingowego mechanizmu drugiego allelu ponad jego istniejącą konsekwencję kodującą.",
        "",
        "## Znaczenie dla hipotezy",
        "",
        "Wynik wspiera dalsze traktowanie obu dokładnych alleli BUB1B jako kandydatów do ukierunkowanej walidacji, ale nie rozstrzyga patogenności drugiego allelu, cis/trans, rozpoznania MVA1 ani synergii z LZTR1. LZTR1 nie był zapytany, ponieważ nie był potrzebny do odpowiedzi o dokładnych allelach BUB1B i poszerzałby zakres bez rozstrzygającej korzyści.",
        "",
        "Priorytet dalszej walidacji pozostaje: faza/segregacja oraz eksperymentalna ocena funkcjonalna dokładnego allelu; sygnał AVI może służyć jedynie do ustalenia kolejności takich badań. Ten wynik sam nie uzasadnia nowego ogólnego przebiegu SV/CNV.",
        "",
        "## Metoda i integralność",
        "",
        "Użyto dostarczonych, tylko-do-odczytu archiwów ZIP_STORED. Człony BGZF udostępniono bez kopiowania przez tymczasowe urządzenia loop, a `tabix` 1.20 wykonał zapytania dokładnych pozycji przy użyciu dostarczonych indeksów `.tbi`. Wszystkie odczytane bloki BGZF zostały poprawnie zdekompresowane; indeksy zawierały 24 kontigi. Nie rozpakowano całych tabel splicing, AVI ani feature importance. Pełnego CRC całych członów nie obliczano, bo członów nie wyodrębniano; deklarowane CRC32 zachowano w manifeście.",
        "",
        "## Ograniczenia naukowe",
        "",
        "- AlphaGenome Atlas jest predykcją in silico, nie walidacją funkcjonalną ani klasyfikacją kliniczną.",
        "- Heurystyka splicing >1 nie jest klinicznie zwalidowanym progiem.",
        "- Niski lub brakujący wynik Atlas nie pozwala uznać allelu za łagodny; tabele obejmują tu SNV, nie wszystkie klasy wariantów.",
        "- Zachowano zgodność GRCh38 oraz dokładnych REF/ALT; nie przenoszono klasyfikacji po samej zmianie białkowej.",
        "",
        "## Źródła publiczne",
        "",
        "- [AlphaGenome / Atlas — Google DeepMind](https://deepmind.google/science/alphagenome/)",
        "- [Oficjalny tutorial oceny splicingu](https://www.alphagenomedocs.com/colabs/splicing_variant_scoring.html)",
        "- [Oficjalny opis variant scoring](https://github.com/google-deepmind/alphagenome/blob/main/docs/source/variant_scoring.md)",
    ]
    report_path = args.out / "report.md"
    report_path.write_text("\n".join(lines) + "\n")
    report_path.chmod(0o600)


if __name__ == "__main__":
    main()
