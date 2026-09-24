#!/usr/bin/env python3
"""Live re-run: re-classify the frozen 500-variant ClinVar sample against the
CURRENT ClinVar via NCBI eutils, proving the frozen labels still match the
live database (drift check). Rate-limited; not part of CI."""
import json, sys, time, urllib.request, urllib.parse

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def esummary(variation_id):
    q = urllib.parse.urlencode({"db": "clinvar", "id": variation_id,
                                "retmode": "json"})
    with urllib.request.urlopen(BASE + "esummary.fcgi?" + q, timeout=30) as r:
        return json.loads(r.read())


def main(sample_path, out_path, limit=None):
    import csv
    rows = list(csv.DictReader(open(sample_path), delimiter="\t",
                               fieldnames=["variation_id", "gene", "hgvs",
                                           "significance", "review", "chrom",
                                           "pos", "ref", "alt"]))
    if limit:
        rows = rows[:limit]
    drift = []
    for i, r in enumerate(rows):
        try:
            d = esummary(r["variation_id"])
            doc = d["result"][r["variation_id"]]
            live = doc.get("germline_classification", {}).get("description", "")
            if live and live != r["significance"]:
                drift.append({"id": r["variation_id"],
                              "frozen": r["significance"], "live": live})
        except Exception as e:
            drift.append({"id": r["variation_id"], "error": str(e)})
        if i % 3 == 2:
            time.sleep(0.4)  # stay under 3 req/s without an API key
    json.dump({"checked": len(rows), "drift": drift},
              open(out_path, "w"), indent=2)
    print(f"checked {len(rows)}, drift {len(drift)}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else None)
