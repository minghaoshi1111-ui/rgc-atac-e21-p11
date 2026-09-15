#!/usr/bin/env python3

import json, os, sys, requests

ROOT = "/mnt/d/ngs_rebuild"
os.chdir(ROOT)
os.makedirs("results", exist_ok=True)
URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

def read(p):
    with open(p) as fh:
        return [l.strip() for l in fh if l.strip()]

bg = read("results/gene_lists/background_genes.txt")
print(f"background: {len(bg)} genes")

COLS = ["source","native","name","p_value","term_size","query_size",
        "intersection_size","precision","recall"]

for stage in ("E21_up", "P11_up"):
    q = read(f"results/gene_lists/{stage}_genes.txt")
    print(f"{stage}: querying {len(q)} genes ...", flush=True)
    r = requests.post(URL, json={
        "organism": "rnorvegicus",
        "query": q,
        "background": bg,
        "domain_scope": "custom",
        "sources": ["GO:BP"],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "no_evidences": True,
    }, headers={"User-Agent": "rgc-atac-rebuild"}, timeout=600)
    r.raise_for_status()
    hits = r.json().get("result", [])
    out = f"results/GO_{stage}.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for h in sorted(hits, key=lambda x: x["p_value"]):
            fh.write("\t".join(str(h.get(c, "")) for c in COLS) + "\n")
    print(f"  {len(hits)} enriched GO:BP terms -> {out}")
