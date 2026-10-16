"""Export the ladder orbits shown on the findings page to data/ladder_orbits.json and print a
Markdown table for the README. Run after web/build_data.js."""
import json, os
here = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(here, "web/ladder_data.json")))
checks = {}
for f in ["results/claims_extra_check.json"]:
    p = os.path.join(here, f)
    if os.path.exists(p):
        for r in json.load(open(p)): checks[r["id"]] = r["res_max"]
raw = json.load(open(os.path.join(here, "results/verification_raw.json")))["result"]["results"][0]["claims"]
import re
for c in raw:
    m = re.search(r"Residual ([0-9.]+e-[0-9]+)", c["detail"]) or re.search(r"max\|Y\(T\)-Y0\| = ([0-9.]+e-[0-9]+)", c["detail"])
    if m: checks[c["id"]] = float(m.group(1))
rows = []
for o in sorted(D["orbits"], key=lambda o: (o["m3"], o["k"], o.get("extra", False))):
    k = o["k"]
    rows.append({
        "id": o["id"], "m3": o["m3"], "k": k, "label": "3'" if o.get("extra") else str(k),
        "v1": o["v1"], "v2": o["v2"], "T": o["T"], "E": o["E"], "Tstar": o["Tstar"],
        "word": o["canon"], "word_compact": f"(Ba)^{k+1}(bA)^{k+1}", "letters": len(o["canon"]),
        "status": o["status"], "note": o.get("note"), "match": o.get("match"),
        "closure_this_code": o["closure"], "closure_independent": checks.get(o["id"]),
        "max_floquet_multiplier": o.get("maxMult"),
    })
# 3' (X3@0.25) crosses the slice again at half period at scale 0.8905: the same orbit from a second start
for r in rows:
    if r["id"] == "X3@0.25":
        r["alternate_start"] = {"v1": 0.15864128667707048, "v2": 0.10037530718943394, "T": 20.768761004887477,
                                "note": "same orbit, reached at T/2 rescaled by 0.890466; count once"}
out = {"description": "Periodic orbits of the ladder family (Ba)^(k+1)(bA)^(k+1) in the isosceles slice m1 = m2 = 1, r1 = (-1,0), r2 = (1,0), r3 = (0,0), v1 = v2 = (v1, v2), v3 = -2 (v1, v2)/m3, G = 1. T* = T |E|^(3/2).",
       "orbits": rows}
json.dump(out, open(os.path.join(here, "data/ladder_orbits.json"), "w"), indent=1)
print("| rung k | m₃ | status | v₁ | v₂ | T | T* | closure (independent) | growth / period |")
print("|---|---|---|---|---|---|---|---|---|")
for r in rows:
    st = "**new**" if r["status"] == "new" else ("catalogued" if (r["match"] or "").startswith("Li") else "known word, new mass")
    ci = f"{r['closure_independent']:.1e}" if r["closure_independent"] else "–"
    mm = f"×{r['max_floquet_multiplier']:.1f}" if r["max_floquet_multiplier"] else "–"
    print(f"| {r['label']} | {r['m3']} | {st} | {r['v1']:.10f} | {r['v2']:.10f} | {r['T']:.8f} | {r['Tstar']:.3f} | {ci} | {mm} |")
