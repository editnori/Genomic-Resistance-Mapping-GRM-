#!/usr/bin/env bash
# KIT10 round 2: PORT's synthetic master root for drv-gates-v4.sh (stub-regress), made on the test box with no data
# from anywhere. PORT's gen-root.py writes fasta/, bvbrc/ and catalog/ (stdlib, fixed seeds); a synthetic layer adds
# what stage-manifest.py reads besides: species/campylobacter_jejuni/genomes.txt, curated/, species-check/ and the
# revision-9 canary (rev9/report/index.json and one dataset). The manifest lists every file. Run 1's master manifest
# (3b178c6f...) came from PORT's own build, which also ran grm, so this one isn't expected to equal it; the last line
# says whether it does.
# usage: gen-port-master.sh GEN_ROOT_PY OUT_DIR MANIFEST
set -eu
G=${1:?gen-root.py} O=${2:?out dir} MAN=${3:?manifest}
PIN=3b178c6f2d95731e786b52188f0e07ee3405e1a15eed7180c8fc17ad5dac7720
[ ! -e "$O" ] || { echo "STOP: $O is already there"; exit 1; }
python3 -I -B "$G" "$O" 'Campylobacter jejuni:60' --length 30000 --clades 6 --seed 20261105
S=campylobacter_jejuni
mkdir -p "$O/species/$S" "$O/curated" "$O/species-check" "$O/rev9/report/datasets"
python3 -I -B - "$O" "$S" << 'EOF'
import json, os, sys
o, s = sys.argv[1:]
ids = sorted(json.loads(l)["genome_id"] for l in open(os.path.join(o, "bvbrc/genomes.jsonl")))
open(os.path.join(o, "species", s, "genomes.txt"), "w").write("\n".join(ids) + "\n")
ds = [{"dataset": "%s__%s" % (s, d), "species": s, "drug": d, "genomes": len(ids), "synthetic": True} for d in ("amikacin", "ciprofloxacin")]
json.dump(ds, open(os.path.join(o, "curated/datasets.json"), "w"), indent=1, sort_keys=True)
json.dump({"synthetic": True, "species": {s: {"checked": len(ids), "fail": 0}}}, open(os.path.join(o, "species-check/species_check.json"), "w"), indent=1, sort_keys=True)
json.dump({"synthetic": True, "revision": 9, "datasets": [d["dataset"] for d in ds]}, open(os.path.join(o, "rev9/report/index.json"), "w"), indent=1, sort_keys=True)
for d in ds:
    json.dump({"synthetic": True, "dataset": d["dataset"], "metrics": {}}, open(os.path.join(o, "rev9/report/datasets", d["dataset"] + ".json"), "w"), indent=1, sort_keys=True)
EOF
(cd "$O" && find . -type f | sed 's#^\./##' | LC_ALL=C sort | xargs sha256sum) > "$MAN"
m=$(sha256sum < "$MAN" | cut -d' ' -f1)
echo "master: $(wc -l < "$MAN") files, $(du -sb "$O" | cut -f1) bytes, manifest $m"
[ "$m" = "$PIN" ] && echo "MASTER REPRODUCED: equals run 1's $PIN" || echo "MASTER NOT REPRODUCED: run 1's was $PIN; this run uses the fresh synthetic master above"
