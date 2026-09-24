#!/bin/bash
cd ~/mega/item01/discovery
while [ $(wc -l < spliceai_scores.tsv) -lt 700 ]; do sleep 20; done
pkill -f "spliceai_run.py sub_ids.txt"; sleep 5
python3 spliceai_run.py cand_ids.txt >> spliceai_cand.log 2>&1
echo cand_done >> spliceai_cand.log
python3 spliceai_run.py sub_ids.txt >> spliceai.log 2>&1
