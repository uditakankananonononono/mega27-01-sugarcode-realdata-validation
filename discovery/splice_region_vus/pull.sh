#!/bin/bash
curl -sL https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz | gunzip -c | awk -F'\t' 'NR==1{print > "header.txt"; next} $17=="GRCh38" && $2=="single nucleotide variant" && $3 ~ /\(.+\):c\.[0-9]+[+-]([3-9]|1[0-9]|20)[ACGT]>[ACGT]/ {print $31"\t"$5"\t"$3"\t"$7"\t"$25"\t"$19"\t"$32"\t"$33"\t"$34}' > splice_region_snvs.tsv
echo done > pull.done
