#!/bin/bash

wd=tmp/001
mkdir -p ${wd}

S=17

etc/mkAluIns ${S} etc/BRCA2-ex3.bed etc/some-alus.fa > ${wd}/var.txt 2> ${wd}/exp.txt
v=$(cat ${wd}/var.txt)
msg="homozygous mut insAlu ${v}"

zot spikein -g ~/data/hg19 -b etc/BRCA2-exons.bed -N 50000 -S $((${S} + 1)) -V 1.0 -f ${wd}/var.txt -z ${wd}/reads

zot alu-finder -g ~/data/hg19 etc/BRCA2-exons.bed ${wd}/reads*.fastq.gz > ${wd}/output.txt

ac=$(cat ${wd}/exp.txt | cut -f 3)
p0=$(cat ${wd}/exp.txt | cut -f 4)
p1=$(($p0 + 1))
se=$(cat ${wd}/exp.txt | cut -f 5)

cat > ${wd}/rules.py << EOF
- "accession == '${ac}'"
- "int(after) == ${p0}"
- "int(before) == ${p1}"
- "'${se}'.startswith(lhsSeq)"
- "'${se}'.endswith(rhsSeq)"
EOF
./bin/check ${wd}/rules.py ${wd}/output.txt > ${wd}/result.txt

echo "$0" $(cat ${wd}/result.txt) ${msg}
