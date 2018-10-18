#!/bin/bash
# See the coverage of each connected component.
# I separate the initial file by chromosome.

rm coverage_results_ID_score_chr22
awk 'BEGIN{OFS="\t"} $1=="chr22"{print $1,$2,$3,$4,$5,NR}' coverage_fichierA_intervalles_ID_score_corrige.bed | while read chr beg end id score no
do
  rm tmpAchr_*.bed tmpBchr_*.bed
  printf "$chr\t$beg\t$end\t$id\t$score\n" > tmpAchr_$id.bed
  awk -v no=$no 'NR==no{print $0}' composante_connexes_seuil12000_intervalles | sed 's/\t/\n/g' | sed 's/;/\t/g' | awk 'BEGIN{OFS="\t"} {print $2,$3,$4 "\n" $5,$6,$7}' > tmpBchr_$id.bed
  bedtools coverage -a tmpAchr_$id.bed -b tmpBchr_$id.bed >> coverage_results_ID_score_chr22
done
