#!/bin/bash
# Create file for each sub-group of the Venn diagramm.
# Careful with the files names and the writing of 'chiapet', 'chic' and
#'insitu' !

#Separate the connected component of each subgroup of the Venn diagramm.
grep 'chiapet' compconnexes_intervalles.txt > tmp
grep 'insitu' tmp > tmpbis
grep 'chic' tmpbis > compconnexes_intervalles_centre.txt
grep -v 'chic' tmpbis > compconnexes_intervalles_ChiapetInsitu.txt
grep 'chic' tmp | grep -v 'insitu' > compconnexes_intervalles_ChiapetChic.txt
grep -v 'insitu' tmp | grep -v 'chic' > compconnexes_intervalles_singletonChiapet.txt

grep 'insitu' compconnexes_intervalles.txt > tmp2
grep -v 'chiapet' tmp2 > tmp2bis
grep -v 'chic' tmp2bis > compconnexes_intervalles_singletonInsitu.txt
grep 'chic' tmp2bis > compconnexes_intervalles_InsituChic.txt

grep 'chic' compconnexes_intervalles.txt > tmp3
grep -v 'insitu' tmp3 | grep -v 'chiapet' > compconnexes_intervalles_singletonChic.txt

rm tmp tmpbis tmp2 tmp2bis tmp3


#Create one consensus interval for each connected component.
python merge_intevals.py compconnexes_intervalles_centre.txt > Venn_nx_centre.txt
python merge_intevals.py compconnexes_intervalles_ChiapetInsitu.txt > Venn_nx_ChiapetInsitu.txt
python merge_intevals.py compconnexes_intervalles_ChiapetChic.txt > Venn_nx_ChiapetChic.txt
python merge_intevals.py compconnexes_intervalles_singletonChiapet.txt > Venn_nx_singletonChiapet.txt

python merge_intevals.py compconnexes_intervalles_singletonInsitu.txt > Venn_nx_singletonInsitu.txt
python merge_intevals.py compconnexes_intervalles_InsituChic.txt > Venn_nx_InsituChic.txt

python merge_intevals.py compconnexes_intervalles_singletonChic.txt > Venn_nx_singletonChic.txt
