#!/bin/bash
# This script generate all the restriction fragments of HindIII enzyme.
# It uses the tool hicup_v0.5.7.
# This script was written by Sylvain.

enzyme="HindIII"
site="A^AGCTT" #restriction site
tmp=`date +%F_%H-%M-%S`; #for the temporary file.

for species in homo_sapiens ; do #If we want to generate the file from different annotation.
version=hg19.gencv19;
echo "#Digesting $version with $enzyme...";
fasta=/work/project/fragencode/data/species/$species/$version/$species.fa;
output=/work/project/dynagen/cmestre/restrictionsfragments_"$enzyme".bed;
tmpdir=/work/project/dynagen/cmestre/tmp$tmp;
mkdir -p $tmpdir;
/work/project/dynagen/cmestre/tools/hicup_v0.5.7/hicup_digester --genome $species --re1 $site,$enzyme --outdir $tmpdir $fasta;
cat $tmpdir/Digest_"$species"_"$enzyme"_*.txt| awk -v OFS="\t" '$1!=prev{n=0}NR>2{print $1,$2-1,$3,$13,$1,$1"_"++n,"0","+"}{pre$1}' > $output && echo "## Done in $output" && rm -Rf $tmpdir;
done
