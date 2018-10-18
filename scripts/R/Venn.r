# Careful with the file names et number ! This was wroten with ChIA-PET RAD21,
#but it is easy to convert to three files only.

#Setting work directory
#!/usr/local/bioinfo/src/R/R-3.3.3/bin/R
library("futile.logger")
library("VennDiagram")

insitu=unlist(read.delim("liste_insitu.txt", sep="\t", h=FALSE))
cHic=unlist(read.delim("liste_cHic.txt", sep="\t", h=FALSE))
emo=unlist(read.delim("liste_chiapetEMO.txt", sep="\t", h=FALSE))
emq=unlist(read.delim("liste_chiapetEMQ.txt", sep="\t", h=FALSE))

venn.plot <- venn.diagram(
x=list( #input data
A=insitu,
B=cHic,
C=emo,
D=emq
),
filename = NULL,
fill=c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4"), #field color
alpha=0.50,
category = c("insitu_HiC","cHiC","ChIA-PET_EMO","ChIA-PET_EMQ"), #legend
fontfamily="arial",
main.cex=2,
main.fontfamily="arial",
sub.cex=1.3,
sub.fontfamily="arial",
cat.fontfamily=c("arial","arial","arial","arial") #font of legend
)
grid.draw(venn.plot)
