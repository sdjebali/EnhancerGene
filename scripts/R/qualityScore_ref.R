# R script: study the distribution of scores for each BEDPE file containing
#spatials interactions.
# If these files do not exist any more, just choose the initial files and write
#on the script the right column number.

# Do not forget to define the work directory
library(data.table)

## capture Hi-C
cHiC=fread("indice_cHiC.tsv", sep="\t", h=FALSE)

summary(cHiC)
# Min.    1st Qu.   Median    Mean    3rd Qu.   Max.
#6.780    9.016     9.644     9.665   10.297    14.226

par(mfrow=c(1,2))
boxplot (cHiC)
hist(cHiC[,1])


## in situ Hi-C
insitu=fread("indice_insituHiC.tsv", sep="\t", h=FALSE)
insituFDR=fread("indice_insituHiC_FDR.tsv", sep="\t", h=FALSE)

#Calculating log(observed/expected)
bl <- vector()
donut <- vector()
h <- vector()
v <- vector()

for (i in (1:9448)){ #+1 in order to avoid infinity values
  bl <- c(bl, log(insitu[i,1]/(insitu[i,2]+1)))
  donut <- c(donut, log(insitu[i,1]/(insitu[i,3]+1)))
  h <- c(h, log(insitu[i,1]/(insitu[i,4]+1)))
  v <- c(v, log(insitu[i,1]/(insitu[i,5]+1)))
}

summary(bl)
# Min.   1st Qu.  Median    Mean     3rd Qu.    Max.
#0.5358  0.7548   0.8896    0.9474    1.0630    5.8450

summary(donut)
# Min.   1st Qu.  Median  Mean    3rd Qu.  Max.
#0.5461  0.7719  0.9050   0.9642  1.0810  4.7520

summary(h)
# Min.    1st Qu.  Median    Mean     3rd Qu.    Max.
#0.3964  0.6429    0.7762    0.8306    0.9446   5.4420

summary(v)
# Min.    1st Qu.  Median    Mean     3rd Qu.    Max.
#0.3931  0.6403   0.7734     0.8178  0.9425     4.7200

summary(insituFDR)
#       bl                  donut               h                   v
#Min.   :0.000e+00   Min.   :0.000e+00   Min.   :0.0000000   Min.   :0.000e+00
#1st Qu.:0.000e+00   1st Qu.:0.000e+00   1st Qu.:0.0000000   1st Qu.:0.000e+00
#Median :1.000e-08   Median :0.000e+00   Median :0.0000036   Median :4.740e-06
#Mean   :9.780e-04   Mean   :6.941e-04   Mean   :0.0038747   Mean   :4.056e-03
#3rd Qu.:1.046e-05   3rd Qu.:5.950e-06   3rd Qu.:0.0008527   3rd Qu.:9.225e-04
#Max.   :9.397e-02   Max.   :9.005e-02   Max.   :0.0980803   Max.   :9.997e-02

par(mfrow=c(2,4))
boxplot (bl)
boxplot (donut)
boxplot (h)
boxplot (v)

hist(bl)
hist(donut)
hist(h)
hist(v)

par(mfrow=c(2,4))
hist(bl)
hist(donut)
hist(h)
hist(v)

hist(insituFDR[,1])
hist(insituFDR[,2])
hist(insituFDR[,3])
hist(insituFDR[,4])

par(mfrow=c(2,2))
plot(bl, insituFDR[,1])
plot(donut, insituFDR[,2])
plot(h, insituFDR[,3])
plot(v, insituFDR[,4])


## ChIA-PET (rad21 files)
emo=read.table("indice_EMO.tsv", sep="\t", h=FALSE)
emq=read.table("indice_EMQ.tsv", sep="\t", h=FALSE)

summary(emo)
# Min.   1st Qu.  Median    Mean    3rd Qu.    Max.
#0.2662  2.2120    2.5970   2.7280  3.1220    5.5170

summary(emq)
# Min.   1st Qu.  Median    Mean  3rd Qu.    Max.
#-1.846   2.353   2.709    2.796   3.151    5.057

par(mfrow=c(2,2))
boxplot (emo)
boxplot (emq)
hist(emo[,1])
hist(emq[,1])
