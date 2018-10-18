#Calcul the Pearson correlation coefficient between the promoter and all DHS peak
#in less than 500kb.
#Usage: python Prom-other_correlation_eco.py DNase_DHS.bed DHS_promoter.bed > prom-other_correlation.txt

import argparse as ap
import scipy.stats
from math import log10

class DHS(object): #For each DHS.
    def __init__(self,fields):
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.ID = fields[3]
        self.counts = map(float,fields[4:]) #All counts are floats.

    def __str__(self):
        return ("%s\t%s\t%s\t%s" % (self.chrom,self.start,self.end,self.ID))

    @property
    def log_counts(self): #Avoid false correlation and NA values
        epsilon = 0.001
        return [log10(i+epsilon) for i in self.counts]

    @property
    def interval(self):
        return [self.chrom,self.start,self.end]

    def is_within_neighbourhood(self,DHS_interval):
        if ((self.chrom == DHS_interval[0])
                 and ((self.start > DHS_interval[2] and self.start <= DHS_interval[2]+500000)
                 or (self.end < DHS_interval[1] and self.end+500000 >= DHS_interval[1]))) or self.interval == DHS_interval:
            return True
        else:
            return False


def DHS_correlation (DHS1, DHS2):
    coefP = scipy.stats.pearsonr(DHS1.log_counts,DHS2.log_counts)[0]
    pvalueP = scipy.stats.pearsonr(DHS1.log_counts,DHS2.log_counts)[1]
    coefS = scipy.stats.spearmanr(DHS1.log_counts,DHS2.log_counts)[0]
    pvalueS = scipy.stats.spearmanr(DHS1.log_counts,DHS2.log_counts)[1]
    return str(coefP),str(pvalueP),str(coefS),str(pvalueS)


parser = ap.ArgumentParser()
parser.add_argument("file_dnase")
parser.add_argument("file_DHSprom")
args = parser.parse_args()

all_DHS = []
with open (args.file_dnase,"r") as f:
    for line in f:
        fields = line.rstrip().split("\t")
        if fields[0] == "#chr":
            header = fields
        else:
            all_DHS.append(DHS(fields))

left_border = 0
print("\t".join(["#prom_chr","prom_start","prom_end","prom_ID","chr","start",
"end","ID","Pearson_corr","p_value","Spearman_corr","p_value"])) #Header line
with open (args.file_DHSprom,"r") as f:
    for line in f:
        fields = line.split("\t")
        prom_DHS = DHS(fields)

        while not prom_DHS.is_within_neighbourhood(all_DHS[left_border].interval):
            left_border += 1

        i = left_border
        while prom_DHS.is_within_neighbourhood(all_DHS[i].interval):
            if prom_DHS.interval != all_DHS[i].interval:
                (coefP,pvalueP,coefS,pvalueS) = DHS_correlation(prom_DHS,all_DHS[i])
                print("\t".join([str(prom_DHS),str(all_DHS[i]),coefP,pvalueP,coefS,pvalueS]))
            i += 1
            if i >= len(all_DHS):
                break
