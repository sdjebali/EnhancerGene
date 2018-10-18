#Calcul the Pearson correlation coefficient between the promoter and all DHS peak
#in less than 500kb.
#Usage: python Prom-other_correlation_eco.py DNase_DHS.bed DHS_promoter.bed > prom-other_correlation.txt

import argparse as ap
import scipy.stats
from math import log10
import numpy as np
import sys
from threading import Thread, RLock

class DHS(object): #For each DHS.
    def __init__(self,fields):
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.ID = fields[3]
        self.counts = map(float,fields[4:]) #All counts are floats.
        self.sigma = np.std(self.log_counts)
        self.moy = np.mean(self.log_counts)
        self.log_counts_prime = np.array([i-self.moy for i in self.log_counts])

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


class Calcul(Thread):
     #Thread charge de lancer le calcul de correlation.
    def __init__(self, DHS1, DHS2):
        Thread.__init__(self)
        self.DHS1=DHS1
        self.DHS2=DHS2

    def run(self):
        #Code a executer pendant l'execution du thread.
        (coefP,coefS) = DHS_correlation(self.DHS1,self.DHS2)
        with verrou:
            sys.stdout.write("\t".join([str(self.DHS1),str(self.DHS2),coefP,coefS,"\n"]))
            sys.stdout.flush()


def DHS_correlation (DHS1, DHS2):
    numerator = np.dot(DHS1.log_counts_prime,DHS2.log_counts_prime)
    coefP = numerator/(10*DHS1.sigma*DHS2.sigma)
    coefS = scipy.stats.spearmanr(DHS1.log_counts,DHS2.log_counts)[0]
    return str(coefP),str(coefS)


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

verrou = RLock()
left_border = 0
print("\t".join(["#prom_chr","prom_start","prom_end","prom_ID","chr","start",
"end","ID","Pearson_corr","Spearman_corr"])) #Header line
with open (args.file_DHSprom,"r") as f:
    dejavu=[]
    for line in f:
        fields = line.split("\t")
        prom_DHS = DHS(fields)

        while not prom_DHS.is_within_neighbourhood(all_DHS[left_border].interval):
            left_border += 1

        i = left_border
        while prom_DHS.is_within_neighbourhood(all_DHS[i].interval):
            if prom_DHS.interval != all_DHS[i].interval and [all_DHS[i],prom_DHS] not in dejavu:
                thread = Calcul(prom_DHS,all_DHS[i]) # Creation des threads
                thread.start() # Lancement des threads
                thread.join() # Attend que les threads se terminent
                dejavu.append([prom_DHS,all_DHS[i]])
            i += 1
            if i >= len(all_DHS):
                break
