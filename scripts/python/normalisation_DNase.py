#Normalize each cell by the column sum.
#Usage: python normalisation.py DNase_data.bed > DNase_data_normalized.bed

import argparse as ap

class Line(object):
    def __init__(self,fields):
        self.chrom = fields[0]
        self.start = fields[1]
        self.end = fields[2]
        self.small_intestine = fields[3]
        self.adrenal_gland = fields[4]
        self.HepG2 = fields[5]
        self.stomach = fields[6]
        self.heart = fields[7]
        self.thymus = fields[8]
        self.K562 = fields[9]
        self.GM12878 = fields[10]
        self.H1_hESC = fields[11]
        self.IMR_90 = fields[12]
    def __str__(self):
        return("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
        (self.chrom,self.start,self.end,self.small_intestine,self.adrenal_gland,
        self.HepG2,self.stomach,self.heart,self.thymus,self.K562,self.GM12878,
        self.H1_hESC,self.IMR_90))


parser = ap.ArgumentParser()
parser.add_argument("file_dnase")
args = parser.parse_args()

output = []
with open (args.file_dnase,"r") as f:
    for line in f:
        fields = line.split("\t")
        if fields[0] == "chr":
            header = fields
        else:
            output.append(Line(fields))

small_intestine = sum([float(i.small_intestine) for i in output])
adrenal_gland = sum([float(i.adrenal_gland) for i in output])
HepG2 = sum([float(i.HepG2) for i in output])
stomach = sum([float(i.stomach) for i in output])
heart = sum([float(i.heart) for i in output])
thymus = sum([float(i.thymus) for i in output])
K562 = sum([float(i.K562) for i in output])
GM12878 = sum([float(i.GM12878) for i in output])
H1_hESC = sum([float(i.H1_hESC) for i in output])
IMR_90 = sum([float(i.IMR_90) for i in output])

print("\t".join(header)[:-1])
for line in output:
    line.small_intestine = int(line.small_intestine)/small_intestine * 1000000
    line.adrenal_gland = int(line.adrenal_gland)/adrenal_gland * 1000000
    line.HepG2 = int(line.HepG2)/HepG2 * 1000000
    line.stomach = int(line.stomach)/stomach * 1000000
    line.heart = int(line.heart)/heart * 1000000
    line.thymus = int(line.thymus)/thymus * 1000000
    line.K562 = int(line.K562)/K562 * 1000000
    line.GM12878 = int(line.GM12878)/GM12878 * 1000000
    line.H1_hESC = int(line.H1_hESC)/H1_hESC * 1000000
    line.IMR_90 = int(line.IMR_90)/IMR_90 * 1000000
    print(line)
