#Give the connexions for "connected_components.py".

import argparse as ap


class Interval(object):
    def __init__(self,chrom,start,end):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
    def __str__(self):
        return ("%s:%d-%d" % (self.chrom,self.start,self.end))

class Interaction(object):
    def __init__(self,left,right,name,indice):
        self.left = left
        self.right = right
        self.name = name
        self.indice = indice
    def __str__(self):
        return("%s,%s,%s" % (self.left,self.right,self.indice))
    @property
    def left_start(self):
        return (self.left.start)
    @property
    def left_end(self):
        return (self.left.end)
    @property
    def right_start(self):
        return (self.right.start)
    @property
    def right_end(self):
        return (self.right.end)


def condition(line1,line2):
    if str(line1.left.chrom)==str(line2.right.chrom):
        if str(line1.left) == str(line2.left):
            if str(line1.right_start) == str(int(line2.right_end)+1):
                return("arrete")
            elif str(int(line1.right_end)+1) == str(line2.right_start):
                return("arrete")
            else:
                return("presque_arrete")
        elif str(line1.right) == str(line2.right):
            if str(line1.left_start) == str(int(line2.left_end)+1):
                return("arrete")
            elif str(int(line1.left_end)+1) == str(line2.left_start):
                return("arrete")
            else:
                return("presque_arrete")
        else:
                if str(line1.left_start) == str(int(line2.left_end)+1):
                    return("presque_arrete")
                elif str(int(line1.left_end)+1) == str(line2.left_start):
                    return("presque_arrete")
                else:
                    return("non_arrete")
    else:
        return("stop")



parser = ap.ArgumentParser()
parser.add_argument("fichier_bedpe_sort")
args = parser.parse_args()
#Bedpe as input: cHic_sort.bedpe

liste=[]
with open (args.fichier_bedpe_sort,"rb") as fichier:
    for line in fichier:
        interaction = line.split("\t")
        if interaction[0] != "#chr":
            left = Interval(interaction[0],interaction[1],interaction[2])
            right = Interval(interaction[3],interaction[4],interaction[5])
            if left.chrom==right.chrom: #If same chromosome
                liste.append(Interaction(left,right,interaction[6],interaction[7]))

arretes = []
for line1 in range(len(liste)):
    line2 = line1+1
    stop = 0
    while stop == 0:
        if line2 == len(liste):
            stop = 1
        else:
            if condition(liste[line1],liste[line2]) == "arrete":
                arretes.append([str(liste[line1]),str(liste[line2])])
            elif condition(liste[line1],liste[line2]) == "non_arrete":
                stop = 1
            elif condition(liste[line1],liste[line2]) == "stop":
                stop = 1
            line2 += 1


for c in arretes:
    print ("\t".join(c))
