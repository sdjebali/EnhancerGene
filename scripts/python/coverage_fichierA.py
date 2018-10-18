# This script take as input a text file from connected component with this format:
# tool;chr:beg-end,chr:beg-end;tool;chr:beg-end,chr:beg-end
# Use: python coverage_fichierA.py compconnexes_seuil12000_intervalle
# Option: --bedpe True (default: False) to have the consensus of the connected component.
import sys
import re
import argparse as ap
from collections import defaultdict

class Interval(object):
    def __init__(self,chrom,start,end):
        self.chrom = chrom
        self.start = start
        self.end = end
    def __str__(self):
        return ("%s:%d-%d" % (self.chrom,self.start,self.end))
    def lenght(self):
        return (self.end - self.start)

class Interaction(object):
    def __init__(self,left,right):
        self.left = left
        self.right = right
    def __str__(self):
        return("%s,%s,%s" % (self.tool,self.left,self.right))
    @property #not a property...
    def chrom(self):
        if self.left.chrom == self.right.chrom:
            return (True)
        else:
            return (False)
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

class Insitu(Interaction):
    def __init__(self,tool,left,right):
        Interaction.__init__(self,left,right)
        self.tool = tool

class Capture(Interaction):
    def __init__(self,tool,left,right):
        Interaction.__init__(self,left,right)
        self.tool = tool

class Chiapet_emo(Interaction):
    def __init__(self,tool,left,right):
        Interaction.__init__(self,left,right)
        self.tool = tool

class Chiapet_emq(Interaction):
    def __init__(self,tool,left,right):
        Interaction.__init__(self,left,right)
        self.tool = tool

class Component(object):
    def __init__(self,list_interactions,Id):
        self.list_interactions = list_interactions
        self.chrom = list_interactions[0].left.chrom
        self.size = len(list_interactions)
        self.Id = Id
    def left_start(self):
        left_start = min([i.left_start for i in self.list_interactions])
        return (left_start)
    def left_end(self):
        left_end = max([i.left_end for i in self.list_interactions])
        return (left_end)
    def right_start(self):
        right_start = min([i.right_start for i in self.list_interactions])
        return (right_start)
    def right_end(self):
        right_end = max([i.right_end for i in self.list_interactions])
        return (right_end)
    def bed(self):
        chrom = self.chrom
        start = self.left_start()
        end = self.right_end()
        Id = self.Id
        nb_interactions = int(self.size)
        tool = defaultdict(int)
        for i in self.list_interactions:
            tool[i.tool] += 1
        tool_str="\t".join([str(tool["cHic"]),str(tool["insitu"]),str(tool["chiapetEMO"]),str(tool["chiapetEMQ"])])
        return ("\t".join([chrom, str(start), str(end),str(Id),str(nb_interactions),tool_str]))
    def bedpe(self):
        chrom = self.chrom
        start1 = self.left_start()
        end1 = self.left_end()
        start2 = self.right_start()
        end2 = self.right_end()
        Id = self.Id
        nb_interactions = int(self.size)
        tool = defaultdict(int) #Si la cle n'existe pas, il la cree automatiquement en lui mettant une valeur de 0.
        for i in self.list_interactions:
            tool[i.tool] += 1
        tool_str="\t".join([str(tool["cHic"]),str(tool["insitu"]),str(tool["chiapetEMO"]),str(tool["chiapetEMQ"])])
        return ("\t".join([chrom, str(start1), str(end1),chrom,str(start2), str(end2),str(Id),str(nb_interactions),tool_str]))




def getregioninfo(interval):
    p = re.compile('(?P<chrom>\w+):(?P<start>\d+)-(?P<end>\d+)')
    m = p.search(interval)
    if not m:
        print("Not a valid region %s" % region)
        exit(1)
    chrom = m.group('chrom')
    start = int(m.group('start'))
    end = int(m.group('end'))
    return (chrom, start, end)

def new_interval(interval_str):
    [chrom,start,end] = getregioninfo(interval_str)
    if end<start:
        print("Not a valid region.")
        exit(1)
    return (Interval(chrom,start,end))

def new_interaction(interaction_str):
    [tool,left_str,right_str] = interaction_str.split(",")
    left = new_interval(left_str)
    right = new_interval(right_str)
    if tool == "insitu":
        return (Insitu("insitu",left,right))
    elif tool == "cHic":
        return (Capture("cHic",left,right))
    elif tool == "chiapetEMO":
        return (Chiapet_emo("chiapetEMO",left,right))
    elif tool == "chiapetEMQ":
        return (Chiapet_emq("chiapetEMQ",left,right))

def tests(line): #Is it a good input ?
    ok = 1
    p = re.compile('(?P<chrom>\w+):(?P<start>\d+)-(?P<end>\d+)*')
    m = p.search(line)
    if not m:
        ok = 0
    return (ok)



parser = ap.ArgumentParser()
parser.add_argument("fichier_compconnexes")
parser.add_argument("--bedpe",default=False)
args = parser.parse_args()

all_components = []
Id = 0
with open (args.fichier_compconnexes,"r") as f:
    for line in f:
        ok = tests(line)
        if ok == 1:
            Id += 1
            same = False
            one_component = []
            components = line.split(";")
            for component in components:
                interaction = new_interaction(component)
                if interaction.chrom:
                    one_component.append(interaction)
                    same = True
            if same:
                all_components.append(Component(one_component,Id))
    if ok == 0:
        print("Not a valid file.")

for component in all_components:
    if args.bedpe :
        print(component.bedpe())
    else:
        print (component.bed())
