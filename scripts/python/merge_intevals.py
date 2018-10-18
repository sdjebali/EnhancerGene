#Merge the connected components
#Input : interval connected components file
#Ouput : a bedpe consensus to save in a new file

import csv
from math import exp, log
import argparse as ap
import sys
import re


class Interval(object):
    def __init__(self,chrom,start,end):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
    def __str__(self):
        return ("%s:%d-%d" % (self.chrom,self.start,self.end))

class Interaction(object):
    def __init__(self,left,right):#,indice):
        self.left = left
        self.right = right
        #self.indice = indice
    def __str__(self):
        #return("%s,%s\t%s\t%s" % (self.left,self.right,self.indice))
        return("%s,%s\t%s" % (self.left,self.right))
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

class Component(object):
    def __init__(self,list_interactions):
        self.list_interactions = list_interactions
        self.chrom_left = list_interactions[0].left.chrom
        self.chrom_right = list_interactions[0].right.chrom
        #self.indice = [i.indice for i in self.list_interactions]
        self.size = len(list_interactions)
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
    def bedpe(self):
        chrom_left = self.chrom_left
        start1 = self.left_start()
        end1 = self.left_end()
        chrom_right = self.chrom_right
        start2 = self.right_start()
        end2 = self.right_end()
        #score = moy_indices(self.indice)
        size = self.size
        return ("\t".join([chrom_left, str(start1), str(end1),chrom_right,str(start2), str(end2),str(size)]))#,str(score),str(size)]))



def getregioninfo(interval):
    p = re.compile('(?P<chrom>\w+):(?P<start>\d+)-(?P<end>\d+)')
    m = p.search(interval)
    if not m:
        print("Not a valid region %s" % interval)
        exit(1)
    chrom = m.group('chrom')
    start = int(m.group('start'))
    end = int(m.group('end'))
    return (chrom, start, end)

def new_interval(interval_str):
    [chrom,start,end] = getregioninfo(interval_str)
    if end<start:
        print("Not a valid region.")
        #exit(1)
    return (Interval(chrom,start,end))

def moy_indices (indice):
    ratio_tot = 0
    for i in range(len(indice)):
        ratio = exp(float(indice[i]))
        ratio_tot += ratio
    ratio_tot = ratio_tot/len(indice)
    return (log(ratio_tot))

def nouvelle_ligne (liste):
    chrom_left = liste[0].left.chrom
    chrom_right = liste[0].right.chrom
    test='a_faire'
    for i in liste:
        if i.left.chrom != chrom_left or i.right.chrom != chrom_right :
            test = 'erreur'
    if test != 'erreur':
        return (Component(liste))
    else:
        return ("erreur")


parser = ap.ArgumentParser()
parser.add_argument("bedpe_file_sort")
args = parser.parse_args()

with open (args.bedpe_file_sort,"rb") as fichier:
    for line in fichier:
        liste=[]
        interactions = line.split(";")
        for interaction_str in interactions:
            intervalle_str = interaction_str.split(",")
            left = new_interval(intervalle_str[1])
            right = new_interval(intervalle_str[2])
            #score = intervalle_str[2]
            liste.append(Interaction(left,right))#,score))
        component = nouvelle_ligne(liste)
        if component != "erreur":
            print(component.bedpe())
