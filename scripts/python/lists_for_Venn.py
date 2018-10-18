# Use: python listes_Venn.py compconnexes.txt
# Careful with the files names and techniques spelling.

import csv
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument("compconnexes")
args = parser.parse_args()

insitu=[]
cap=[]
chiapet=[]
with open (args.compconnexes,"rb") as compco:
	for line in compco:
		groupes = line.split("\t")
		for i in groupes:
			if i[:4]=="chic":
				cap.append(groupes[0])
			elif i[:6]=="insitu":
				insitu.append(groupes[0])
			elif i[:7]=="chiapet":
				chiapet.append(groupes[0])

with open ("liste_insitu.txt","wb") as insitu_ecriture:
	ecriture=csv.writer(insitu_ecriture)
	for i in insitu:
		ecriture.writerow([i])

with open ("liste_cHic.txt","wb") as cap_ecriture:
	ecriture=csv.writer(cap_ecriture)
	for i in cap:
		ecriture.writerow([i])

with open ("liste_chiapet.txt","wb") as chiapet_ecriture:
	ecriture=csv.writer(chiapet_ecriture)
	for i in chiapet:
		ecriture.writerow([i])
