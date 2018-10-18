# Change an ID on interval thanks to all the files imported.
# Work only if the initial file have all the ID in the right order.
# Careful with the files names and number and with the technique spelling.

import csv

cluster=[]
with open ("composante_connexes_seuil12000","rb") as file:
    contenu=csv.reader(file, delimiter="\t")
    for ligne in contenu:
        cluster.append(ligne)

chic=[]
with open ("capture_hic_sort_ID.bedpe","rb") as file:
    contenu=csv.reader(file, delimiter="\t")
    for ligne in contenu:
        chic.append(ligne)
insitu=[]
with open ("insitu_hic_sort_ID.bedpe","rb") as file:
    contenu=csv.reader(file, delimiter="\t")
    for ligne in contenu:
        insitu.append(ligne)
emo=[]
with open ("rad21_chiapet_EMO_sort_ID.bedpe","rb") as file:
    contenu=csv.reader(file, delimiter="\t")
    for ligne in contenu:
        emo.append(ligne)
emq=[]
with open ("rad21_chiapet_EMQ_sort_ID.bedpe","rb") as file:
    contenu=csv.reader(file, delimiter="\t")
    for ligne in contenu:
        emq.append(ligne)


for i in range (len(cluster)):
    for j in range (len(cluster[i])):
        ind=cluster[i][j]
        if ind[:4]=="cHic":
            if chic[int(ind[4:])-1][len(chic[0])-1] != ind:
                print("false, chic")
            else:
                output = ','.join(("cHic",chic[int(ind[4:])-1][0]))
                output = ':'.join((output,chic[int(ind[4:])-1][1]))
                output = '-'.join((output,chic[int(ind[4:])-1][2]))
                output = ','.join((output,chic[int(ind[4:])-1][3]))
                output = ':'.join((output,chic[int(ind[4:])-1][4]))
                output = '-'.join((output,chic[int(ind[4:])-1][5]))

                cluster[i][j]=output

        elif ind[:6]=="insitu":
            if insitu[int(ind[9:])-1][len(insitu[0])-1] != ind:
                print("false, insitu")
            else:
                output = ','.join(("insitu",insitu[int(ind[9:])-1][0]))
                output = ':'.join((output,insitu[int(ind[9:])-1][1]))
                output = '-'.join((output,insitu[int(ind[9:])-1][2]))
                output = ','.join((output,insitu[int(ind[9:])-1][3]))
                output = ':'.join((output,insitu[int(ind[9:])-1][4]))
                output = '-'.join((output,insitu[int(ind[9:])-1][5]))

                cluster[i][j]=output

        elif ind[:10]=="chiapetEMO":
            if emo[int(ind[10:])-1][len(emo[0])-1] != ind:
                print("false,emo")
            else:
                output = ','.join(("chiapetEMO",emo[int(ind[10:])-1][0]))
                output = ':'.join((output,emo[int(ind[10:])-1][1]))
                output = '-'.join((output,emo[int(ind[10:])-1][2]))
                output = ','.join((output,emo[int(ind[10:])-1][3]))
                output = ':'.join((output,emo[int(ind[10:])-1][4]))
                output = '-'.join((output,emo[int(ind[10:])-1][5]))

                cluster[i][j]=output

        elif ind[:10]=="chiapetEMQ":
            if emq[int(ind[10:])-1][len(emq[0])-1] != ind:
                print("false,emq")
            else:
                output = ','.join(("chiapetEMQ",emq[int(ind[10:])-1][0]))
                output = ':'.join((output,emq[int(ind[10:])-1][1]))
                output = '-'.join((output,emq[int(ind[10:])-1][2]))
                output = ','.join((output,emq[int(ind[10:])-1][3]))
                output = ':'.join((output,emq[int(ind[10:])-1][4]))
                output = '-'.join((output,emq[int(ind[10:])-1][5]))

                cluster[i][j]=output

for c in range (len(cluster)):
    i=0
    output = ''
    while i<len(cluster[c]):
        if output == '':
            output = cluster[c][i]
        else:
            output = ';'.join((output,cluster[c][i]))
        i+=1
    print(output)
