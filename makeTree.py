# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 19:22:08 2020
this cript will make a subset of the lactobacilli species infected by each prophage cluster
@author: Nnamdi
"""
import os
import ete3
from ete3 import Tree

lBNames = "../data/rawData/bigData.tsv"
tree = ete3.Tree("../data/rawData/lgc_protein_rooted.nw")
path = "../data/generatedData/treeSubsets/input/"


treeReps = {}
leava = []
for leaf in tree:
    x = (leaf.name).replace("'","")
    leava.append(x)
''' change cluster GCA/GCF numbers to lactobacilli names'''
with open(lBNames,'r') as text:
    for genome in text:
        number = genome.split("\t")
        genomeNum = number[0]
        LBName = number[1]
        for item in leava:
            # item = (leaf.name).replace("'","")
            if (item.strip() == genomeNum.strip()):
                treeReps[genomeNum]=LBName
                
nameCluster = []
nameReps = []
common = []
keys = []
treeCluster= {}
with open(lBNames,'r') as text:
    for genome in text:
        number = genome.split("\t")
        genomeNum = number[0]
        LBName = number[1]
        ## manually change the input file to match name of each prophage cluster
        with open("../data/generatedData/trees/inputETE3/cluster2.tsv") as file:
            for line in file:
                if (line.strip() == genomeNum.strip()):
                    treeCluster[genomeNum]=LBName
for i in treeCluster.values():
    nameCluster.append(i)
# print (nameCluster)
            
for i in treeReps.values():
    nameReps.append(i)
# # print (nameReps)

for item in nameCluster:
    for item2 in nameReps:
        if (item == item2):
            common.append(item)
for k in treeReps.keys():
    v = treeReps.get(k)
    for i in common:
        if (i==v):
            keys.append("'" + k + "'")
tree.prune(keys)
# take note in saving the output file to have same cluster name
print (tree)  
tree.write(format=1, outfile="../data/generatedData/trees/outputETE3/cluster2.nw") 


