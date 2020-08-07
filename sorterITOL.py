# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 22:56:29 2020

@author: Nnamdi
"""
''' this script will sort the prophageNames to match the tip labels as shown in iTol.
The output is useful for assigning color ranges to the tree labels in iTol'''
dicti = {}
with open ("C:/Users/Nnamdi/Desktop/filer.txt",'r') as file1:
    for line in file1:
        prophage = line.strip().lower()
        with open("../data/generatedData/trees/phageHosts/cluster19.tsv.csv",'r') as file2:
            for line in file2:
                pick = line.split(',')
                prophageName = pick[0].strip().lower()
                LBName = pick[1]
                if (prophage == prophageName):
                    dicti[prophage] = LBName
with open("../data/generatedData/trees/phageHosts/cluster19Sorted.csv","a") as f:     
    for key in  dicti.keys():
        f.write("%s,%s"%(key, dicti[key])) 
      
                