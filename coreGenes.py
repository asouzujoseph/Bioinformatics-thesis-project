#-*- coding: utf-8 -*-
"""
Created on Sun Jun 21 14:57:34 2020

@author: Nnamdi
"""
import pandas as pd
import os

'''
this script counts the number of core genes in the orthogroups
we applied 10% exception rule i.e the core genes must be common in at least 
90% of the prophages. Also, the prophages with more than one core genes are discarded
'''

'''this function will match the columns of the clusters 
 to the column of the orthogroup file from Orthofinder.'''
def matchColumns (file,path):
    df = pd.read_csv(file,index_col=0) # convert the csv file to a data frame
    listOfSubsets = []  #empty list to collect the data frame subsets
    for filename in os.listdir(path): ### loop through all files in the directory
        if filename.endswith(".txt"): 
            file = path+filename
            columnNames = []
            with open(file,'r') as handle: ## open the clusters folder
                for line in handle:
                    phage = line.replace('"','')
                    colls = list(df.columns)  # list of the data frame columns names
                    for name in colls:
                        if (name.strip() == phage.strip()): 
                            columnNames.append(name)
                            subsetDF = df[columnNames]   # make a subset of the matched columns
                listOfSubsets.append(subsetDF)   ## take note of tab position 
    # print (len(listOfSubsets))
    return listOfSubsets

'''this function will sum the core genes that passed the exception rule and apply 10% exception rule
## in each row of every cluster. The sum helps us to compare the best linkage method
## with more core genes which will generate better phylogenetic trees. We will use only orthogroups 
that passed the 10 % exception rule'''            
def extractCoreGenes (listOfSubsets):  
    collector = [] 
    for subset in listOfSubsets: # for each subset in the list    
        clusterGenes = [] 
        rows = range(len(subset.index))  # range covering number of rows in each subset i.e orthogroups
        rolls = list(subset.index)
        cols = range(len(subset.columns)) # range covering number of columns in each subset i.e prophages
        colls = list(subset.columns)
        for i in rows: # for each orthogroup
            count = list() # make a list for each row of the data frame
            for x in cols: ##for every column in the row i.e for each prophage per the current row
                cell = subset.iloc[i,x] #loop through each cell in the row
                if pd.notna(cell):  # if the value of the cell is not null i.e if the prophage contains genes
                    if (',' in cell)== False: # if the cell has only one value i.e if the prophage has only one gene
                        gene = cell.split('#')
                        geneName = gene[0] # name of prophage gene
                        count.append([geneName,colls[x],rolls[i]]) # append gene name, prophage name and orthogroup name as a nested list
            if len(count)>= int(0.9*len(subset.columns)): #if total number of prophages with single copy genes per row > number of columns *0.9
                clusterGenes.append(count)
        collector.append(clusterGenes)
    counter = 1
    for item in collector:
        with open ("../data/generatedData/pythonFiles/coreGenes2/cluster"+str(counter)+".tsv",'a') as file:
            file.write ('%s\n' %item)
        counter = counter+1
    
    return collector

'''
This function will swap the prophage name with the lactobacilli genome ID and make a file containing 
the core gene name, the genome ID and orthogroup name for each orthogroup in the gene family. This file will
be used as an input for progenomics tool.'''
def coreGeneID (path):
    for filename in os.listdir(path): ### loop through all files in the directory
        if filename.endswith(".tsv"):
            file = path+filename
            handle = open(file,'r')
            content = handle.readlines()
            for iii in content:
                oooo = iii.split('],')
                for item in oooo:
                    duck = item.split(',')
                    coreGene = duck[0].replace('[','').replace("'","").strip()
                    prophage = duck[1].replace("'","").strip()
                    orthogroup = duck[2].replace(']','').replace("'","").strip()
                    with open("../data/generatedData/pythonFiles/progenomicsInputFiles/clusters2/"+filename,'a') as output:
                        output.write('%s\t%s\t%s\n' % (coreGene,prophage,orthogroup))
            
            handle.close()
    

'''
This function will make a text file containing the path to the prophages containing core genes for 
each cluster. This result will be used as an input for progenomic tool to make a multiple sequence
alignment. 
'''
def pathToProtSeq (path1,path2):
    for filename in os.listdir(path1): ### loop through all files in the directory
        if filename.endswith(".faa"):       ## if file ends with ".faa" 
            pathToFile = path1+filename     ## this is the path to protein sequence file
            pathToFile = pathToFile.replace('C:',"/mnt/c")
            prophageName = filename.split('.')
            prophage = (prophageName[0]).strip()
            for cluster in os.listdir(path2):
                if cluster.endswith(".tsv"):
                    file = path2+cluster
                    with open(file,"r") as handle:
                        for line in handle:
                            content = line.split('\t')
                            prophage2 = content[1].strip().lower().replace("'","").strip()
                            if (prophage == prophage2):
                                with open("../data/generatedData/pythonFiles/progenomicsInputFiles/path/"+cluster+".txt",'a') as output:
                                    output.write('%s\n' % pathToFile)
           
              
file = "../data/generatedData/orthofinderResults/Orthogroups/Orthogroups.csv" # file obtained from OrthoFinder results
pathWard = "../data/generatedData/clustersWard/"   # path to files obtained using ward linkage method of hierarchical clustering
matchNames = matchColumns(file, pathWard)  # change to path for average linkage 
extractor = extractCoreGenes(matchNames)
coreGenomeIDs = coreGeneID("../data/generatedData/pythonFiles/coreGenes2/") 
pathToProtSeq = pathToProtSeq ("C:/Users/Nnamdi/Desktop/thesisMain/data/generatedData/Proteome/","../data/generatedData/pythonFiles/progenomicsInputFiles/clusters/" )