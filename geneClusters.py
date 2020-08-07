#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 10:37:56 2020

This script generates datasets which are visualized as results and also as the phylogenetic trees

@author: Nnamdi
"""
import os
import re
import numpy as np
import pandas as pd



'''this function will extract the genome accession number from the original PHASTER results'''
def genome_name(filename):
    find = filename.split("_")  # this just picks a file without reading it. Essentially, we are modifying only the file name.
    rename = str(find[0]+'_'+find[1])   # extracts only the accession number i.e GCA_xxxxxx | GCF_xxxxx
    return rename
    
''' this function extracts the coordinates and contig number of all intact phages found in each genome, 
which are stored in a dictionary with genome accession number as key'''
def phageCoords(path):
    dictionary = {}
    for filename in os.listdir(path): ### loop through all files in the directory
        if filename.endswith("_genomic"):       ## if file ends with ".fna" 
            genome = genome_name(filename).strip()  ## extract accession number from filename
            file = path + filename          
            input = open (file,'r') # the raw file downloaded from PHASTER
            text=input.read()   # read all the file content
            #Spliting data at \n
            split=text.split('\\n')
            ## useful data starts from index 32 till the second last element. The last element is '}'.
            useful=split[32:]
            newsplit=useful[:-1] # this is because the last element is "}" and it is not useful
            for item in newsplit:  # for each element in newsplit
                seg = item.split()    
                try:
                    score = seg[2]
                    position = seg[4]
                except IndexError:
                    score = 'null'        # catch exceptions
                    position = 'null'
                seqPos = position.split(',')
                if ("intact" in score):         # check for only intact prophages
                    coordinates = seqPos[-1]
                    contig = seqPos[0]
                    dictionary.setdefault(genome, []).append([coordinates,contig])
    # print (len(dictionary))
    return dictionary 


'''this function will create unique prophage names for each of the coordinates and contigs'''
def uniqName (dictionary):
    test = []   # empty list to collect coordinates associated with each key of the dictionary
    genes = {}    # make an empty dictionary 
    for key, values in dictionary.items(): # the keys are the genome accession numbers
        for item in values:
            test.append(item)
    for counter, values in enumerate(test):
        phageID = "prophage"+str(counter)+".faa"  # create a unique identifier
        genes[phageID] = values    # make a dictionary of unique phage names and gene coordinates
    #print (len(genes))
    return genes

def exploratoryAnalysis (dictionary, genes):
    explore = {} # create an empty dictionary to collect data for exploratory anaylsis
    for genome in dictionary.keys():
        values = dictionary.get(genome)   # this is a list of coordinates and contig info
        for item in values:
            for name in genes.keys():    # unique prophage name
                coordinates = genes.get(name)   # this is a list of coordinates and contig info
                if (item == coordinates):
                    explore[name] = genome
    ### generate dataset for exploratory analysis
    ### this will show the number of intact prophage genes in each lactobacilli genome
    uniq = []
    for k,v in explore.items():         
        #for item in v:     
        uniq.append(v) 
    uniqList=list(set(uniq)) 
    matrix=np.zeros((len(uniqList),len(explore.keys())))
    ### fill the matrix
    col=0  # starting at the first column ( which will hold the names of ratings )
    for v in explore.values(): #    for each phage(key) in the dictionary
        x = 1
        index=uniqList.index(v)  # assign the position of that rating in the list to a variable ---> index
        matrix[index][col]= np.cumsum(x)  # assign 1 to a cell of the matrix if the rating present at a specified row 
        x=x+1
        col=col+1 
        
    ### convert to pandas data frame and write to excel file
    dRow = [i for i in uniqList]  #    Rows are populated by elements of the list ---> uniqScore
    dCol = [j for j in list(explore)]     # columns are populated by keys of the dictionary ----> specieScore
    df = pd.DataFrame(matrix)     
    df.index = dRow     # since the matrix is transposed, then the columns becomes the rows
    df.columns = dCol
    df.to_excel(excel_writer="../data/generatedData/exploratoryData/allLBIntactPhages.xlsx")                
    return explore

'''this function will make a dictionary that combines genome number, contig info,
unique phage names and gene coordinates'''
def comboDict (dictionary, genes):
    mix = {}    # create an empty dictionary
    for genome in dictionary.keys():
        values = dictionary.get(genome)   # this is a list of coordinates and contig info
        for item in values:
            for name in genes.keys():    # unique prophage name
                coordinates = genes.get(name)   # this is a list of coordinates and contig info
                if (item == coordinates):  #list of coordinates and contig info
                    mix [name] = [genome,coordinates] # make a dictionary. check coordinates
    return mix
	
                                
'''this function will extract the accession number of the most common phage in each contig'''
def phageAccessionNum (path):
    phageAccNum = {}
    for filename in os.listdir(path): ### loop through all files in the directory
        if filename.endswith("_genomic"):       ## if file ends with ".fna" 
            genome = genome_name(filename).strip()  ## extract accession number from filename
            file = path + filename          
            input = open (file,'r') # the files downloaded from PHASTER
            text=input.read()   # read all the file content
            #Spliting data at \n
            split=text.split('\\n')
            ## useful data starts from index 32 till the second last element. The last element is '}'.
            useful=split[32:]
            newsplit=useful[:-1] # this is because the last element is "}" and it is not useful
            for item in newsplit:  # for each element in newsplit
                seg = item.split()    
                try:
                    mostCommon=[]
                    score = seg[2]
                    phage = seg[13]
                    mostCommon.append(phage.split(',')[0]) #pick only the most common phage i.e the first Phage
                except IndexError:
                    score = 'null'        # catch exceptions
                    phage = 'null'
                if ("intact" in score):         # check for only intact prophages
                    for phage in mostCommon:
                        x=(re.sub('^(.*)(?=NC+.)',"",phage))    #this will generate list of accession number of extracted prophages . ##useful for taxonomy classification  
                        accessionNum=(re.sub(r'\([^()]*\)', '',x))  # this code will remove the parenthesis attached to the accession number. 
                        phageAccNum.setdefault(genome,[]).append(accessionNum)
    
    return phageAccNum


'''this function will make a dictionary of prophage name and cluster ID  ''' 
def phageClusters (path):
    ### extract the prophage names found in each prophage cluster / species
    clusters = {}
    for filename in os.listdir(path):
        if filename.endswith (".txt"):
            name = filename.split('.')
            clusterID = name[0]
            file = path + filename
            with open(file,'r') as handle:
                for line in handle:
                    uniqProphageName = line.replace('"','').lower().strip()
    # make a dictionary using prophage name as key and cluster ID of the prophage as value
                    clusters[uniqProphageName]=clusterID 
    
    return clusters



''' this function will replace the prophage names in each lactobacilli genome with
the name of the prophage cluster to which it belongs. Hence, we can know the number of
prophage clusters in each lactobacilli genome'''
def phageGenomeID (mix,clusters): #taxonomyFile,phageAccNum,lactobacilliReps):   
    genomePhageClusters = {}
    for key,value in mix.items():
        phage = key.split('.')
        uniqPhageName = phage[0]
        genomeID = value[0]
        for prophageName,clusterID in clusters.items():
            if (prophageName == uniqPhageName):
                genomePhageClusters.setdefault(genomeID,[]).append(clusterID) ### write to a csv file for tree ordered heatmap
    uniq = []
    for k,v in genomePhageClusters.items():         
        for item in v:     
            uniq.append(item) 
    uniqList=list(set(uniq)) 
    matrix=np.zeros((len(uniqList),len(genomePhageClusters.keys())))
   
    # ### fill the matrix
    col=0  # starting at the first column of the matrix
    for k,v in genomePhageClusters.items(): #    for each phage(key) in the dictionary
        for item in v:
            index=uniqList.index(item)  # use the position of the item as the index
            matrix[index][col]= 1 # assign 1 to a cell of the matrix if the item is present  
        col=col+1 # move the next column
        
    ### convert to pandas data frame and write to excel file
    dRow = [i for i in uniqList]  #    Rows are populated by elements of the list 
    dCol = [j for j in list(genomePhageClusters)]     # columns are populated by keys of the dictionary 
    df = pd.DataFrame(matrix).T   # transpose the matrix  
    df.index = dCol     # since the matrix is transposed, then the columns becomes the rows
    df.columns = dRow
    df.to_excel(excel_writer="../data/generatedData/exploratoryData/prophageClusters.xlsx")
    return genomePhageClusters

'''this function will change the prophage cluster name to corresponding prophage family name'''
def phageFamily (taxonomyFile,genomePhageClusters, phageAccDict):
    phageFamily = {}
    with open (taxonomyFile) as taxo:
        for line in taxo:
            separator = line.split('\t')
            phageAccNum = separator[0].strip()
            phageInfo = separator[1].split()
            phageFam = phageInfo[-3]
            # for genome, clusters in genomePhageClusters.items():
            for LBgenome,accessionNum in phageAccDict.items():
                for phageNCNum in accessionNum:
                    # if (genome == LBgenome):
                    if (phageAccNum == phageNCNum.strip()):
                        phageFamily.setdefault(LBgenome,[]).append(phageFam)
    
    return phageFamily

''' this function will change host bacteria genome accession number to the corresponding 
Lactobacilli species name'''
def speciesName (lactobacilliReps,genomePhageClusters,phageFamily ):
    genomeNumToName = {}
    genomeNumToName2 ={}
    with open(lactobacilliReps,'r') as text:
        for genome in text:
            number = genome.split("\t")
            genomeIdentityNum = number[0]    #accession number
            LBName = number[1]   # lactobacillus species name
            
            for key3,value3 in phageFamily.items(): ## use for family level analysis
                if (key3.strip().lower() == genomeIdentityNum.strip().lower()):
                    genomeNumToName2.setdefault(LBName,[]).append(value3)
            
            for genomeID, value in genomePhageClusters.items(): # use for cluster level analysis
                if (genomeID.strip().lower() == genomeIdentityNum.strip().lower()):
                    genomeNumToName.setdefault(LBName,[]).append(value)
                    
    
    ## generate a dataset for visualization 
    ## alternate between both dictionaries and rename the output file accordingly
    uniq = []
    for k,v in genomeNumToName.items():
        for items in v: 
            for x in items:
                uniq.append(x) 
    uniqList=list(set(uniq)) 
    matrix=np.zeros((len(uniqList),len(genomeNumToName.keys())))
   
    ### fill the matrix
    col=0  # starting at the first column ( which will hold the names of ratings )
    for k,v in genomeNumToName.items(): #    for each phage(key) in the dictionary
        #x = 1
        for item in v:
            for cluster in item:
                index=uniqList.index(cluster)  # assign the position of that rating in the list to a variable ---> index
                matrix[index][col]= 1 #np.cumsum(x)    # assign 1 to a cell of the matrix if the rating present at a specified row 
            #x=x+1
        col=col+1 
        
    ### convert to pandas data frame and write to excel file
    dRow = [i for i in uniqList]  #    Rows are populated by elements of the list ---> uniqScore
    dCol = [j for j in list(genomeNumToName)]     # columns are populated by keys of the dictionary ----> specieScore
    df = pd.DataFrame(matrix).T     
    df.index = dCol     # since the matrix is transposed, then the columns becomes the rows
    df.columns = dRow
    df.to_excel(excel_writer="../data/generatedData/exploratoryData/clus.xlsx")
    return genomeNumToName
    
    
''' this function will give a lifestyle to each LB genome based on the ref file'''
def lifestyle(speciesName,lifeStyleFile): 
    clustersLifeStyle = {}
    with open (lifeStyleFile) as file:
        for line in file:
            text = line.split(",")
            species = text[0].lower().strip()
            life = text[1]
            
            for key,value in speciesName.items():
                LBName = key.lower().strip()
                if (species == LBName):
                    clustersLifeStyle[LBName] = life
    ggg = []
    with open("../data/generatedData/exploratoryData/clusterLife.csv", 'a') as f:     
        for key in  clustersLifeStyle.keys():
            ggg.append(key)
            f.write("%s,%s"%(key, clustersLifeStyle[key]))
    hhh = list(set(ggg))
    print (len(hhh))
    return  clustersLifeStyle


#### objective 2 : PHYLOGENETIC TREES

'''This function will get the name of lactobacilli species infected by a prophage'''
def phageHosts (mix, lactobacilliReps):
    treeFile = {}
    with open(lactobacilliReps,'r') as text:
        for genome in text:
            number = genome.split("\t")
            genomeIdentityNum = number[0]
            LBName = number[1]
            for key,value in mix.items():
                phage = key.split('.')
                uniqPhageName = phage[0]
                genomeID = value[0]
                if (genomeID == genomeIdentityNum):
                    treeFile[uniqPhageName] = LBName
    return treeFile


'''This function will make a file for each cluster and match each prophage to its host '''
def phageHosts2 (path, treeFile):
    for filename in os.listdir(path):
        phageHosts = {}
        if filename.endswith (".tsv"): ## loop the through the core genes in each cluster
            file = path + filename
            with open(file,'r') as handle:
                for line in handle:
                    content = line.split('\t')
                    prophage = content[1].lower().strip() # extract the prophage name 
                    for k,v in treeFile.items():
                        prophage2 = k.lower().strip()
                        if (prophage == prophage2):
                            phageHosts[prophage]=v
        # print (len(phageHosts))
        with open("../data/generatedData/trees/phageHosts/"+filename+".csv","a") as f:     
            for key in  phageHosts.keys():
                f.write("%s,%s"%(key, phageHosts[key])) 
                            
'''this function will retrieve the accession number of all LB genomes infected by phages in 
a cluster for subsetting with ete3 package'''
def treeSubSet (explore, path):
    for filename in os.listdir(path):
        if filename.endswith (".tsv"): ## loop the through the core genes in each cluster
            file = path + filename
            with open(file,'r') as handle:
                for line in handle:
                    content = line.split('\t')
                    prophage = content[1].lower().strip() # extract the prophage name
                    for key in explore.keys(): # loop through the prophages in this dictionary
                        value = explore.get(key)  ## accession number of infected host
                        phage = key.split('.')
                        prophage2 = phage[0].strip().lower()
                        if (prophage == prophage2): #if both prophages match
                            with open("../data/generatedData/trees/inputETE3/"+filename,"a") as output:
                                output.write('%s\n' %value) #write the accession number


### files
clustersFolder = "../data/generatedData/clustersWard/"
lgcReps = "../data/rawData/updatedNames.tsv"
speciesLifeStyle = "../data/rawData/latestLifeStyles.txt"
taxonomyFile = "../data/rawData/taxoGenomes.txt"
coregenes = "../data/generatedData/pythonFiles/progenomicsInputFiles/clusterCoreGenes/"

### uncomment to call functions
coordinates = phageCoords("../data/phasterData/downloadedResults/")
uniqName = uniqName(coordinates)
explore = exploratoryAnalysis(coordinates,uniqName)
comboDict = comboDict (coordinates, uniqName)
phageAccNum = phageAccessionNum("../data/phasterData/downloadedResults/")
clusters = phageClusters(clustersFolder)
phageGenomeID = phageGenomeID(comboDict,clusters)
phageFamily = phageFamily(taxonomyFile,phageGenomeID, phageAccNum)
speciesName = speciesName(lgcReps, phageGenomeID,phageFamily )                      
lifestyle = lifestyle(speciesName, speciesLifeStyle)  
tree = phageHosts(comboDict,lgcReps)
tree2 = phageHosts2 (coregenes,tree )
subset = treeSubSet (explore,coregenes )



