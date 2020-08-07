#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 10:37:56 2020

This script will extract the amino acid sequences corresponding to the prophage genes using their posiions
    The generated files will be used in Orthofinder to generate Orthogroups.

@author: Nnamdi
"""
import os
import numpy as np
import pandas as pd
import re

## this function will extract the genome accession number from the original PHASTER results
def genome_name(filename):
     # pick a file without reading it. Essentially, we are modifying only the file name.
    find = filename.split("_") 
    # extracts only the accession number i.e GCA_xxxxxx | GCF_xxxxx
    rename = str(find[0]+'_'+find[1])  
    return rename
    
''' this function extracts the coordinates and contig number of all intact phages 
found in each genome, which are stored in a dictionary with genome accession number as key'''
def phageCoords(path):
    dictionary = {}
    completeGenomes = []
    bag = []
    for filename in os.listdir(path): ## loop through all files in the directory
        collector = []
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
            completeGenomes.append(genome)
            for item in newsplit:  # for each element in newsplit
                seg = item.split()    
                try:
                    score = seg[2]
                    position = seg[4]
                except IndexError:
                    score = 'null'  # catch exceptions for files empty files (i.e without phages)
                    position = 'null'
                seqPos = position.split(',')
                if ("intact" in score):         # check for only intact prophages
                    coordinates = seqPos[-1]
                    contig = seqPos[0]
                    dictionary.setdefault(genome, []).append([coordinates,contig])
                if ("intact" in score) or ("incomplete" in score) or ("questionable" in score): 
                    collector.append(genome)
        bag.append(collector)
    bag[:] = [item for item in bag if item != []] # remove empty elements
    
    # a list of lactobacilli genomes with either intact, incomplete or questionable phages
    genomesWithPhages=[]    
    for i in bag:
        for b in i:
            genomesWithPhages.append(b)
    # alist of all lactobacilli genomes        
    allGenomes = []       
    for k in completeGenomes:
        allGenomes.append(k)
    # a list of genomes without phages 
    noPhages = []
    noPhages =  [list(set(allGenomes).difference(genomesWithPhages))]
    print (noPhages)
    print ("the number of genomes without phages is :", [len([v for v in i ])for i in noPhages])
    return dictionary 


'''this function will create unique prophage names for each of the coordinates '''
def uniqName (dictionary):
    test = []   # empty list to collect coordinates associated with each key of the dictionary
    genes = {}    # make an empty dictionary 
    for key, values in dictionary.items(): # the keys are the genome accession numbers
        for item in values:
            test.append(item)
    for counter, values in enumerate(test):
        phageID = "prophage"+str(counter)+".faa"  # create a unique identifier
        genes[phageID] = values    # make a dictionary of unique phage names and gene coordinates
    
    return genes
                    

''' this function will make a dictionary that combines genome number, contig info,
unique phage names and gene coordinates'''
def comboDict (dictionary, genes):
    mix = {}    # create an empty dictionary (for downstream analysis)
    for genome in dictionary.keys():
        values = dictionary.get(genome)   # this is a list of coordinates and contig info
        for item in values:
            for name in genes.keys():    # unique prophage name
                coordinates = genes.get(name)   # this is a list of coordinates and contig info
                if (item == coordinates):
                    mix [name] = [genome,coordinates] # make a dictionary. check coordinates    
    return mix
 	

''' this function will extract the genome accession number from name of each file 
obtained from Prodigal. '''
def genome_name2(filename):
    find = filename.split("_")  # this just picks a file without reading it. 
    rename = str(find[0]+'_'+find[1])   # extracts only the accession number i.e GCA_xxxxxx | GCF_xxxxx
    return rename


''' this function will use the coordinates linked to each unique phage name to extract the gene 
sequences for that unique prophage, then it will write the list of amino acid sequence for each 
prophage gene to a file in the specified folder. These files are used for analysis in Orthofinder '''
def extractSeq (uniqName, path):
    for filename in os.listdir(path): ### loop through all files in the directory
        if filename.endswith(".faa"):       ## if file ends with ".faa" 
            genome = filename.split('.')
            genomeNum = (genome[0]+'.'+genome[1]).strip()
            file = path + filename 
            handle = open(file,'r')     # open file for reading of data
            text = handle.read()    # read all the content at a go.
            sections = text.split('>')  # split by FASTA symbol
            for segment in sections:
                sequence = '>'+segment      # reformat by appending FASTA symbol
                coords = sequence.split('#')
                try:                        # catch exceptions for empty files
                    firstCoordinate = int(coords[1])
                    contigFile = coords[0]
                except IndexError:
                    firstCoordinate = 0
                    contigFile ='null'
                for uniqName in comboDict.keys(): # for each unique prophage name
                    value = comboDict.get(uniqName) # get the values associated with the keys
                    genome = value[0]       # genome accession number
                    coords_contig = value[1]       # list containing contig info and coordinates
                    contig = coords_contig[1]   
                    coordinates = coords_contig[0]
                    separator = coordinates.split(':')  
                    pos = separator[1]
                    position = pos.split('-')
                    startPosition = int(position[0])    # start position of amino acid sequences
                    endPosition = int(position[1])      # end position of amino acid sequences
                    ###ensure that we are working on same genome containing intact phages
                    if (genome == genomeNum): 
                        # verify contig info
                        if (contig.lower().strip() in contigFile.lower().strip()): 
                            #if coordinate is within the range
                            if (startPosition <= firstCoordinate <= endPosition):
                                # create a folder to write out the AA seq of prophage genes
                                f = open("Proteome/"+str(uniqName),'a') 
                                    # write the amino acid sequence in FASTA format
                                f.write(sequence.encode().decode('unicode-escape').strip() + '\n') 
                                f.close()
                          

### uncomment to call functions
coordinates = phageCoords("../data/phasterData/downloadedResults/")
uniqName = uniqName(coordinates)
comboDict = comboDict (coordinates, uniqName)
if not os.path.exists("Proteome"):
      os.makedirs("Proteome")
extractSeq(comboDict, "../data/generatedData/prodigalResults/")


