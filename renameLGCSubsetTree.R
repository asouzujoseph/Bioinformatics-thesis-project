## this script will rename the tree tips of phylogenetic tree of lactobacili species infected by a cluster of prophages
library(tidyverse)
library(ape)
library(RColorBrewer)
library(ggtree)
library(treeio)
library(Ecfun)
library(phytools)
library(ggplot2)

rm(list=ls())  # clear all variables in global environment before next analysis



rm(list=ls()) 
#### subsets ##
newTree <- read.tree("C:/Users/Nnamdi/Desktop/thesisMain/data/generatedData/trees/outputETE3/cluster19.nw")
tree = midpoint.root(newTree)
treeTips = as.vector(newTree$tip.label)  #convert tree tips to a vector for easy manipulation


lab = read.delim("C:/Users/Nnamdi/Desktop/thesisMain/data/rawData/updatedNames.tsv",stringsAsFactor=FALSE) #read data
colnames(lab) <- c('GC-Number', 'Lactobacilli host')  # create headers for the file
gcNum = as.vector(lab$`GC-Number`)  # make a vector of prophage gene names
LBspecies = as.vector(lab$`Lactobacilli host`)  # make a vector of lactobacilli host name

#
# match the two data frames
same = match(treeTips,gcNum)  # match the two vectors
collect = LBspecies[same]  # use the index of matched elements in the vectors to get the corresponding lactobacilli host name

# make data frame for ggtree plot
df <- data.frame(label = newTree$tip.label, species = collect) #create a data frame for renaming the tree tips with host species name


rn <-300
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))



#visualize tree
require(ggtree)
ggtree(tree, layout = "rectangular") %<+% df + theme_tree2() +  xlim(0,2.5) +
  geom_tiplab(aes(label=species, parse = TRUE,col = species))  + scale_fill_manual(values=col_vector,aesthetics = "fill") 

# write to file
tr2 = rename_taxa(newTree, df, label, species) #rename the tree tips and write to file
write.tree(tr2, file = "C:/Users/Nnamdi/Desktop/thesisMain/data/generatedData/trees/renamedSubsetTrees/cluster19subset.nw")


##### origunal lactobacilli phylogeny###
# read the prophage cluster tree
# change the folder and  file name, then run all codes till end. 
#tree = read.tree("C:/Users/Nnamdi/Desktop/thesisMain/data/rawData/lgc_protein_rooted.tree")
#removeNATips = c("GCA_002337935.1", "GCF_000407005.1" ,"GCF_002897515.1","GCF_001885765.1", "GCF_000407025.1",
               # "GCF_001021755.1", "GCF_001886105.1", "GCF_000178235.1","GCF_001815025.1", "GCF_001047395.1" ,
               # "GCF_001878755.1", "GCF_001885805.1" ,"GCF_900104595.1","GCF_000288775.1", "GCF_001591725.1", "GCF_001832885.1" ,
               # "GCF_900167335.1", "GCF_000498295.1","GCF_900163795.1" ,"GCF_000313915.1", "GCA_003264085.1", "GCF_000315445.1" ,
               # "GCF_000374325.1","GCF_900129085.1", "GCF_900112455.1","GCF_900169305.1", "GCF_000425865.1" ,"GCF_000429585.1",
               # "GCF_000317975.2" ,"GCF_900111355.1" ,"GCF_000745125.1" ,"GCF_900113675.1", "GCF_900115825.1",
               # "GCA_003521095.1" ,"GCF_000183205.1" ,"GCF_000518205.1" ,"GCF_900110675.1" ,"GCF_000301035.1",
               # "GCF_000160075.2", "GCF_900167405.1", "GCF_900100125.1", "GCF_002884575.1", "GCF_000421665.1",
               # "GCF_000526675.1", "GCF_900111135.1", "GCF_001975685.1" ,"GCF_900129415.1")
#newTree = drop.tip(tree,removeNATips)


