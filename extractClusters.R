
## this script will cluster the prophages based on thier gene counts

library(dendextend)
library(vegan)
library(factoextra)
library(ggplot2)
library(ggdendro)
library(dynamicTreeCut)

dat = read.csv("Orthogroups.GeneCount.csv", row.names=1)  #read file obtained from Orthofinder
data = t(dat)  # transpose file
fix(data) # produce a matrix of the data in a new window for visual inspection
distMatrix = vegdist(data) # use vegdist function for calculating distance matrix between the gene counts
hier = hclust(distMatrix, method = "ward.D2")  # hierarchical clustering
save(distMatrix,hier,file="matrices.rda")  # save file to disk for ease of reading 
load(file="matrices.rda")  # reload saved file
groups = cutreeHybrid(dendro = hier, distM = as.matrix(distMatrix), deepSplit = 1) #use dynamicTreeCut package to determine cut height
labels = groups$labels
max(labels)

## function to extract members of each cluster
extractMembers = function(data, labels){
  k = 19  # number of seeds or clusters
  for (i in 1:k){ # loop over each seed
    members = data.frame(data[labels == i,]) # pick a cluster
    names=row.names(members)  # extract names of cluster members
    # for each cluster, write the members to a text file
    write.table(names, file = paste("C:/Users/Nnamdi/Desktop/thesisMain/data/generatedData/clusters/","clusters",toString(i),".txt", sep=""), sep = "\t",
                row.names = FALSE, col.names = FALSE) 
  }
}
extractMembers(data,labels) #call function











