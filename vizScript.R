### this script creates a visualization of the datasets obtained from geneClusters.py
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)

##make a function to make a list of the number of phages detected in genomes
# this function will be useful if the genomes are huge and we cannot count manually
numPhages = function (x){
  y = unique(x)        # retrieves only unique values so we know the actual number of phages (for example, highest number of phages)
  z = c()             # make an empty list and add these numbers to this list as shown in the for loop
  for(i in 1:length(y)) {
    z[i] = y[i]
  }
  print(z)
}


## exploratory analysis ---> number of intact prophages
## import data
data = read.csv("C:/Users/Nnamdi/Desktop/thesisMain/data/generatedData/exploratoryData/allLBIntactPhages.csv",row.names=1)
j = numPhages(data$total)  # call the function
# # loop across each value in the specified column and sum the specified number of phages, then create a new column in each loop
sumPhages= for(i in j){
  r = data %>% group_by(total)%>% summarise(ty = sum(total==i))
  for(i in r){
    m = max(r[2])
  }
  print(m) # this values match the indices in j (each index in j represent number of phages)
}
print(j)
# since we had 1438 intact prophages then 0 represents the non-intact prophages
labels = c(0,1,2,3,4,5,6)   # the indices in j
values = c(1022,825,430,128,41,10,4)  # obtained from sumPhages and carefully arranged to match labels i.e values of j
barplot(values,names.arg=labels, col = "sky blue",border = "black", ylim = c(0,1200), ylab = "lactobacilli genomes",xlab = "Intact prophages")


## exploratory analysis ---> prophage clusters
data2 = read.csv("C:/Users/Nnamdi/Desktop/thesisMain/data/generatedData/exploratoryData/prophageClusters.csv",row.names=1)
j = numPhages(data2$total)  # call the function
# # loop across each value in the specified column and sum the specified number of phages, then create a new column in each loop
sumPhages= for(i in j){
  r = data2 %>% group_by(total)%>% summarise(ty = sum(total==i))
  for(i in r){
    m = max(r[2])
  }
  print(m) # this values match the indices in j (each index in j represent number of phages)
}
print(j)
labels = c(1,2,3,4,5)   # the indices in j
values = c(875,427,118,16,2)  # obtained from sumPhages and carefully arranged to match labels i.e values of j
barplot(values,names.arg=labels, col = "red",border = "black", ylim = c(0,1000), ylab = "lactobacilli genomes",xlab = "Prophage clusters")


## exploratory analysis ---> prophage families
data3 = read.csv("C:/Users/Nnamdi/Desktop/thesisMain/data/generatedData/exploratoryData/taxofam.csv",row.names=1)
j = numPhages(data3$total)  # call the function
# # loop across each value in the specified column and sum the specified number of phages, then create a new column in each loop
sumPhages= for(i in j){
  r = data3 %>% group_by(total)%>% summarise(ty = sum(total==i))
  for(i in r){
    m = max(r[2])
  }
  print(m) # this values match the indices in j (each index in j represent number of phages)
}
print(j)
labels = c(1,2,3,5)   # the indices in j
values = c(120,45,5,1)  # obtained from sumPhages and carefully arranged to match labels i.e values of j
barplot(values,names.arg=labels, col = "pink",border = "black", ylim = c(0,140), ylab = "lactobacilli genomes",xlab = "Number of Prophage families")


## Lifestyle (cluster level)
data4 = read.csv("C:/Users/Nnamdi/Desktop/thesisMain/data/generatedData/exploratoryData/groups.csv",row.names=1)
## Bar plot --> Lifestyle of host bacteria associated with Phage clusters
clusterLife = data4 %>% 
  group_by(lifestyle) %>% 
  summarise(clusters17=sum(clusters17), clusters3=sum(clusters3),clusters19=sum(clusters19),clusters1=sum(clusters1),
            clusters5=sum(clusters5), clusters10=sum(clusters10), clusters13=sum(clusters13),clusters2=sum(clusters2),
            clusters14=sum(clusters14), clusters9=sum(clusters9),clusters12=sum(clusters12),clusters4=sum(clusters4),
            clusters7=sum(clusters7), clusters8=sum(clusters8), clusters11=sum(clusters11), clusters6=sum(clusters6),
            clusters15=sum(clusters15),clusters18=sum(clusters18), clusters16=sum(clusters16)) 

label = c("cluster17","cluster3","cluster19","cluster1", "cluster5", "cluster10", "cluster13","cluster2",
          "cluster14","cluster9","cluster12","cluster4",  "cluster7", "cluster8", "cluster11","cluster6", "cluster15","cluster18","cluster16")
colors = c("coral4","darkturquoise","orange","violetred2","forestgreen","burlywood","darkturquoise",
           "chartreuse2","chocolate3", "coral4","cornsilk4","forestgreen","deeppink3","mediumorchid",
           "sienna3", "violetred2","forestgreen")
lifestyle= c("free living","insect","nomadic","unknown","vertebrate")
fix(clusterLife)
values = as.matrix(clusterLife[,2:20])
par(mar=c(8, 4, 2, 2) + 0.1)
barplot(values,names.arg=label,ylab="proportion based on lifestyle ",col= colors, las=2, border = 1, cex.lab=1, cex.axis=1, font=1,col.axis="black", ylim=c(0,80))
legend("topright", lifestyle, cex = 0.8, fill = colors, title = "Lifestyle")


## Lifestyle (Family level)

famLife = data3 %>% 
  group_by(lifestyle) %>% 
  summarise(Siphoviridae=sum(Siphoviridae), DNA=sum(DNA),Herelleviridae=sum(Herelleviridae),Podoviridae=sum(Podoviridae),
            Myoviridae=sum(Myoviridae)) 

label = c("Siphoviridae","DNA","Herelleviridae","Podoviridae", "Myoviridae")
colors = c("coral4","darkturquoise","orange","violetred2","forestgreen")
lifestyle= c("free living","insect","nomadic","unknown","vertebrate")
dim(famLife)
values = as.matrix(famLife[,2:6])
par(mar=c(8, 4, 2, 2) + 0.1)
barplot(values,names.arg=label,ylab="proportion based on lifestyle ",col= colors, las=2, border = 1, cex.lab=1, cex.axis=1, font=1,col.axis="black", ylim=c(0,150))
legend("topright", lifestyle, cex = 0.8, fill = colors, title = "Lifestyle")









