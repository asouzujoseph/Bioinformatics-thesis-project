# thesisProject
These codes were used to analyze genetic data and show evidence of Co-evolution between Lactobacillus bacteria and Bacteriophages

WORK FLOW
The codes used in this project follow this order:

(a) upload.sh ---> uploads genome sequences to PHASTER and downloads the results from server.

(b) prodigalScript.sh ---> this script is used to predict genes in lactobacilli genomes with prodigal software. 

(c) extractGenes.py --->  This script will extract the amino acid sequences corresponding to the gene coordinates of prophages in each genome. 

(d) orthoscript.pbs ---> The extracted prophage genes will be analyzed with Orthofinder in HPC cluster to generate Orthogroups. 

(e) extractClusters.R ---> Perform clustering on "Orthogroups.GenesCount.csv" which is the result from Orthofinder and select cut off height using dynamicTreeCut package.

(f) geneClusters.py ---> The gene clusters will be analyzed in respect to host lifestyle and taxonomy to produce datasets for visualization.

(g) vizScript.R ---> for visualization of results. 

(h) make_heatmap.R ---> creates heatmap of prophage and prophage family distribution

(i) coreGenes.py ---> produces core genes of orthogroups in each clusteer

(j) makeTrees.py ---> makes a phylogenetic tree of the lactobacilli species infected by each cluster of prophages

(k) renameSubSet.R ---> converts the accession number of the lactobacilli genomes in each subset phylogenetic tree to the name of the corresponding lactobacillus species. 

(l) sorterITOl.py ---> useful for assigning color codes to prophages in a cluster using iTOl website.



