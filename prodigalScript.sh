#! /bin/bash

# this script is used to predict genes with prodigal tool

### move to the directory where raw genome sequences are stored.
## unzip files if they are zipped
cd ../data/rawData/wholeDataSet

for doc in *.fna.gz
do 
	gunzip $doc			# unzip each file
done

for file in *.fna
do
	prodigal -i $file -a "$file.faa"		# loop through each file and analyze with prodigal
		
	mv "$file.faa" /mnt/c/Users/Nnamdi/Desktop/thesisMain/data/generatedData/genePredResults  # move to a new folder
done

