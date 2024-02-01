**### Step 3 - Import input files**

\###

\### The following input data is assumed to exist within the data subdirectory:

\###

\### **- raw data-files: **

\### Fastq files containing the raw metagenomic sequencing data for each

\### individual sample(step 4 - step 5). 

\###

\### **-metadata-file: **

\### File with clinical phenotypes (columns) per individual (rows). Used to 

\### test for associations with, or for confounder analysis. Individuals are labeled

\### SRR/DRR number. (step 4 - step 13)

\###

\### **-genome Annotation-files**

\### GTF file containing genome annotations 

\### of over 3000 representative species. (step 5)

\###

\### **-metabolomic-file**: 

\### Input data matrix for abundance of 325 metabolites per individual.

\###  Note that no additional normalization is done in this script, so 

\### data is assumed to be comparable in these regards. (step 13)

\###  

\### **-SV profile file: **

\### The provided CSV file comprises normalized coverage depth data 

\### pertaining to structural variants within DSGV (Diverse Structural Genome Variants) 

\### and VSGV (Variant Structural Genome Variants), obtained as output from SGVFinder, 

\### a structural variant detection tool. In this dataset, individual rows correspond 

\### to distinct structural variation types, while individual columns represent individual

\### subjects or samples.

 

\### All demonstration files are tab-delimited text files, but other formats would 

\### equally work after modifying the respective file-import commands in the script.

\###

\###

\```R

options (stringsAsFactors = FALSE)

\### **- metadata.tab**

phenotypes = read.table (file = "data/metadata.txt", header = T, sep = "\t")

ctrl = rownames (subset (phenotypes, Disease_Stage == "HC")) ### set of control samples for analyses restricted to health individuals.

crc = rownames (subset (phenotypes, Disease_Stage == "CRC")) ### set of control samples for analyses restricted to CRC individuals.

 

\### ** - metabolomic.tab**

metabolomic = read.table (file = "data/metabolomic.tab", row.names = 1, header = T, sep = "\t")

 

\### - taxonomy.tab

mgs_taxonomy = read.table (file = "data/MGS_taxonomy.csv", sep = "\t", row.names = 1, header = T)

 

\### -DSGV_VSGV.CSV

 df_dSV= read.table ("data/ALL.dsgv.csv", sep = ",", header=T,row.names=1)

df_vSV= read.table ("data/ALL.vsgv.csv", sep = ",", header=T,row.names=1)

\```