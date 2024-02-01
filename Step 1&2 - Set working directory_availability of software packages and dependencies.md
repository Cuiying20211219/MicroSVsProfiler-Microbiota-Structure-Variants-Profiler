

\### **Step 1 - Set working directory**

\###

\### The code is designed to search for input files and store output files

\### within subdirectories relative to a central directory. On your computer, 

\### establish a directory for the analysis. To ensure generality, we 

\### will refer to this directory as "top/." Subsequently, create subdirectories

\### to achieve the following directory structure:

\###

\###  top/ 

\###   code/  // containing this script and any others.

\###   data/   // containing input files.

\###   results/ // location of output files.

\###

\###

\###

**### Step 2 - Ensure availability of software packages and satisfaction of dependencies**

\###

\### The code in this section simply makes all script-wise dependencies

\### and settings available in the present workspace. The analysis was 

\### tested using CentOS Linux release 7.9.2009 (Core).

\### (64*Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz，526GB memory)

\### The code may execute under another version of packages and 

\### softwares. However, it is essential to acknowledge that the specifics

\### of implementation of package and software may introduce minor 

\### variations in the results. To ensure the robust functionality of the 

\### following software applications and packages, along with their 

\### indirect dependencies necessary for compilation and operation 

\### – some of which may need to be installed via Bioconductor – 

\### meticulous attention must be given to the installation and 

\### loading processes.

 

\### **local software environment:**

\### perl 5, version 26, subversion 2 (v5.26.2) 

\### Platform: x86_64-linux-thread-multi

\### R, version version 4.0.5 (2021-03-31)

\### Running under: CentOS Linux 7 (Core)

\### Compiler: gcc version 4.1.2 20080704 (Red Hat 4.1.2-54)

\### The script was testing using the following version of packages:

 

\### **Linux Software & Installation:**

\### KneadData tool (version 0.7.5)

\### seqtk tool (version 1.3-r117-dirty)

\### GEM Mapper (http://gemlibrary.sourceforge.net)

\### ICRA & SGVF

\### Kraken2 (version 2.1.2)

\### Braken

\### The HMP Unified Metabolic Analysis Network 3.0 (HUMAnN 3.0)

\### Bowtie 64-bit

\### DIAMOND (version 0.8.38.100)

 

\### **2.1 GEM Mapper version 2 Install4  **

\```shell

bzip2 -dc GEM-binaries.tbz2

\#added into PATH variables：

export PATH="~/SVGFinder/GEM-binaries/:$PATH"

source ~/.bashrc

\```

\### **2.2 Python environment creation**

\```shell

conda create -n python3.8.1 python=3.8.1 

\```

\### Python modules recommend:

\### numpy (tested with 1.14.2)

\### biopython (tested with 1.68)

\### ujson (tested with 1.35)

\### pandas (tested with 0.23.4)

\### scipy (tested with 1.1.0)

\### bokeh (tested with 0.12.6)

\### Cpp11

\### Cython, it can also been installed with Anaconda

 

\### **2.3 SGVFinder Download and Install **

\```shell

\# Download the tools and databases, and install the tools

wget https://zenodo.org/record/3237975/files/DataFiles.tar.gz

gzip -dc DataFiles.tar.gz > DataFiles

git clone https://github.com/segalab/SGVFinder.git

mv DataFiles SGVFinder/

python setup.py build

\```

\### **2.4 kneaddata Install ** 

\```shell

creat envirmont and install kneaddata

conda create -n kneaddata

conda activate kneaddata

conda install -c biobakery kneaddata

\# download database

mkdir kneaddata_database

cd kneaddata_database/

kneaddata_database --download human_genome bowtie2 ./

#kneaddata_database --download mouse_C57BL bowtie2 ./

kneaddata_database --download human_transcriptome bowtie2 ./

kneaddata_database --download ribosomal_RNA bowtie2 ./

\```

\### **2.5 Seqtk Install **

\```shell

git clone https://github.com/lh3/seqtk.git

cd seqtk

make

\```

\### **2.6 HUMAnN 3.0 Install**

\```shell

\# Build a virtual environment and install python3.7

conda create --name mpa -c bioconda python=3.7

\# install humann by conda

source activate mpa

conda install humann

\# run the following command to test whether the installation was successful

humann_test

\```

 

\### **2.7 Kraken2  (version 2.1.2) install**

\```shell

conda install -y kraken2

\```

\### **2.8 DIAMOND install**

\```shell

wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz

tar xzf diamond-linux64.tar.gz

\# creating a diamond-formatted database file

./diamond makedb --in reference.fasta -d reference

\```

\### **2.9 DIAMOND install**

\```shell

 conda install graphlan

\```



\### **2.10 R package install**:

\```R

\### Fast read of large files into R

library (data.table)

\### Clustering software

library (flashClust) 

\### Partial Spearman correlations, for confounder analysis. Previously reported work done using v1.0

library (ppcor) 

\### Plotting

library (gplots)

\### Plotting; to arrange several plots on the same page

library (cowplot)

\### Plotting

library (ggplot2) 

\### Data transformations

library (plyr) 

\```

\### locale:

\### [1] da_DK.UTF-8/da_DK.UTF-8/da_DK.UTF-8/C/da_DK.UTF-8/da_DK.UTF-8

\### 

\### attached base packages:

\### [1] stats   graphics grDevices utils   datasets methods  base   

\### 

\### other attached packages:  

\### 

\### loaded via a namespace (and not attached):

\### [1] "UpSetR"    "plyr"     "rfPermute"  "rfUtilities" "AUCRF"

\### [6] "pROC"     "caret"    "randomForest" "FactoMineR"  "forcats"

\### [11] "stringr"   "dplyr"    "purrr"    "readr"    "tidyr"

\### [16] "tibble"    "tidyverse"  "wesanderson" "cluster"   "fpc"

\### [21] "vegan"    "lattice"   "permute"   "meta"     "gggenes"

\### [26] "ggpubr"    "reshape2"   "pheatmap"   "ggsci"    "RColorBrewer"

\### [31] "ggplot2"   "stats"    "graphics"   "grDevices"  "utils"

\### [36] "datasets"   "methods"   "base"   “fpc”  ”MicrobiomeProfiler“  ”clusterProfiler“

 \### [37] "org.Hs.eg.db"  

\###

\###

\###

 

 

