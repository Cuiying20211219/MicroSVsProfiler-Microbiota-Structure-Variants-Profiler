

### Procedure for MicroSVsProfiler Framework analysis pipeline



The exploration of structural variations in the gut microbiota is crucial for understanding their role in the context of colorectal cancer. This script encompasses a comprehensive analysis of metagenomic structural variation data, utilizing a combination of shell scripting, R, Perl, and Python and aims to unveil the intricate landscape of structural variations within the gut microbiota and explore their potential as diagnostic markers for gastrointestinal cancer. These scripts accompany the 'Characterization and Diagnostic Potential of Structural Variations  in the Gut Microbiota for Colorectal Cancer' , and is intended to allow adaptations for analogous work on other datasets.

The sections in this script corresponds to procedural sections of '*Structural variation in the gut microbiome associates with host health*' and '*Characterization of gut microbial structural variations as determinants of human bile acid metabolism*'.



#### Getiing started

The data required for this study are next generation metagenomic sequencing data and related metadata. The code in the pipeline simply makes all script-wise dependencies and settings available in the present workspace. The analysis was tested using CentOS Linux release 7.9.2009 (Core) (64*Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz，526GB memory)

The code may execute under another version of packages and softwares. However, it is essential to acknowledge that the specifics of implementation of package and software may introduce minor  variations in the results. To ensure the robust functionality of the  following software applications and packages, along with their  indirect dependencies necessary for compilation and operation – some of which may need to be installed via Bioconductor –  meticulous attention must be given to the installation and loading processes.

##### local software environment:

\### perl 5, version 26, subversion 2 (v5.26.2) 

\### Platform: x86_64-linux-thread-multi

\### R, version version 4.0.5 (2021-03-31)

\### Running under: CentOS Linux 7 (Core)

\### Compiler: gcc version 4.1.2 20080704 (Red Hat 4.1.2-54)



The script was testing using the following version of packages:

##### Linux Software Install:

\### KneadData tool (version 0.7.5)

\### seqtk tool (version 1.3-r117-dirty)

\### GEM Mapper (http://gemlibrary.sourceforge.net)

\### ICRA & SGVF（）

\### Kraken2 (version 2.1.2)

\### Braken

\### The HMP Unified Metabolic Analysis Network 3.0 (HUMAnN 3.0)

\### Bowtie 64-bit

\### DIAMOND (version 0.8.38.100)



##### R package install:

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



Zeevi, D., T. Korem, A. Godneva, N. Bar, A. Kurilshikov, M. Lotan-Pompan, A. Weinberger, J. Fu, C. Wijmenga, A. Zhernakova, et al., *Structural variation in the gut microbiome associates with host health.* Nature, 2019. ***\*568\****(7750). doi: 10.1038/s41586-019-1065-y.

Wang, D., M. Doestzada, L. Chen, S. Andreu-Sanchez, I.C.L. van den Munckhof, H.E. Augustijn, M. Koehorst, A.J. Ruiz-Moreno, V.W. Bloks, N.P. Riksen, et al., *Characterization of gut microbial structural variations as determinants of human bile acid metabolism.* Cell Host Microbe, 2021. ***\*29\****(12). doi: 10.1016/j.chom.2021.11.003

Last version 01/31/2024.



