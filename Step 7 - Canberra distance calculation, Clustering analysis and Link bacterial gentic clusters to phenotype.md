**### Step 7 - Canberra distance calculation, Clustering analysis and Link bacterial gentic clusters to phenotype**

 

\### This script built in R and Perl language involves merging 

\### vSV and sSV data, calculating the Canberra distance, 

\### standardizing distance matrices, using computational tools

\### to assess genetic relatedness among microbial communities, 

\### and visualizing research results. The calculation and visualization 

\### were carried out using the R package 'vegan', and 'ade4', and 

\### statistical analysis of covariate effects using the 'adonis' function 

\### from the 'vegan' package. 





\### preprocessing data of SVs

\### Canberra distance calculation and Visualization

\### Canberra Distance Metric Calculation

\### The vSV (variant structural variants) and dSV (deletion structural 

\### variants) profiles of 51 species were combined.

\### Subsequently, the Canberra distance was computed across all samples.

\### This distance metric quantifies dissimilarity between samples. Specifically,

\### distance matrices for all samples were calculated using the vegdist() function

\### from the R package vegan (version 2.5-6). This function facilitates the 

\### computation of dissimilarity metrics.

 

\```R

shared_sv_dis<-function(inDf){

 \#inDf<-all_sv

 inList <- split(as.matrix(inDf), seq(nrow(inDf)))

 

 shared_n_func<-function(inVec,Vec_i){

  \#inVec<-inList[[3]]

  Vec_i<-Vec_i

  insvdf<-data.frame(inVec,Vec_i) %>% na.omit

  sv_dis<- vegdist(t(insvdf), method = "canberra")

  \#length(na.omit(Vec_i+inVec))

  return(sv_dis)

 }

 

 marker_n_mat<-matrix(NA, nrow = nrow(inDf), ncol = nrow(inDf))

 for (i in 1:nrow(inDf)) { #nrow(inDf)

  \#i<-2

  Vec_i<-inList[[i]]

  shared_n_i<-sapply(inList, shared_n_func,Vec_i = Vec_i)

  marker_n_mat[i,]<-shared_n_i

 }

 rownames(marker_n_mat)<-rownames(inDf)

 colnames(marker_n_mat)<-rownames(inDf)

 

 return(marker_n_mat)

}

\```



\### Principal Coordinates Analysis (PCoA) & Data Visualization

\### The assessment was performed at the SV level using Principal 

\### Coordinates Analysis (PCoA). This multivariate technique helps 

\### visualize the dissimilarity patterns in microbial communities. 

\### To visualize the results, PCOAs were plotted using the 'ade4' package,

\### enhancing the interpretability of dissimilarity patterns among microbial communities.

\### This code creates a plot for visualizing differentiated SVs between different Group and Study.

\### Assessing Covariate Effects

\### The pipeline then involved evaluating the effect size and statistical 

\### significance of various covariates on the structure of microbial communities.



\### adonis Function

\### This assessment was conducted using the 'adonis' function within 

\### the 'vegan' package. The 'adonis' function employs a permutational 

\### multivariate analysis of variance (PERMANOVA) approach to 

\### assess the influence of covariates on community structure.

\```R

p_all_sv_pcoa_DSGV<-ggplot(all_sv_avg_dist_pcoa_DSGV,aes(X1,X2, color = map_sv_rmna$Group,shape=map_sv_rmna$Study))+

 geom_point(size = 2,alpha = 0.5)+

 stat_ellipse(aes(group = map_sv_rmna$Group, fill = map_sv_rmna$Group, color = map_sv_rmna$Group) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend

= F)+

 xlab(paste("PCo1=",round(all_sv_avg_dist_mds_DSGV$eig[1],digits = 2),"%",sep = ""))+

 ylab(paste("PCo2=",round(all_sv_avg_dist_mds_DSGV$eig[2],digits = 2),"%",sep = ""))+

 scale_color_manual(values = c("#b2182b","#4393c3"))+

 theme(plot.subtitle = element_text(vjust = 1),

​    plot.caption = element_text(vjust = 1),

​    axis.line.x = element_line(),

​    axis.line.y = element_line(),

​    legend.position = 'bottom',

​    legend.title = element_blank(),

​    legend.key = element_rect(fill = NA),

​    panel.grid.major = element_line(colour = NA),

​    panel.grid.minor = element_line(colour = NA),

​    panel.background = element_rect(fill = NA))

\# boxplot

p_all_sv_boxplot_DSGV<-ggExtra::ggMarginal(p_all_sv_pcoa_DSGV, type = "boxplot", groupColour = TRUE, groupFill = TRUE,xparams = list(color = 'white'), yparams = list(color = 'white'))

source("/publicgp/cuiying/03project_public6/02_CRC_SV_analysis/CRC_SV-2021_12_19_script/Groningen-Microbiome-master/Projects/SV_BA/functions.R")

df<-read.table("ALL.vsgv.dsgv.HC_CRC.txt",sep="\t",header=T,row.names=1)

metadata<-read.table("CRC_HC_CRA.pheno.all.dicovery.HC_CRC.txt",sep="\t",header=T,row.names=1)

metadata$ID<-as.character(metadata$ID)

df<-df[metadata$ID,]

sv<-read.table("DE",sep="\t",header=F)

df<-df[,sv$V1]

df[is.na(df)]<-0

dsgv<-df

map<-metadata

\# Create data frame containing distances between each sample after removing empty rows

all_shared_sv_dis_DSGV<-shared_sv_dis(dsgv)

all_msv_dist_std_avg_DSGV<-all_shared_sv_dis_DSGV

all_msv_dist_std_avg_rmna_DSGV<-dist_rmna(all_msv_dist_std_avg_DSGV)

map_sv_rmna<- map[match(rownames(all_msv_dist_std_avg_rmna_DSGV), rownames(map)),]

\# Create multidimensional scaling to compare before and after

all_sv_avg_dist_mds_DSGV<-cmdscale(all_msv_dist_std_avg_rmna_DSGV, k=5, eig = T)

\## Create a primary PCoA plot using the PCoA data frame

all_sv_avg_dist_pcoa_DSGV <- data.frame(all_sv_avg_dist_mds_DSGV$points)

\#

p_all_sv_pcoa_DSGV<-ggplot(all_sv_avg_dist_pcoa_DSGV,aes(X1,X2, color = map_sv_rmna$Group,shape=map_sv_rmna$Study))+

 geom_point(size = 2,alpha = 0.5)+

 stat_ellipse(aes(group = map_sv_rmna$Group, fill = map_sv_rmna$Group, color = map_sv_rmna$Group) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend

= F)+

 xlab(paste("PCo1=",round(all_sv_avg_dist_mds_DSGV$eig[1],digits = 2),"%",sep = ""))+

 ylab(paste("PCo2=",round(all_sv_avg_dist_mds_DSGV$eig[2],digits = 2),"%",sep = ""))+

 scale_color_manual(values = c("#b2182b","#4393c3"))+

 theme(plot.subtitle = element_text(vjust = 1),

​    plot.caption = element_text(vjust = 1),

​    axis.line.x = element_line(),

​    axis.line.y = element_line(),

​    legend.position = 'bottom',

​    legend.title = element_blank(),

​    legend.key = element_rect(fill = NA),

​    panel.grid.major = element_line(colour = NA),

​    panel.grid.minor = element_line(colour = NA),

​    panel.background = element_rect(fill = NA))

\# Create a boxplot to overlay on the above plot 

p_all_sv_boxplot_DSGV<-ggExtra::ggMarginal(p_all_sv_pcoa_DSGV, type = "boxplot", groupColour = TRUE, groupFill = TRUE,xparams = list(color = 'white'), yparams = list(color = 'white'))

\# adonis test

adonis(as.dist(all_shared_sv_dis_DSGV)~.,data = as.data.frame(map[,"Group"])) 

\```



\### SVs Clustering Technique

\### To group the samples based on their genetic dissimilarity, 

\### K-medoids clustering was employed. This technique partitions a dataset

\### into K clusters, with each cluster represented by a medoid, which is 

\### the object closest to the median of the cluster by 

\### Partitioning Around Medoids (PAM) was utilized. Samples were assigned 

\### to distinct clusters, with a specified cluster number 'k'. The determination of

\### the best cluster number 'k' was carried out using the estimated optimum 

\### average silhouette width. This evaluation was done using the 'pamk'

\### function within the 'fpc' package.

 

\### Chi-Square Tests

\### Once the samples were clustered, Chi-Square tests were performed. 

\### These tests aimed to explore potential associations and dependencies 

\### between the profiles of structural variants (SVs) of species and the . 

\### The Chi-Square tests are associated with a host phenotype of interest. 

\### In the example work we describe here, this was disease status.

\### Specifically, CRC stages were categorized into two groups: 

\### Early stage (TNM stage I-II) and Advanced stage (TNM stage III-IV).

 

\### This code predominantly consolidates the structural variation features 

\### of individual bacteria into separate files.

\```perl

#!/usr/bin/perl -w
use strict;

#####Process VSGV.bacteria file

open F, "VSGV.bacteria";
while (<F>) {
    chomp;
    my @l = split(/\t/, $_);

open O, ">$l[0].txt";

#####Open ALL.vsgv.dsgv.name_1.T.txt file

open IN, "ALL.vsgv.dsgv.name_1.T.txt";
while (<IN>) {
    chomp;
    my @l1 = split(/\t/, $_);
    my @l2 = split(/:/, $l1[0]);

if ($. == 1) {
    print O "$_\n";
    next;
}

#####Check for matching IDs

if ($l[0] eq $l2[0]) {
    print O "$_\n";
}

}

close(IN);
close(O);

}

close(F);

#!/usr/bin/perl -w
use strict;

#####Process DSGV.bacteria file

open F, "DSGV.bacteria";
while (<F>) {
    chomp;
    my @l = split(/\t/, $_);

open O, ">$l[0].txt";

#####Open ALL.vsgv.dsgv.name_1.T.txt file

open IN, "ALL.vsgv.dsgv.name_1.T.txt";
while (<IN>) {
    chomp;
    my @l1 = split(/\t/, $_);
    my @l2 = split(/:/, $l1[0]);

if ($. == 1) {
    print O "$_\n";
    next;
}

#####Check for matching IDs

if ($l[0] eq $l2[0]) {
    print O "$_\n";
}

}

close(IN);
close(O);

}

close(F);

#!/usr/bin/perl -w
use strict;

my $n;

#####Process DSGV_1121098.PRJDB570.dist.txt file and generate DSGV_1121098.PRJDB570.dist.rst1.txt

open O, ">DSGV_1121098.PRJDB570.dist.rst1.txt";
for (my $i = 1; $i <= 689; $i++) {
    open IN, "DSGV_1121098.PRJDB570.dist.txt";
    while (<IN>) {
        chomp;
        my @l = split(/\t/, $_);

if ($. == 1) {
    $n = $l[$i - 1];
    next;
}

print O "$l[0]\t$n\t$l[$i]\tDSGV_1121098.PRJDB570\n";

}

close(IN);

}

close(IN);
close(O);

#!/usr/bin/perl -w
use strict;

my (%h1, %h2, %h3);

#####Process NT.bray_dis.all.txt file

open IN, "NT.bray_dis.all.txt";
while (<IN>) {
    chomp;
    my @l = split(/\t/, $_);
    my $key1 = "$l[0]\t$l[1]";
    my $key2 = "$l[1]\t$l[0]";

#####Check for unique pairs

if ((not exists $h2{$key1}) && (not exists $h2{$key2})) {
    print "$_\n";
    $h1{$key1} = 0;
    $h2{$key2} = 0;
}

}

close(IN);

#!/usr/bin/perl -w
use strict;

my (%h1, %h2, %h3);

#####Process all.data.metadata.4cohort.0201.txt file

open F, "all.data.metadata.4cohort.0201.txt";
while (<F>) {
    chomp;
    my @l = split(/\t/, $_);
    $h1{$l[0]} = $l[1];
    $h2{$l[0]} = $l[5];
    $h3{$l[0]} = $l[6];
}

close(F);

#####Process Dist.SV.all.rst file and print formatted output

open IN, "Dist.SV.all.rst";
while (<IN>) {
    chomp;
    my @l = split(/\t/, $_);

if ($. == 1) {
    print "S1\tS2\tCb\tBa\tG1\tG2\tSt\n";
}

if ((exists $h1{$l[0]}) && (exists $h1{$l[1]}) && (exists $h2{$l[0]}) && (exists $h2{$l[1]}) && ($l[2] ne "NA") && ($l[0] ne $l[1])) {
    print "$_\t$h1{$l[0]}-$h1{$l[1]}\t$h2{$l[0]}-$h2{$l[1]}\t$h3{$l[0]}-$h3{$l[1]}\n";
}

}

close(IN);

\```

 

\### This code primarily conducts centroid clustering on the structural variation 

\### features of single bacteria, considering all their structural variation 

\### characteristics.

\```R

#####Load required libraries

library(fpc)
library(cluster)
library(graphics)
library(wesanderson)
library(tidyverse)
library(FactoMineR)
library(factoextra)

#####Load distance matrix from file "Dist.SV.435591.txt"

df <- read.table("Dist.SV.435591.txt", sep="\t", header=TRUE, row.names=1)

#####Remove rows and columns with missing values

df1 <- df[rowSums(is.na(df)) == 0, rowSums(is.na(df)) == 0]

#####Determine the optimal number of clusters using PAM (Partitioning Around Medoids) method

fit <- pamk(df1)
num_clusters <- fit$nc

#####Perform PAM clustering

otu_pam <- pam(df1, num_clusters)

#####Read metadata from "all.data.metadata.4cohort.0201.txt"

meta <- read.table("all.data.metadata.4cohort.0201.txt", sep="\t", header=TRUE, row.names=1)

#####Extract metadata for samples in each cluster

meta1 <- meta[fit$pamobject$clustering, ]

#####Display the contingency table of Type1 and cluster assignments

table(as.data.frame(cbind(meta1['Type1'], fit$pamobject$clustering)))

#####Perform a chi-squared test on the contingency table

chisq_result <- chisq.test(table(as.data.frame(cbind(meta1['Type1'], fit$pamobject$clustering))))

#####Print the chi-squared test results

print(chisq_result)

#####Iterate through filenames in "filename2" file

file <- read.table("filename2", sep="\t", header=FALSE)

for (i in file$V1) {
    print(i)

#####Read data from current filename

​    df <- read.table(i, sep="\t", header=TRUE, row.names=1)

#####Remove rows and columns with missing values

df1 <- df[rowSums(is.na(df)) == 0, rowSums(is.na(df)) == 0]

#####Determine optimal number of clusters using PAM

fit <- pamk(df1)
num_clusters <- fit$nc

#####Perform PAM clustering

otu_pam <- pam(df1, num_clusters)

#####Extract metadata for samples in each cluster

meta1 <- meta[row.names(as.data.frame(fit$pamobject$clustering)), ]

#####Perform chi-squared test on Type2 and cluster assignments

chisq_result <- chisq.test(table(as.data.frame(cbind(meta1['Type2'], fit$pamobject$clustering))))

#####Print the contingency table and chi-squared test results

print(table(as.data.frame(cbind(meta1['Type2'], fit$pamobject$clustering))))
print(chisq_result)

#####Write results to a file "all.rst.Type2"

cat(i, fit$n, chisq_result$p.value, as.numeric(chisq_result$statistic), "\n", sep="\t", file="all.rst.Type2", append=TRUE)

}

#####Visualization using factoextra and wesanderson palette

for (i in file$V1) {
    print(i)
    df <- read.table(i, sep="\t", header=TRUE, row.names=1)
    df1 <- df[rowSums(is.na(df)) == 0, rowSums(is.na(df)) == 0]
    fit <- pamk(df1)
    otu_pam <- pam(df1, fit$nc)

#####Visualization and saving to PDF

pdf(file=paste0("plot_", i, ".pdf"))
fviz_cluster(otu_pam, ellipse.type="norm", geom="point", ggtheme=theme_minimal(),
             palette="Set2", pointsize=4, ellipse.level=0.87, main="1121115")
dev.off()

}

\```









 