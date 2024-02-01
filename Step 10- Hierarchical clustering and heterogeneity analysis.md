**### Step 10-** **Hierarchical clustering and heterogeneity analysis**

##### 10.1 ANOVA-Type Analysis

\### To evaluate the extent of heterogeneity in microbial traits across the 

\### four distinct cohorts, we conducted an analysis similar to ANOVA 

\### (Analysis of Variance). Within this analysis, we decomposed 

\### the total variance observed in the abundance or structural variants 

\### (SVs) associated with a specific microbial species. We compared 

\### this variance to the portion explained by disease status and the 

\### portion explained by the cohort factor, similar to a linear model.

\### This framework included both the disease status parameter 

\### and the cohort factor as explanatory variables, aiming to understand

\### their contributions to the abundance of species or SVs characteristics.



\````R

\### Load required libraries

library(tidyverse)

\### Load data

df <- read.table("ALL.vsgv.dsgv.name_1.T.txt", sep = "\t", header = TRUE, row.names = 1)
df[is.na(df)] <- 0
df <- t(df)
metadata <- read.table("all.data.metadata.4cohort.0201.txt", sep = "\t", header = TRUE)

\### Subset data based on column sums

df1 <- df[, colSums(df != 0) > 200]
df2 <- df1[metadata$ID,]

\### Calculate sample similarity indices for disease and study

ss.disease <- apply(t(df2), 1, FUN = function(x, label) {
  rank.x <- rank(x) / length(x)
  ss.tot <- sum((rank.x - mean(rank.x))^2) / length(rank.x)
  ss.o.i <- sum(vapply(unique(label), function(l) {
    sum((rank.x[label == l] - mean(rank.x[label == l]))^2)
  }, FUN.VALUE = double(1))) / length(rank.x)
  return(1 - ss.o.i / ss.tot)
}, label = metadata %>% pull(Group))

ss.study <- apply(t(df2), 1, FUN = function(x, label) {
  rank.x <- rank(x) / length(x)
  ss.tot <- sum((rank.x - mean(rank.x))^2) / length(rank.x)
  ss.o.i <- sum(vapply(unique(label), function(l) {
    sum((rank.x[label == l] - mean(rank.x[label == l]))^2)
  }, FUN.VALUE = double(1))) / length(rank.x)
  return(1 - ss.o.i / ss.tot)
}, label = metadata %>% pull(Study))

\### Create and write results to a file

plot1 <- cbind(ss.study, ss.disease, num = colSums(df2 != 0))
write.table(plot1, "SV.feature.txt", sep = "\t", quote = FALSE)

\### Define custom theme for ggplot

cy_theme <- theme(
  panel.grid = element_blank(),
  panel.background = element_rect(color = 'black', fill = 'transparent', size = 1),
  axis.text.x = element_text(size = 16, color = "black"),
  axis.text.y = element_text(size = 14, color = "black"),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  plot.margin = margin(t = 3, r = 3, b = 3, l = 3, unit = "cm"),
  plot.title = element_text(size = 16, face = "bold"),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16)
)

\### Create ggplot

ggplot(as.data.frame(plot1), aes(x = ss.disease, y = ss.study), fill = "gray", col = "black") +
  geom_point(aes(size = num, fill = num), shape = 21, col = alpha(c('black'), alpha = 0.4)) +
  scale_fill_viridis_c() +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0.01, slope = 1) +
  geom_abline(intercept = -0.01, slope = 1) +cy_theme
  

##### 10.2 Hierarchical Clustering Analysis for Structural Variations:

\### Clustering Analysis

\### Within our study, we employed the pheatmap package to conduct 

\### hierarchical clustering analysis.

\### The primary goal of this analysis was to assess dissimilarities 

\### between early colorectal cancer patients and healthy individuals,

\### as well as between advanced colorectal cancer patients and healthy 

\### individuals within each cohort.

\### By following these analytical steps, we aimed to comprehensively 

\### evaluate the heterogeneity of microbial traits across cohorts and 

\### assess structural variations within specific subgroups related to colorectal cancer.

 

\```R

\### Load required library

library(pheatmap)

\### Read data and annotations

df <- read.table("dsgv.vsgv.pheatmap.not.all.txt", sep = "\t", header = TRUE, row.names = 1)
Canno <- read.table("dsgv.vsgv.pheatmap.not.all.col.anno.txt", sep = "\t", header = TRUE, row.names = 1)
Ranno <- read.table("dsgv.vsgv.pheatmap.not.all.row.anno.txt", sep = "\t", header = TRUE, row.names = 1)

\### Create a PDF file for the heatmap

pdf("dsgv.vsgv.heat.not.all.0315.pdf", height = 18, width = 65)

\###  Generate heatmap with clustering

pheatmap(df, cluster_cols = TRUE, annotation_row = Canno, annotation_col = Ranno[c('Type')],
         cutree_row = 7, show_colnames = FALSE, border_color = "white")

\### Close the PDF file

dev.off()

\### Generate hierarchical clustering and save it to a file

heatmapobj <- pheatmap(df)
geneClust <- heatmapobj$tree_row

\### Perform clustering

geneCluster2Three <- cutree(geneClust, k = 7)

\### Create a new PDF for clustering results

pdf("dsgv.vsgv.a.not.all7.pdf", height = 100, width = 100)

\###  Print clustering results

print(geneCluster2Three)

\### Plot clustering with rectangular boxes

plclust(geneClust)
rect.hclust(geneClust, k = 7)

\### Close the PDF file

dev.off()

\```

 