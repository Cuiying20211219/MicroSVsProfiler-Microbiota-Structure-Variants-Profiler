**### Step 13 - Correlation analysis between SV features of significant differences with metabolites**

\### **01 differential Metablites**

\### Read and prepare data 

df<-read.table("df.csv",sep=",",header=T,row.names=1)

\### Initialize vectors

Pvalue12<-c(rep(1,nrow(df)))
Pvalue13<-c(rep(1,nrow(df)))
Pvalue23<-c(rep(1,nrow(df)))

HC_mean<-c(rep(1,nrow(df)))
Early_mean<-c(rep(1,nrow(df)))
Advanced_mean<-c(rep(1,nrow(df)))

\### Compute Wilcoxon rank sum test and save results

for(i in 3:ncol(df)){
    y12=try(wilcox.test(as.numeric(subset(df,Group=="Healthy")[,i]),as.numeric(subset(df,Group=="Early")[,i]),paired=FALSE))
     y13=try(wilcox.test(as.numeric(subset(df,Group=="Healthy")[,i]),as.numeric(subset(df,Group=="Advanced")[,i]),paired=FALSE))
     y23=try(wilcox.test(as.numeric(subset(df,Group=="Early")[,i]),as.numeric(subset(df,Group=="Advanced")[,i]),paired=FALSE))

try(Pvalue12[i]<-y12$p.value)
try(Pvalue13[i]<-y13$p.value)
try(Pvalue23[i]<-y23$p.value)

HC_mean[i]<-as.numeric(mean(as.numeric(subset(df,Group=="Healthy")[,i])))
Early_mean[i]<-as.numeric(mean(as.numeric(subset(df,Group=="Early")[,i])))
Advanced_mean[i]<-as.numeric(mean(as.numeric(subset(df,Group=="Advanced")[,i])))

}

\### Output results

out<-cbind(colnames(df),Pvalue12,Pvalue13,Pvalue23,HC_mean,Early_mean,Advanced_mean) 
write.table(out,file="Met.wilcox.txt",quote=FALSE,sep="\t",row.names=FALSE)



##### #### Metablites Complex heatmap plot

\### Load necessary library

library("ComplexHeatmap")

df <- read.table("Complex.Heatmap.new72.txt", sep = "\t", header = TRUE, row.names = 1)

\### Scale the data

df_scaled <- t(scale(t(df)))

\### Read additional data files

dfp <- read.table("Complex.Heatmap.pvalue.new72.txt", sep = "\t", header = TRUE, row.names = 1)
dfm <- read.table("Complex.Heatmap.Metabolism2.txt", sep = "\t", header = TRUE)
dfo <- read.table("Complex.Heatmap.Organismal2.txt", sep = "\t", header = TRUE)
dfh <- read.table("Complex.Heatmap.Human2.txt", sep = "\t", header = TRUE)
dfe <- read.table("Complex.Heatmap.Enver2.txt", sep = "\t", header = TRUE)
dfg <- read.table("Complex.Heatmap.Genetic2.txt", sep = "\t", header = TRUE)

\### Define color palette for annotations

col_heat <- list(
  Metabolism = c("Amino acid metabolism" = "#e31a1c", "Biosynthesis of other secondary metabolites" = "#a6cee3", "Carbohydrate metabolism" = "#33a02c", "Lipid metabolism" = "#fb9a99", "Metabolism of cofactors and vitamins" = "#ff7f00", "Metabolism of other amino acids" = "#cab2d6", "Nucleotide metabolism" = "#ffff99", "Other" = "#bdbdbd", "Other Metabolism" = "#737373"),
  Organismal = c("Digestive system" = "#8dd3c7", "Endocrine system" = "#ffffb3", "Environmental adaptation" = "#bebada", "Sensory system" = "#fb8072", "Other" = "#bdbdbd"),
  Human.Diseases = c("Cancer: overview" = "#80b1d3", "Drug resistance: antineoplastic" = "#fdb462", "Endocrine and metabolic disease" = "#b3de69", "Other" = "#bdbdbd"),
  Environmental.Information.Processing = c("Membrane transport" = "#fccde5", "Signal transduction" = "#d9d9d9", "Other" = "#bdbdbd"),
  Genetic.Information.Processing = c("Folding, sorting and degradation" = "#bc80bd", "Translation" = "#ccebc5", "Other" = "#bdbdbd")
)

\### Output PDF file for the heatmap

pdf("Figure6A.plot.pdf", height = 17, width = 11)

\### Create the heatmap using ComplexHeatmap

Heatmap(df_scaled,
        cluster_columns = FALSE,
        circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), c('#2c7bb6', '#abd9e9', '#FFFFFF', '#fdae61', '#fdae61')),
        column_names_gp = grid::gpar(fontsize = 20),
        row_names_gp = grid::gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Spearman's rho", direction = 'horizontal', at = c(-1, -0.5, 0, 0.5, 1)),
        cell_fun = function(j, i, x, y, w, h, col) {grid.text(dfp[i, j], x, y)},
        row_km = 5,
        left_annotation = rowAnnotation(
          Metabolism = dfm[, "Metabolism"],
          Organismal.Systems = dfo[, "Organismal.Systems"],
          Human.Diseases = dfh[, "Human.Diseases"],
          Environmental.Information.Processing = dfe[, "Environmental.Information.Processing"],
          Genetic.Information.Processing = dfg[, "Genetic.Information.Processing"],
          col = col_heat
        )
)

\### Turn off the PDF device

dev.off()



##### #### circlize &  metabloites-SVs heatmap Visualization

\# Visualization

\```R

\### Load necessary libraries

library(circlize)
library(viridis)
library(reshape2)
library(tidyverse)
library(wesanderson)

\### Read data for VSGV plot

df_vsgv <- read.table("VSGV.diff.plot1.test3.anno.txt", sep = "\t", header = TRUE)

\### Read order files for VSGV

vsgv_bacter_order <- read.table("vsgv.bacter.order", sep = "\t", header = FALSE)
vsgv_metab_order <- read.table("vsgv.metab.order", sep = "\t", header = FALSE)

met_vsgv <- vsgv_metab_order$V1
bac_vsgv <- vsgv_bacter_order$V1

\###Output PDF file for VSGV circos plot

pdf("VSGV.cluster_ba.circos.pdf", width = 14, height = 14)

\### Create chord diagram for VSGV

chordDiagram(df_vsgv, annotationTrack = "grid", grid.col = c(wes_palette("Darjeeling1", 14, type = "continuous"), rep('grey', 19)),
             order = order1,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df_vsgv))))))

\### Add sector names

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")

  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), col = "black")
}, bg.border = NA)

dev.off()

\### Read data for DSGV plot

df_dsgv <- read.table("DSGV.test5.anno.txt", sep = "\t", header = TRUE)

\### Read order files for DSGV

dsgv_bacter_order <- read.table("dsgv.bacter.order", sep = "\t", header = FALSE)
dsgv_metab_order <- read.table("dsgv.metab.order", sep = "\t", header = FALSE)

\###Output PDF file for DSGV circos plot

pdf("DSGV.cluster_ba.circos.pdf", width = 14, height = 14)

\### Create chord diagram for DSGV

chordDiagram(df_dsgv, annotationTrack = "grid", grid.col = c(wes_palette("Darjeeling1", 14, type = "continuous"), rep('grey', 38)),
             order = c(dsgv_metab_order$V1, dsgv_bacter_order$V1),
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df_dsgv))))))

\### Add sector names

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")

  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), col = "black")
}, bg.border = NA)

dev.off()

\### Read data for correlation plot for VSGV

df1 <- read.table("heatmap.rst.DSGV", sep = "\t", header = TRUE, row.names = 1)
f1 <- read.table("Dsgv.cy", sep = "\t", header = FALSE)
df_vsgv_cor <- df1[f1$V1, ]

\### Output PDF file for VSGV correlation plot

pdf("Vsgv.cor.pdf")

\### Create correlation heatmap for VSGV

pheatmap(as.matrix(df_vsgv_cor), col = colorRampPalette(c("#02CBFB", "#FFFFFF", "#F20C34"))(15),
         display_numbers = TRUE, number_format = "%.2f", cellwidth = 30, cellheight = 30, fontsize_number = 8)

dev.off()

\### Read data for correlation plot for DSGV

df2 <- df1[f$V1, ]
f2 <- read.table("DSGV.cy", sep = "\t", header = FALSE)
df_dsgv_cor <- df2[f2$V1, ]

\###Output PDF file for DSGV correlation plot

pdf("Dsgv.cor.pdf", width = 20, height = 30)

\### Create correlation heatmap for DSGV

pheatmap(as.matrix(df_dsgv_cor), col = colorRampPalette(c("#02CBFB", "#FFFFFF", "#F20C34"))(15),
         display_numbers = TRUE, number_format = "%.2f", cellwidth = 30, cellheight = 30, fontsize_number = 8)

dev.off()

\```

 

 

 

 

 

 