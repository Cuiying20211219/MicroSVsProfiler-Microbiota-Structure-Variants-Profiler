**### Step 12 - Microbial SVs functional annotation**

\### The reference genomes were downloaded from progenome 

\### (http://progenomes1.embl.de/) and annotated using the web-based 

\### genome annotation service provided by PATRIC (version 3.6.6, https://www.patricbrc.org/) .

\### This code utilizes the enrichment analysis function of MicrobiomeProfiler, based on the obtained structural variations, to identify the gene functions associated with these variations. The next step involves visualizing the enriched functions.

 \```R

\### Load required library

library(MicrobiomeProfiler)

\### Read the data from the file

df <- read.table("02.Advanced.SV.anno.02AUS.txt", sep = "\t", header = TRUE)

\### Subset data based on Type

df1 <- subset(df, Type == "VSGV")
df2 <- subset(df, Type == "DSGV")

\### Perform enrichment analysis for VSGV

rst1 <- enrichKO(df1$Knum, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2)

\### Perform enrichment analysis for DSGV

rst2 <- enrichKO(df2$Knum, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2)

\### Write results to separate output files

write.table(rst1, "02.Advanced.SV.anno.02AUS.KEGG.VSGV.txt", sep = "\t", quote = FALSE)
write.table(rst2, "02.Advanced.SV.anno.02AUS.KEGG.DSGV.txt", sep = "\t", quote = FALSE)

\```

\### Visualization

\```R

\### Load necessary libraries

library(ggplot2)
library(reshape2)
library(plyr)

\### Read data from file 'Early.sum.pathway.plot1205.tax-2.txt'

L2 <- read.table('Early.sum.pathway.plot1205.tax-2.txt', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")

\### Define color palette for the plot

col.hm <- c('#550A46', '#902044', '#CB414B', '#D16F7C', '#E56B46', '#F4A862', '#F6DB86', '#DFE899', '#A4D5A0', '#62BA9F', '#3681AD', '#5A4C97', "#9970AB", "#40004B", '#9D9B9C')

\### Read additional data from file 'Early.kegg'

df1 <- read.table("Early.kegg", sep = "\t", header = FALSE)

\### Set levels of the 'Pathway' factor based on 'df1'

L2$Pathway <- factor(L2$Pathway, levels = df1$V1)

\### Output PDF file for the plot

pdf("Early.jun_anno.pdf", height = 6, width = 13)

\### Create the bar plot using ggplot2

ggplot(L2, aes(x = as.factor(Tax), y = Num, fill = Pathway)) +
  geom_bar(stat = "identity", width = 0.8, col = 'black', position = 'fill') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(hjust = 1, angle = 90),
        panel.background = element_rect(fill = NULL, colour = 'white')) +
  scale_fill_manual(values = col.hm) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_flip()

\### Turn off the PDF device

dev.off()

\```

 

\#pheatmap

\```R

\### Load necessary libraries

library(tidyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

df <- read.table("Advanced.jun_annorst4_rst.txt", sep = "\t", header = FALSE)

\### Spread the data using tidyr's spread function

df1 <- spread(df, V5, V1)

\### Write the spread data to a new file

write.table(df1, "Advanced.jun_annorst4_rst.wide.txt", sep = "\t", quote = FALSE, row.names = FALSE)

\### Specify the path to the data file for the heatmap

heatmap_file <- "/publicgp/cuiying/03project_public6/03New_CRC_SV_analysis_220808/all_PerFile_221122/01proGenomes/rst_01_Early_advanced/02analsis_plot/03chohort-gongneng/Early.DSGV.plot1.txt"

\### Read data from the heatmap file

df_heatmap <- read.table(heatmap_file, sep = "\t", header = TRUE, row.names = 1)

\### Output PDF file for the heatmap

pdf("Early.DSGV.heatmap.pdf")

\### Create the heatmap using pheatmap

pheatmap(df_heatmap, cluster_rows = TRUE, cluster_cols = FALSE, cellwidth = 25, cellheight = 25, 
         color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrBr")))(100)), 
         display_numbers = TRUE, fontsize = 8)

\### Turn off the PDF device

dev.off()

\```

 \###KEGG enrich plot

\```R

\### Load necessary library

library(ggplot2)

\### Read data from file "Early.KEGG.plot"

df <- read.table("Early.KEGG.plot", sep = "\t", header = FALSE)

\### Output PDF file for the dot plot

pdf("Early.dot.plot.importance.pdf", height = 6, width = 7)

\### Create the dot plot using ggplot2

ggplot(df, aes(V1, V2)) +
  geom_point(aes(fill = V3, size = V4), shape = 21) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text = element_text(color = 'black', size = 10)
  ) +
  scale_fill_gradient(low = "purple", high = "yellow") +
  labs(x = NULL, y = NULL)

\### Turn off the PDF device

dev.off()

\```

 

 