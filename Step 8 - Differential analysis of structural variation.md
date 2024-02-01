**### Step 8 - Differential analysis of structural variation**

\### In this step, this pipeline encompasses the calculation of a dissimilarity 

\### matrix, K-medoids clustering using the PAM algorithm, determination 

\### of the optimal cluster number, and Chi-Square tests to investigate 

\### associations between bacterial SVs clusters and CRC progression 

\### stages. This comprehensive approach allows for the exploration of 

\### potential relationships in the dataset.

 

\### Through the difference between vSV and dSV data forms, this 

\### process performs a difference analysis on the two types of data 

\### in different groups. Our grouping method is mainly based on disease 

status similar to step7.

\## library packages

 

\```R

library(tidyverse)

library(ggpubr)

library(GGally)

library(ggmosaic)

library(cowplot)

library(ggsci)

library(ggthemes)

library(gplots)

library(VennDiagram)

library(circlize)

library(viridis)

library(wesanderson)

library(ggalluvial)

library(reshape2)

library(Hmisc)

library(vegan)

library(microbiome)

library(mediation)

 

##### 8.1 Convert SV name

\```R

changeSVname<-function(SVrawid){

 testname   <- SVrawid

 species_name <- as.character(taxonomy$organism[match(str_replace_all(testname, "\\..*",""), taxonomy$X)])

 region    <- str_replace_all(testname, ".*\\:","")

 region_list <- str_split(region,";")

 

 region_suf  <- NULL

 i<-1

 for (i in c(1:length(region_list))){

  if(length(region_list[[i]]) <= 2){

   region_suf[i] <- paste(":",region[i],sep = "")

  }else{

   region_suf[i] <- paste(":",region_list[[i]]," and ",length(region_list[[i]])-1," segments", sep = "")

  }

  i <- i+1

 }

 paste(species_name,region_suf,sep = "")

}

\```

 

##### 8.2 VSGV Differential Analysis

\### The code calculates the Wilcoxon rank sum test for vSV between two 

\### groups: HC (Healthy Controls) and Early (early stages of a condition) for VSGV.

\### HC-Early vSV Differential Analysis

\### Load and prepare data 

vsgv <- read.table("ALL.vsgv.name_1.csv", sep = ",", header = TRUE, row.names = 1)
meta <- read.table("all.data.metadata1211_Early_Advanced0131.txt", sep = "\t", header = TRUE)

\### Subset data

meta_all <- subset(meta, Type2 != "CRA" & Type2 != "Advanced")
meta_all$ID <- as.character(meta_all$ID)
vsgv <- vsgv[as.character(meta_all$ID),]
vsgv <- t(vsgv)
meta_all$ID == colnames(vsgv)

\###  Initialize vectors

Pvalue <- rep(1, nrow(vsgv))
HC_mean <- rep(1, nrow(vsgv))
Early_mean <- rep(1, nrow(vsgv))

\###  Compute Wilcoxon rank sum test and save results

for (i in 1:nrow(vsgv)) {
  y <- try(wilcox.test(as.numeric(vsgv[i, subset(meta_all, Type2 == "HC")[,"ID"]]),
                       as.numeric(vsgv[i, subset(meta_all, Type2 == "Early")[,"ID"]]), paired = FALSE))
  try(Pvalue[i] <- y$p.value)
  HC_mean[i] <- as.numeric(summary(as.numeric(vsgv[i, subset(meta_all, Type2 == "HC")[,"ID"]]))[4])
  Early_mean[i] <- as.numeric(summary(as.numeric(vsgv[i, subset(meta_all, Type2 == "Early")[,"ID"]]))[4])
}

\### Output results

out <- cbind(rownames(vsgv), Pvalue, HC_mean, Early_mean, rep("all", nrow(vsgv)))
write.table(out, file = "vsgv.wilcox.HC_Early.all0131.txt", quote = FALSE, sep = "\t", row.names = FALSE)

\# HC-Advanced vSV Differential Analysis

vsgv<-read.table("ALL.vsgv.name_1.csv",sep=",",header=T,row.names=1)

meta<-read.table("all.data.metadata1211_Early_Advanced0131.txt",sep="\t",header=T)

meta_all<-subset(meta,Type2!="CRA" & Type2!="Early")

meta_all$ID<-as.character(meta_all$ID)

vsgv<-vsgv[as.character(meta_all$ID),]

vsgv<-t(vsgv)

meta_all$ID==colnames(vsgv)

 

Pvalue<-c(rep(1,nrow(vsgv)))

HC_mean<-c(rep(1,nrow(vsgv)))

Advanced_mean<-c(rep(1,nrow(vsgv)))

for(i in 1:nrow(vsgv)){

 y=try(wilcox.test(as.numeric(vsgv[i,subset(meta_all,Type2=="HC")[,"ID"]]),as.numeric(vsgv[i,subset(meta_all,Type2=="Advanced")[,"ID"]]),paired=FALSE))

 try(Pvalue[i]<-y$p.value)

 HC_mean[i]<-as.numeric(summary(as.numeric(vsgv[i,subset(meta_all,Type2=="HC")[,"ID"]]))[4])

 Advanced_mean[i]<-as.numeric(summary(as.numeric(vsgv[i,subset(meta_all,Type2=="Advanced")[,"ID"]]))[4])

}

out<-cbind(rownames(vsgv),Pvalue,HC_mean,Advanced_mean,rep("all",nrow(vsgv)))

write.table(out,file="vsgv.wilcox.HC_Advanced.all.txt",quote=FALSE,sep="\t",row.names=FALSE)

\```

 

\### HC-Advanced vSV Differential Analysis

\### Read and prepare data 

vsgv <- read.table("ALL.vsgv.name_1.csv", sep = ",", header = TRUE, row.names = 1)
meta <- read.table("all.data.metadata1211_Early_Advanced0131.txt", sep = "\t", header = TRUE)

\### Subset data

meta_all <- subset(meta, Type2 != "CRA" & Type2 != "Early")
meta_all$ID <- as.character(meta_all$ID)
vsgv <- vsgv[as.character(meta_all$ID),]
vsgv <- t(vsgv)
meta_all$ID == colnames(vsgv)

\### Initialize vectors

Pvalue <- rep(1, nrow(vsgv))
HC_mean <- rep(1, nrow(vsgv))
Advanced_mean <- rep(1, nrow(vsgv))

\### Compute Wilcoxon rank sum test and save results

for (i in 1:nrow(vsgv)) {
  y <- try(wilcox.test(as.numeric(vsgv[i, subset(meta_all, Type2 == "HC")[,"ID"]]),
                       as.numeric(vsgv[i, subset(meta_all, Type2 == "Advanced")[,"ID"]]), paired = FALSE))
  try(Pvalue[i] <- y$p.value)
  HC_mean[i] <- as.numeric(summary(as.numeric(vsgv[i, subset(meta_all, Type2 == "HC")[,"ID"]]))[4])
  Advanced_mean[i] <- as.numeric(summary(as.numeric(vsgv[i, subset(meta_all, Type2 == "Advanced")[,"ID"]]))[4])
}

\### Output results

out <- cbind(rownames(vsgv), Pvalue, HC_mean, Advanced_mean, rep("all", nrow(vsgv)))
write.table(out, file = "vsgv.wilcox.HC_Advanced.all.txt", quote = FALSE, sep = "\t", row.names = FALSE)

 

##### 8.3 DSGV differential analysis

 

\### This code calculates the sum of each rows withing HC_1, CRC_1, 

\### HC_0, CRC_0, and all, before running a chisq or a fisher test on 

\### each row, depending on the values found. The type of the test run 

\### on each row is also defined, along with the p-value.

 

\```R

\### HC-Early dSV Differential Analysis

\### Read in files, subset data, and transpose the dsgv matrix

dsgv <- read.table("ALL.dsgv.name_1.csv", sep = ",", header = TRUE, row.names = 1)
meta <- read.table("all.data.metadata1211_Early0131_all.txt", sep = "\t", header = TRUE)

\### Subset data

meta_all <- subset(meta, Type2 != "CRA" & Type2 != "Advanced")
meta_all$ID <- as.character(meta_all$ID)
dsgv <- dsgv[as.character(meta_all$ID),]
dsgv[is.na(dsgv)] <- 0
dsgv <- t(dsgv)
meta_all$ID == colnames(dsgv)

\### Initialize vectors

HCn <- nrow(subset(meta_all, Type2 == "HC"))
CRCn <- nrow(subset(meta_all, Type2 == "Early"))
alln <- nrow(meta_all)

\### Summary and checks

meta_all[1:HCn, 6]
meta_all[(HCn + 1):alln, 6]
summary(meta_all[1:HCn, 6] == "HC")
summary(meta_all[(HCn + 1):alln, 6] == "Early")
colnames(dsgv[, 1:HCn]) == meta_all[1:HCn, "ID"]
colnames(dsgv[, (HCn + 1):alln]) == meta_all[(HCn + 1):alln, "ID"]

\### Calculate sums and create result matrix

HC_1 <- rowSums(dsgv[, 1:HCn])
CRC_1 <- rowSums(dsgv[, (HCn + 1):alln])
HC_0 <- HCn - rowSums(dsgv[, 1:HCn])
CRC_0 <- CRCn - rowSums(dsgv[, (HCn + 1):alln])
rst <- cbind(HC_1, CRC_1, HC_0, CRC_0, all = HC_1 + CRC_1 + HC_0 + CRC_0)

\### Summary of results

summary(rst[, 3] < 5)
summary(rst[, 4] < 5)

\### Generate counts

Pvalue <- rep(1, nrow(rst))
Type <- rep(1, nrow(rst))

\### Perform hypothesis tests

for (i in 1:nrow(rst)) {
  if ((rst[i, 1] >= 5) && (rst[i, 2] >= 5) && (rst[i, 3] >= 5) && (rst[i, 4] >= 5) && (rst[i, 5] >= 40)) {
    y = chisq.test(as.data.frame(rbind(c(rst[i, 1], rst[i, 2]), rbind(c(rst[i, 3], rst[i, 4])))), correct = FALSE)
    x = "chisq"
  } else if ((rst[i, 1] < 1) || (rst[i, 2] < 1) || (rst[i, 3] < 1) || (rst[i, 4] < 1) || (rst[i, 5] < 40)) {
    y = fisher.test(as.data.frame(rbind(c(rst[i, 1], rst[i, 2]), rbind(c(rst[i, 3], rst[i, 4])))))
    x = "fisher"
  } else if ((((rst[i, 1] >= 1) && (rst[i, 1] < 5)) || ((rst[i, 2] >= 1) && (rst[i, 2] < 5)) || 
              ((rst[i, 3] >= 1) && (rst[i, 3] < 5)) || ((rst[i, 4] >= 1) && (rst[i, 4] < 5))) && (rst[i, 5] >= 40)) {
    y = chisq.test(as.data.frame(rbind(c(rst[i, 1], rst[i, 2]), rbind(c(rst[i, 3], rst[i, 4])))), correct = TRUE)
    x = "chisq corrected"
  }
  try(Pvalue[i] <- y$p.value)
  Type[i] <- x
}

\### Output results

out <- cbind(rownames(dsgv), rst, Pvalue, rep("all:HC_Early", nrow(rst)), Type)
write.table(out, file = "dsgv.chisq.HC_Early0131.txt", quote = FALSE, sep = "\t", row.names = FALSE)

\### HC-Advanced

\### Repeat similar steps for HC-Advanced

\```

 

##### 8.4 Upset plot

\```R

\### Load required libraries

library(UpSetR)

\### Plot for DSGV

df_dsgv <- read.table("upset.DE.dsgv.plot.txt", sep = "\t", header = TRUE, row.names = 1)
pdf("DSGV.upset.pdf", height = 6, width = 12)
upset(df_dsgv, point.size = 3, line.size = 1, mainbar.y.label = "Count of Intersection", sets.x.label = "Differential genus in different sites")
dev.off()

\### Plot for VSGV

df_vsgv <- read.table("ALL.DE.dsgv.sort.upset.plot1.txt", sep = "\t", header = TRUE, row.names = 1)
pdf("VSGV.upset.pdf", height = 6, width = 12)
upset(df_vsgv, point.size = 3, line.size = 1, mainbar.y.label = "Count of Intersection", sets.x.label = "Differential genus in different sites")
dev.off()

\```

 

 