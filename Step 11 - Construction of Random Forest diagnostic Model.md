**### Step 11 - Construction of Random Forest diagnostic Model**

\### Construction and evaluation of the CRC diagnostic model based on microbial 

\### abundance composition and SVs signatures. The study begins with a thorough 

\### investigation of potential microbial markers for diagnosing colorectal cancer. 

\### This investigation involves analyzing microbial signatures, including microbial

\### abundance and SVs.

 

\### Model Construction and Evaluation: Several steps are taken to construct and 

\### evaluate the diagnostic model:

\### **Cross-Validation**: Cross-validation models are constructed and evaluated to

 \### ensure robustness.

\### **Cohort-to-Cohort Comparisons**: The study involves comparing different 

\### patient cohorts to assess the model's performance across different populations.

\### **Leave-One-Cohort-Out (LOCO) Evaluation**: LOCO evaluation is conducted, 

\### where one cohort's data is left out for testing while the others are used for training.

\### **Independent Validation**: Independent validation is performed to validate the

\### model's performance on new data.

 

\### The random forest algorithm is selected for this study due to its superior 

\### performance and ability to rank feature importance.

 

\### Cohort Definitions: The study defines several cohorts based on the year of 

\### sampling. Some cohorts are used for discovery, while others are for validation.

 

\### Microbial Characteristics: Model construction is guided by specific microbial 

\### characteristics, which include differential species abundance composition 

\### and structural variation between early-stage CRC patients and healthy individuals.

 

\### Input Characteristics: Three sets of input characteristics are used for the 

\### discovery cohorts, including differential species composition, 

\### differential structural variation characteristics, and the integration of both.

 

\### Feature Selection: The top 20 important features, ranked based on the 

\### importance index of the random forest algorithm, are selected as biomarkers 

\### for constructing a model of microbial community characteristics for each cohort.

 

\### Model Construction: Random forest is used to construct a classified 

\### model for microbial community characteristics. To avoid overfitting, a 

\### 5-fold cross-validation approach is employed.

 

\### Model Evaluation: Model performance is assessed using various metrics, 

\### including area under the curve (AUC), accuracy, sensitivity, and specificity.

 

\### Classifier Application: The classifier constructed using data from all 

\### populations within the queue is applied to independent test cohorts to 

\### assess its classification and correction performance. Receiver Operating 

\### Characteristic (ROC) curves are generated to evaluate performance.

 

\### Stability and Reliability Assessment: The study aims to assess the 

\### stability and reliability of microbial-based classifiers across different populations 

\### and technical variations. This involves inter-study cross-validation, where the

\### classifier is trained on one dataset and evaluated on others, and 

\### leave-one-country-out (LOCO) validation, where one country's data serves as 

\### the test set.

 

\### Generalizability: Through rigorous validation approaches, the study aims 

\### to determine the generalizability and robustness of the classifiers in 

\### different population settings and technical conditions.

 

\### Overall, this process aims to develop a robust and accurate diagnostic 

\### model for colorectal cancer based on microbial characteristics and 

\### evaluate its performance across diverse cohorts and populations.

 

##### 4.1 AUCRF

\```R

library(AUCRF)

\### Read metadata

meta <- read.table("metadata_Early.no10GER.txt", sep = "\t", header = TRUE)
meta$ID <- as.character(meta$ID)
meta$Study <- as.character(meta$Study)
meta$Stage <- as.character(meta$Type2)
meta01 <- meta

\### Read data

data.all <- read.table("ALL.dsgv.name.Early.group.txt", header = TRUE, row.names = 1, sep = "\t")
data.all <- data.all[meta01[, 1], ]

\### Read DE data

DE_HC_CRC <- read.table("dsgv.chisq.HC_Early.allno10GER.DE.txt", sep = "\t", header = FALSE)

\### Subset data based on DE results

data.all1 <- data.all[, c("Type2", as.character(DE_HC_CRC[, 1]))]
data.HC_CRC <- subset(data.all1, Type2 != "CRA")
data.HC_CRC[is.na(data.HC_CRC)] <- 0
data.HC_CRC$Type2 <- as.factor(data.HC_CRC$Type2)

\### Fit AUCRF model

fit <- AUCRF(Type2 ~ ., data = data.HC_CRC)

\### Write results to files

write.table(as.data.frame(fit$AUCcurve), "all.DSGV.AUCRF.AUCcurve.Early.txt", sep = "\t", quote = FALSE)
write.table(as.data.frame(fit$ranking), "all.DSGV.AUCRF.ranking.Early.txt", sep = "\t", quote = FALSE)

\```

 

##### 4.2 5-Fold cross validation

 

\```R

\### Load necessary libraries

library(tidyverse)
library(rfUtilities)
library(rfPermute)
library(randomForest)
library(pROC)
library(caret)

\### Read metadata

meta <- read.table("metadata_Early.no10GER.txt", sep = "\t", header = TRUE)
meta$ID <- as.character(meta$ID)
meta$Study <- as.character(meta$Study)
meta$Stage <- as.character(meta$Type2)
meta01 <- meta

\### Load data

data.all <- read.table("ALL.vsgv.dsgv.name.Early.group.txt", header = TRUE, row.names = 1, sep = "\t")
data.all <- data.all[meta01[, 1], ]

\### Read DE_HC_CRC

DE_HC_CRC <- read.table("dsgv.chisq.HC_Early.allno10GER.DE.txt", sep = "\t", header = FALSE)

\### Subset data

data.all1 <- data.all[, c("Type2", as.character(DE_HC_CRC[, 1]))]
data.HC_CRC <- subset(data.all1, Type2 != "CRA")
data.HC_CRC[is.na(data.HC_CRC)] <- 0

\### Check dimensions

rownames(data.HC_CRC) == meta01$ID
dim(meta)
dim(data.HC_CRC)

\### Build random forest model

RF <- randomForest(as.factor(Type2) ~ ., data = data.HC_CRC, ntree = 1000, proximity = TRUE, importance = TRUE)

\### Plot random forest

plot(RF)

\### Perform 5-fold cross-validation

folds <- createFolds(y = data.HC_CRC[, 1], k = 5)
max = 0
num = 0
fc = numeric()
mod_pre = numeric()

for (i in 1:5) {
  fold_test <- data.HC_CRC[folds[[i]], ]
  fold_train <- data.HC_CRC[-folds[[i]], ]

  model <- randomForest(as.factor(Type2) ~ ., data = fold_train, ntree = 1500, proximity = TRUE, importance = TRUE)
  model_pre <- predict(model, newdata = fold_test, type = "prob")

  fc <- c(fc, as.numeric(as.factor(fold_test$Type2)))
  mod_pre <- c(mod_pre, model_pre[, 1])
}

\### Combine results and plot ROC curve

df <- data.frame(fc, as.numeric(mod_pre))

x <- plot.roc(df[, 1], df[, 2],
              smooth = FALSE,
              lwd = 2,
              ylim = c(0, 1),
              xlim = c(1, 0),
              legacy.axes = TRUE,
              main = "",
              col = "red",
              print.auc = TRUE)

roc.rst <- roc(df[, 1], df[, 2])
roc.rst$auc

\### Write results

write.table(as.data.frame(RF$importance), "all.DE.importance.Early.txt", sep = "\t", quote = FALSE)
write.table(df_max, "all.df_max.Early.txt", sep = "\t", quote = FALSE)

\### Perform rfPermute

richness_rfP <- rfPermute(as.factor(Group) ~ ., data = data.HC_CRC, ntree = 3000, na.action = na.omit, nrep = 100, num.cores = 20)

\### Plot ROC curve and write results

roc <- roc(data.HC_CRC$Group, richness_rfP$rf$votes[, 1], percent = TRUE)
cat("AUC", "Run", gname, "\n", file = name_roc, sep = "\t", fill = FALSE, labels = NULL, append = TRUE)
cat(roc$auc, i, gname, "\n", file = name_roc, sep = "\t", fill = FALSE, labels = NULL, append = TRUE)

Extract and write importance results

imp <- as.data.frame(importance(richness_rfP))
imp1 <- imp[which(imp$MeanDecreaseAccuracy.pval < 1), c("MeanDecreaseAccuracy", "MeanDecreaseAccuracy.pval")]
imp1_rst1 <- data.frame(imp1, Run = rep(i, times = nrow(imp1)), Group_type = rep(gname, nrow(imp1)))

\```

 

##### 4.3 Random Sampling cross validation

\```R

library(tidyverse)

library("randomForest")

library("rfUtilities")

library("rfPermute")

library("pROC")

library(randomForest)

library(caret)

library(pROC)

library(caret)

 

\# Load metadata, library and table

meta<-read.table("metadata_Early.no10GER.txt",sep="\t",header=T)

meta$ID<-as.character(meta$ID)

meta$Study<-as.character(meta$Study)

meta$Group<-as.character(meta$Type2)

data.all<-read.table("ALL.vsgv.dsgv.name.Early.group.txt",header=T,row.names=1,sep="\t")

data.all<-data.all[meta[,1],]

DE_HC_CRC<-read.table("Early.DE2.test",sep="\t",header=F)

data.all1<-data.all[,c("Type2",as.character(DE_HC_CRC[,1]))]

data.HC_CRC<-subset(data.all1,Type2!="CRA")

data.HC_CRC[is.na(data.HC_CRC)] <- 0

rownames(data.HC_CRC)==meta$ID

dim(meta)

dim(data.HC_CRC)

 

\# Create cross validation 

library(plyr)

library(randomForest)

CVgroup <- function(k,datasize,seed){

 cvlist <- list()

 set.seed(seed)

 n <- rep(1:k,ceiling(datasize/k))[1:datasize] # Divide data into K parts and generate a complete dataset n

 temp <- sample(n,datasize) # Shuffle n

 x <- 1:k

 dataseq <- 1:datasize

 cvlist <- lapply(x,function(x) dataseq[temp==x]) # Randomly generate k ordered data sequences in dataseq

 return(cvlist)

}

 



\# Start training random forest models 

ptm <- proc.time() # Start running timer

gname<-"VSGV.HC_CRC.Discovery"

name_roc=paste(gname,".ROC.DE2.test.txt",sep="")

name_imp1=paste(gname,".imp1.DE2.test.txt",sep="")

name_confu=paste(gname,".confu.DE2.test.txt",sep="")

\# Create the ROC.DE2.test.txt, imp1.DE2.test.txt andconfu.DE2.test.txt files 

 

\# Create the model 

L=1

a1<-sample(ncol(data.HC_CRC[,2:21]),16,replace=F) # Choose measurements randomly

a1<-a1+1

data.HC_CRC1<-data.HC_CRC[,c(1,a1)]

data <- data.HC_CRC1 

pred <- data.frame() 

k<-5 # Number of training processes 

m <- seq(1000,1500,by = 100) # Sequence of number of trees used 

model1 <-randomForest(as.factor(Type2)~.,data = data,ntree = 1200) # Build random forest using given parameters   

for(j in m){ 

 progress.bar <- create_progress_bar("text") # Show progress

 progress.bar$init(k) # Show progress

 for (i in 1:k){

 train <- data[-cvlist[[i]],] # Train data using cross-validation 

 test <- data[cvlist[[i]],] # Test data

 model <-randomForest(as.factor(Type2)~.,data = train,ntree = j) # Build random forest accordingly 

 prediction <- predict(model,subset(test,select = -Type2)) # Predict Type2 of test data 

 randomtree <- rep(j,length(prediction)) # Store the number of trees used in each prediction 

 kcross <- rep(i,length(prediction)) # Store the number of cross-validations used per prediction 

\# Produce a confusion Matrix of the predictions made by the Random Forest model, using the test data set 

confu <- table(subset(test, select = Type2), prediction)

\# Write the matrix to an external txt file 

write.table(confu, name_confu, sep="\t",quote=F)

proc.time()-ptm

\```

 

 