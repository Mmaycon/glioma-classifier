
# Packages ----------------
library(dplyr)
library(data.table)
library(stringr)
library(matlab)
library(pheatmap)
library(caret)
library(randomForest)
library(doMC)
library(e1071)
library(openxlsx)
library(readr)



### Prepare the RF input ### 
### Prepare the RF input ###
### Prepare the RF input ### 


# 1. load DNAmtx with MVF
DNAmtx <- readRDS("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DNAmtx_glioma_TCGA_proc_NAreplaced_MostVariableFeature.rds")
dim(DNAmtx) # 28934   691
DNAmtx %>% dim() #28934   691
rownames(DNAmtx) <- DNAmtx$probeID
DNAmtx$probeID <- NULL
DNAmtx$Variance <- NULL
DNAmtx <- as.matrix(DNAmtx)
colnames(DNAmtx) <- str_split_fixed(as.character(colnames(DNAmtx)), "[-][0-9][0-9][A-Z]", 2)[,1] #to match to metadata sample ID


# 2. laod metadata
load("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/from_local_mc/metadata.rda")
dim(metadata) #878  52
metadata[1:4, 1:4]


# 3. correct samples intersected in both metadata and DNAmtx
table(duplicated(colnames(DNAmtx)))
# FALSE  TRUE
#   658    31
DNAmtx <- DNAmtx[, !duplicated(colnames(DNAmtx))]
table(duplicated(colnames(DNAmtx))) # no duplcates anymore
# Adjusting samples on both metada and DNAmtx data
metadata <- metadata[metadata$Case %in% colnames(DNAmtx), ]
length(metadata$Case) #641
DNAmtx <- DNAmtx[, colnames(DNAmtx) %in% metadata$Case]
length(colnames(DNAmtx)) #641


# 4. load the probset chosen
probeset <- readRDS("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/Extreme_upset_cutoff_Glioma_probes_vector.rds")

# Mount RF input data --------------

# Select 70% of the data to train the model (use 'Ceccarelli_six_subtypes' label to do that)
library(caret)
set.seed(666)
meta_70per <- createDataPartition(metadata$Ceccarelli_six_subtypes, p=0.70, list=FALSE, times=1)
meta_70per <- metadata[meta_70per,]

meta_30per <- metadata[!metadata$Case %in% meta_70per$Case, ]

# Filter DNAmtx matrix over those 70% samples
DNAmet_70per <- DNAmtx[, colnames(DNAmtx) %in% meta_70per$Case]
DNAmet_30per <- DNAmtx[, !colnames(DNAmtx) %in% meta_70per$Case]

# Keep only your probeset features
DNAmet_70per <- DNAmet_70per[rownames(DNAmet_70per) %in% probeset, ]
DNAmet_30per <- DNAmet_30per[rownames(DNAmet_30per) %in% probeset, ]

# Replace NA on DNAmtx - if it is necessary
table(is.na(DNAmet_70per))
# FALSE
# 167692
# library(impute)
# DNAmet_70per <- impute.knn(as.matrix(DNAmet_70per), k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)[[1]]
# table(is.na(DNAmet_70per))


table(is.na(DNAmet_30per))
# FALSE
# 70119
# library(impute)
# DNAmet_30per <- impute.knn(as.matrix(DNAmet_30per), k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)[[1]]
# table(is.na(DNAmet_30per))





# Add label columns into the DNAmet matrix
trainingdata <- t(DNAmet_70per)
trainingdata <- merge(trainingdata, meta_70per[,c("Case", "IDH.specific.DNA.Methylation.Cluster", "Supervised.DNA.Methylation.Cluster", 'Ceccarelli_six_subtypes','IDH.status')], by.x=0,by.y="Case", all.x=T)
rownames(trainingdata) <- as.character(trainingdata$Row.names)
trainingdata <- trainingdata[,-1]
trainingdata <- droplevels(trainingdata)

testingdata <- t(DNAmet_30per)
testingdata <- merge(testingdata, meta_30per[,c("Case", "IDH.specific.DNA.Methylation.Cluster", "Supervised.DNA.Methylation.Cluster", 'Ceccarelli_six_subtypes','IDH.status')], by.x=0,by.y="Case", all.x=T)
rownames(testingdata) <- as.character(testingdata$Row.names)
testingdata <- testingdata[,-1]
testingdata <- droplevels(testingdata)


# Replace NA - if it is necessary
table(is.na(trainingdata)); dim(trainingdata)
# FALSE
# 169500
# dim 452 375
# trainingdata <- na.omit(trainingdata) # If necessary: use na.omit to take out samples without label. It is removing samples because it's tranposed
table(is.na(trainingdata)); dim(trainingdata)
# FALSE
# 169500
# dim 452 375


table(is.na(testingdata)); dim(testingdata)
#FALSE
#70875
#dim  189 375
# testingdata <- na.omit(testingdata) # If necessary: use na.omit to take out samples without label. It is removing samples because it's tranposed
table(is.na(testingdata)); dim(testingdata)
# FALSE
# 70875
# dim 189 375


# saveRDS(trainingdata, file="/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/RF/trainingdata_70per_378pset_KnnStringestProbeset.rds")

# saveRDS(testingdata, file="/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/RF/testingdata_30per_378pset_KnnStringestProbeset.rds")







### Train RF ### 
### Train RF ###
### Train RF ### 

# Run it on terminal 
start_time <- Sys.time()

library(caret)
library(randomForest)
library(doMC)
library(e1071)
# Load the data input
trainingdata <- readRDS("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/RF/trainingdata_70per_378pset_KnnStringestProbeset.rds")
trainingdata[1:3, 1:3]
# delet any column but all the features you want plus the predictor variable you want (in this case it is 'Ceccarelli_six_subtypes')
trainingdata$IDH.specific.DNA.Methylation.Cluster <- NULL
trainingdata$Supervised.DNA.Methylation.Cluster <- NULL
trainingdata$IDH.status <- NULL
head(trainingdata)


# register cores for doMC
registerDoMC(cores = 4) # in Tathi's code it was set for 10
# set up k-fold cross validation
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)
# Obs: you may additionally, if you wish use a different method, for validating your model parameters, such as OOB (Out of Bag). OOB is faster.


# Set your seed so your work is repeatable
set.seed(42)

# set values for mtry
# mtry is the "Number of variables randomly sampled as candidates at each split"
# traditionally for classification you use the sqrt of the number of variables
# but here we try a range of mtry values to find the best parameters for our model - Maycon: we could work on this parameter to get better results 
mtryVals <- floor(c(seq(100, 2000, by=100),
                    sqrt(ncol(trainingdata))))
mtryGrid <- data.frame(.mtry=mtryVals)
# Confirm seed again
set.seed(666)
# Set number of cores
registerDoMC(cores = 4)
# Run Training
stemsig <- train(Ceccarelli_six_subtypes ~ ., # variable to be trained on
                 data = trainingdata, # Data we are using
                 method = "rf", # Method we are using
                 trControl = fitControl, # How we validate
                 # We created this object above
                 ntree = 5000, # number of trees
                 # is dependent on training data size
                 importance = TRUE, # calculate varible importance
                 # can be omitted to speed up calc
                 tuneGrid = mtryGrid # set mtrys
                 #subset = myTrain # define training set #comment when train with 100% of samples OR when your trainingdata is ready to go as it is in this case
)
end_time <- Sys.time()
runtime = end_time - start_time
print(runtime)
# Saving the entire environment cause we're running it on the terminal
# save(list=ls(),file="yout_path/RF_70per_378pset_KnnNAreplacedData_MVFplusDiffMeanProbeSelection.Rda")



### Back to Rstudio 

# Plot a confusion matrix and check it's AUC
# Load the Random Forest output (trained model)
load("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/RF/RF_70per_378pset_KnnStringestProbeset.Rda") #RF output objects (stemsig, runtime)
# Load the testing data - it can be either the remaining 30% out of the training data or an independent dataset
testingdata <- readRDS("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/RF/testingdata_30per_378pset_KnnStringestProbeset.rds") #testingdata

# Create a 'RFtype' column using the RF output (stemsig) to predict the label 'Ceccarelli_six_subtypes' into the testing data
testingdata$RFtype <- predict(stemsig, testingdata)

# Plot confusion matrix to get the acuraccy and other statisticla metrics 
confusionMatrix(data = as.factor(testingdata$RFtype), reference = as.factor(testingdata$Ceccarelli_six_subtypes)) # AUC = 0.9101 

table(trainingdata$Ceccarelli_six_subtypes)
table(testingdata$Ceccarelli_six_subtypes)



### GLAS gliomas ###
### GLAS gliomas ###
### GLAS gliomas ###

# GLAS gliomas
# Load DNAmtx and metadata
load("/media/hd/maycon/MESTRADO/GLASS_methy.matrix&pData/beta.picked.Rda") # GLASS gliomas which we should work with
beta.picked[1:4, 1:4]
dim(beta.picked) # 452453    143
colnames(beta.picked)
table(duplicated(colnames(beta.picked)))
# FALSE
# 143

table(is.na(beta.picked))
beta.picked <- beta.picked[rownames(beta.picked) %in% colnames(testingdata), ]
table(is.na(beta.picked))
# FALSE  TRUE
# 52654   399

# Replace probe NA
library(impute)
beta.picked <- impute.knn(as.matrix(beta.picked), k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)[[1]]

print(load("/media/hd/maycon/MESTRADO/GLASS_methy.matrix&pData/pData.GLASS.SI_nsc.gsc.neuron.Rda"))
pData.GLASS %>% head()
dim(pData.GLASS)
pData.GLASS$glassID
pData.GLASS$barcode <- pData.GLASS$glassID
meta_picked <- pData.GLASS[pData.GLASS$barcode %in% colnames(beta.picked), ]

meta_picked$Ceccarelli_six_subtypes <- meta_picked$TCGA.subtype
meta_picked[meta_picked$Ceccarelli_six_subtypes %in% "PA-like", ]$Ceccarelli_six_subtypes <- "LGm6-GBM"
meta_picked$Ceccarelli_six_subtype
meta_picked <- droplevels(meta_picked)

# Structure DNAmtx as in testingdata
# Eg:
#              cg27323784 cg27638126 cg27659841 subtype
# TCGA-06-1804 0.02649859 0.02640942  0.7308063 codel
# TCGA-06-5413 0.02888786 0.03410674  0.4588765 classic
meta_to_merge <- meta_picked %>%
  select(c('Ceccarelli_six_subtypes', 'barcode'))

testingdata_GLASS <- t(beta.picked)
testingdata_GLASS <- data.frame(testingdata_GLASS)
testingdata_GLASS$barcode <- rownames(testingdata_GLASS)
testingdata_GLASS <- merge(testingdata_GLASS, meta_to_merge, by.y = 'barcode')
rownames(testingdata_GLASS) <- testingdata_GLASS$barcode
testingdata_GLASS$barcode <- NULL
testingdata_GLASS <- testingdata_GLASS[, colnames(testingdata_GLASS) %in%  colnames(testingdata)]
dim(testingdata_GLASS) #143 372 (missing 4 probes from TCGA)

table(is.na(testingdata_GLASS))
# FALSE
# 53196

# Apply RF output (Glioma Classifier) on GLASS
testingdata_GLASS$RFtype <- predict(stemsig, testingdata_GLASS)

# Putting column levels in the same order (refactoring levels order)
testingdata_GLASS$RFtype <- factor(testingdata_GLASS$RFtype, 
                                   levels = c(names(table(testingdata_GLASS$Ceccarelli_six_subtypes))))

# Plot confusion matrix
confusionMatrix(data = as.factor(testingdata_GLASS$RFtype), reference = as.factor(testingdata_GLASS$Ceccarelli_six_subtypes))
# AUC = 0.9021







