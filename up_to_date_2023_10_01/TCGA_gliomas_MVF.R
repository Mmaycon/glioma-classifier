# Part 1
#**TCGA Gliomas DNAmet 450k - most variable features selection**#

# Obs:
# 1. Saving processed TCGA gliomas DNAmet samples within most variable feature
# 2. I processed these files most in my masters because the processed files from TCGA were with a lot of NAs on their probes 

# load metadata
load("/media/hd/maycon/MESTRADO/SAVE.point/matrix.methy_pdata_global/dataset.predicao.25.08.2021/pData.prediction_4801_atualizado_30.11.2021.Rda")
dim(pData.prediction_4801) #4801   15

# load matrix meth
load("/media/hd/maycon/MESTRADO/SAVE.point/matrix.methy_pdata_global/dataset.predicao.25.08.2021/matrix.meth.global_tcgaNormalized.TCGA_ID_recovered.24.01.2022.GBM.tcga.CORRETO.Rda")
dim(matrix.meth.global_tcgaNormalized_4801) #452832   4801

# keep only TCGA gliomas 
glioma_TCGA_proc <- pData.prediction_4801[pData.prediction_4801$Class.added %in% c('LGG_tcga.proc',
                                                                                   'GBM_tcga.proc'), ]
meta_glioma_TCGA_proc <- glioma_TCGA_proc[, c('barcode', 'Class.added')]

DNAmtx_glioma_TCGA_proc <- matrix.meth.global_tcgaNormalized_4801[, colnames(matrix.meth.global_tcgaNormalized_4801) %in% glioma_TCGA_proc$barcode]
dim(DNAmtx_glioma_TCGA_proc) #452832    689

# save DNAmet matrix into another object
beta_values <- DNAmtx_glioma_TCGA_proc

# remove mask/chr probes
load("/media/hd/maycon/Glioma_classifier/hm450.anno.Rda") #genome annotation object
probes_retain <- subset(hm450.anno, !chrm_A %in% c("chrX","chrY", "chrM") & MASK_general == FALSE)$probeID
beta_values = beta_values[rownames(beta_values) %in% probes_retain, ]

# Replacing NAs
library(impute)
beta_values <- impute.knn(as.matrix(beta_values), k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)[[1]] # it takes some time (like 5 min)
# saveRDS(beta_values, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DNAmtx_glioma_TCGA_proc_NAreplacedKNN.rds') #save it in along the process becasue it takes some time to run

# Calculate the variance of beta values for each CpG site
beta_values <- readRDS('/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DNAmtx_glioma_TCGA_proc_NAreplacedKNN.rds')
dim(beta_values) #383940    689
variance_values <- apply(beta_values, 1, var)

# Create a dataframe with CpG sites and their variance values
variance_data <- data.frame(CpG_Site = rownames(beta_values), Variance = variance_values)

# Visualize all probes on the elbow plot
library(ggplot2); theme_set(theme_classic())
elbow_plot <- variance_data %>%
  arrange(desc(Variance)) %>%
  mutate(Rank = row_number()) %>%
  ggplot(aes(x = Rank, y = Variance)) +
  geom_line() +
  geom_point() +
  labs(title = "Elbow Plot for Variance Cutoff",
       x = "Number of CpG Sites",
       y = "Variance") +
  theme_minimal(); elbow_plot # 0.125 seems a good cutoff 

variance_data[variance_data$Variance >=  0.125, ] %>% dim() #310 probes 
variance_data[variance_data$Variance >=  0.100, ] %>% dim() #1821    probes
variance_data[variance_data$Variance >=  0.050, ] %>% dim() #28934     probes. It's a decent amount of probes to move on

# Select the top N variable features (e.g., top 100)
library(dplyr)
top_n <- 28934
selected_features <- variance_data %>%
  arrange(desc(Variance)) %>%
  slice(1:top_n)

dim(selected_features)
head(selected_features)
names(selected_features) <- c('probeID', 'Variance')
beta_values <- beta_values[rownames(beta_values) %in% rownames(selected_features), ]
beta_values <- as.data.frame(beta_values)
beta_values$probeID <- rownames(beta_values)
DNAmtx <- merge(beta_values, selected_features, by.y = 'probeID')
DNAmtx[1:4, 1:4]

# saveRDS(DNAmtx, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DNAmtx_glioma_TCGA_proc_NAreplaced_MostVariableFeature.rds')



# Part 2
#** Fine tune the MostVariableFeature (MVF) by DiffMean **#
library(dplyr)
library(stringr)
library(UpSetR)

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

#** NOTE **#
# All of these MVF seems to be significant different accordingly to dmpFinder() - from {minfi} - & this set of probes made a great clustering on t-SNE over the 6 subtypes. That's why we're moving foward without a "Differential DNA methylation" 




# 4. Find DiffMean group by group (subtype 1 vs subtype 2; subtype 1 vs subtype 3 so on) in order to get markers for each subtype

# Automating the DiffMean calculation

# testing the loop 
# i <- 1
# condition_one <- 'Classic-like'
# data <- DNAmtx
# metadata <- metadata
# threshold <- 0.3

# Run this part for the function work out
group_names <- names(table(metadata$Ceccarelli_six_subtypes))
calculate_diff_mean <- function(data, 
                                metadata, 
                                condition_one, 
                                threshold) { #considering all probes are statistical signifcant for us (in this case it is). in other cases it would need a "statistical test" to filter the probes/features. I didn't put all in one function because sometimes people prefer other kind o statistical tests than which I use.
  list_diffmean_dfs <- list()
  data <- as.data.frame(data)
  condition_all_but_one <- group_names[!group_names %in% condition_one]
  ###condition_all_but_one <- 'Mesenchymal-like' #TESTING only. 
  for(i in 1:length(condition_all_but_one)) {
    print(paste0("Processing group: ", condition_all_but_one[i]))
    
    try({ #it ignores any error that might prevent the code of going forward
      # But it is turned off in this script on 'silent = FALSE')
      loop_data <- data
      condition_1_sampleID <- metadata[metadata$Ceccarelli_six_subtypes %in% condition_one, ]$Case
      condition_2_sampleID <- metadata[metadata$Ceccarelli_six_subtypes %in% condition_all_but_one[i], ]$Case
      
      print(paste0("Length of condition_1_sampleID: ", length(condition_1_sampleID)))
      print(paste0("Length of condition_2_sampleID: ", length(condition_2_sampleID)))
      
      loop_data$meanM1 <- apply(loop_data[, condition_1_sampleID], 1, mean, na.rm = TRUE)
      loop_data$meanM2 <- apply(loop_data[, condition_2_sampleID], 1, mean, na.rm = TRUE)
      loop_data$DiffMean <- loop_data$meanM1 - loop_data$meanM2
      loop_data$Comparison <- paste0(condition_one, '_', condition_all_but_one[i], '_', condition_one, '_', 'Orientation')
      
      
      
      loop_data$DNAmet_orientation <- NA
      #loop_data$DNAmet_orientation <- as.character(loop_data$DNAmet_orientation) #this time it's not necessary
      if(dim(loop_data[loop_data$DiffMean > threshold, ])[1] > 0) {
        loop_data[loop_data$DiffMean > threshold, ]$DNAmet_orientation <- 'hyper'
      } else {
        # do nothing
      }
      
      if(dim(loop_data[loop_data$DiffMean < -threshold, ])[1] > 0) {
        loop_data[loop_data$DiffMean < -threshold, ]$DNAmet_orientation <- 'hypo'
      } else {
        # do nothing
      }
      
      loop_data[loop_data$DNAmet_orientation %in% NA, ]$DNAmet_orientation <- 'not_diff'
      
      
      
      loop_data$probeID <- rownames(loop_data)
      
      list_diffmean_dfs[[i]] <- loop_data[, c('DiffMean', 'Comparison', 'DNAmet_orientation', 'probeID')]
      print(paste0(list_diffmean_dfs[[i]]$Comparison[1], ' has been stored into the list.'))
    },  silent = FALSE)
  }
  return(list_diffmean_dfs) 
}


# Classic-like DNAmet markers ---------------
list_diffmean_dfs <-calculate_diff_mean(
  data = DNAmtx,
  metadata = metadata,
  condition_one = group_names[1], #"Classic-like" 
  threshold = 0.3)

Classiclike_markers <- do.call('rbind', list_diffmean_dfs)
head(Classiclike_markers)
length(names(table(Classiclike_markers$Comparison))) # 5 (it should be 5) 
# Attention 2: it's okay to have more probes than the ~28k MVF because in this dataframe should be 5 comparisons. So it's the sum of 5*~28k probes at total.
to_upsetplot <- Classiclike_markers[Classiclike_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  Classic_Codel = to_upsetplot[to_upsetplot$Comparison %in% 'Classic-like_Codel_Classic-like_Orientation', ]$probeID,
  Classic_GCIMPhigh = to_upsetplot[to_upsetplot$Comparison %in% 'Classic-like_G-CIMP-high_Classic-like_Orientation', ]$probeID,
  Classic_GCIMPlow = to_upsetplot[to_upsetplot$Comparison %in% 'Classic-like_G-CIMP-low_Classic-like_Orientation', ]$probeID,
  Classic_LGm6GBM = to_upsetplot[to_upsetplot$Comparison %in% 'Classic-like_LGm6-GBM_Classic-like_Orientation', ]$probeID,
  Classic_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'Classic-like_Mesenchymal-like_Classic-like_Orientation', ]$probeID
)

library(UpSetR)
upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

Classic_probeset_70p <- x1[ rowSums(x$New_data) == 5]


# Codel DNAmet markers ---------------
list_diffmean_dfs <- calculate_diff_mean(
  data = DNAmtx,
  metadata = metadata,
  condition_one = group_names[2], #"Codel" 
  threshold = 0.3)

Codel_markers <- do.call('rbind', list_diffmean_dfs)
head(Codel_markers)
length(names(table(Codel_markers$Comparison))) # 5 (it should be 5)

to_upsetplot <- Codel_markers[Codel_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  Codel_Classic = to_upsetplot[to_upsetplot$Comparison %in% 'Codel_Classic-like_Codel_Orientation', ]$probeID,
  Codel_GCIMPhigh = to_upsetplot[to_upsetplot$Comparison %in% 'Codel_G-CIMP-high_Codel_Orientation', ]$probeID,
  Codel_GCIMPlow = to_upsetplot[to_upsetplot$Comparison %in% 'Codel_G-CIMP-low_Codel_Orientation', ]$probeID,
  Codel_LGm6GBM = to_upsetplot[to_upsetplot$Comparison %in% 'Codel_LGm6-GBM_Codel_Orientation', ]$probeID,
  Codel_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'Codel_Mesenchymal-like_Codel_Orientation', ]$probeID
)

upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

Codel_probeset_147p <- x1[ rowSums(x$New_data) == 5]




# G-CIMP-high DNAmet markers ---------------
list_diffmean_dfs <- calculate_diff_mean(
  data = DNAmtx,
  metadata = metadata,
  condition_one = group_names[3], #"G-CIMP-high" 
  threshold = 0.3)

GCIMPhigh_markers <- do.call('rbind', list_diffmean_dfs)
head(GCIMPhigh_markers)
length(names(table(GCIMPhigh_markers$Comparison))) # 5 (it should be 5)

to_upsetplot <- GCIMPhigh_markers[GCIMPhigh_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  GCIMPhigh_Classic = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-high_Classic-like_G-CIMP-high_Orientation', ]$probeID,
  GCIMPhigh_Codel = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-high_Codel_G-CIMP-high_Orientation', ]$probeID,
  GCIMPhigh_GCIMPlow = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-high_G-CIMP-low_G-CIMP-high_Orientation', ]$probeID,
  GCIMPhigh_LGm6GBM = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-high_LGm6-GBM_G-CIMP-high_Orientation', ]$probeID,
  GCIMPhigh_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-high_Mesenchymal-like_G-CIMP-high_Orientation', ]$probeID
)

upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

CGIMPhigh_probeset_13p <- x1[ rowSums(x$New_data) == 5]


# G-CIMP-low DNAmet markers ---------------
list_diffmean_dfs <- calculate_diff_mean(
  data = DNAmtx,
  metadata = metadata,
  condition_one = group_names[4], #"G-CIMP-low" 
  threshold = 0.3)

GCIMPlow_markers <- do.call('rbind', list_diffmean_dfs)
head(GCIMPlow_markers)
length(names(table(GCIMPlow_markers$Comparison))) # 5 (it should be 5)

to_upsetplot <- GCIMPlow_markers[GCIMPlow_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  GCIMPlow_Classic = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-low_Classic-like_G-CIMP-low_Orientation', ]$probeID,
  GCIMPlow_Codel = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-low_Codel_G-CIMP-low_Orientation', ]$probeID,
  GCIMPlow_GCIMPhigh = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-low_G-CIMP-high_G-CIMP-low_Orientation', ]$probeID,
  GCIMPlow_LGm6GBM = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-low_LGm6-GBM_G-CIMP-low_Orientation', ]$probeID,
  GCIMPlow_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-low_Mesenchymal-like_G-CIMP-low_Orientation', ]$probeID
)

upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

CGIMPlow_probeset_129p <- x1[ rowSums(x$New_data) == 5]




# LGm6-GBM DNAmet markers ---------------
list_diffmean_dfs <- calculate_diff_mean(
  data = DNAmtx,
  metadata = metadata,
  condition_one = group_names[5], #"LGm6-GBM" 
  threshold = 0.3)

LGm6GBM_markers <- do.call('rbind', list_diffmean_dfs)
head(LGm6GBM_markers)
length(names(table(LGm6GBM_markers$Comparison))) # 5 (it should be 5)

to_upsetplot <- LGm6GBM_markers[LGm6GBM_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  LGm6GBM_Classic = to_upsetplot[to_upsetplot$Comparison %in% 'LGm6-GBM_Classic-like_LGm6-GBM_Orientation', ]$probeID,
  LGm6GBM_Codel = to_upsetplot[to_upsetplot$Comparison %in% 'LGm6-GBM_Codel_LGm6-GBM_Orientation', ]$probeID,
  LGm6GBM_GCIMPhigh = to_upsetplot[to_upsetplot$Comparison %in% 'LGm6-GBM_G-CIMP-high_LGm6-GBM_Orientation', ]$probeID,
  LGm6GBM_LGm6GBM = to_upsetplot[to_upsetplot$Comparison %in% 'LGm6-GBM_G-CIMP-low_LGm6-GBM_Orientation', ]$probeID,
  LGm6GBM_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'LGm6-GBM_Mesenchymal-like_LGm6-GBM_Orientation', ]$probeID
)

upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

LGm6GBM_probeset_12p <- x1[ rowSums(x$New_data) == 5]





# Mesenchymal-like DNAmet markers ---------------
list_diffmean_dfs <- calculate_diff_mean(
  data = DNAmtx,
  metadata = metadata,
  condition_one = group_names[6], #"Mesenchymal-like" 
  threshold = 0.3)

Mesenchymallike_markers <- do.call('rbind', list_diffmean_dfs)
head(Mesenchymallike_markers)
length(names(table(Mesenchymallike_markers$Comparison))) # 5 (it should be 5)

to_upsetplot <- Mesenchymallike_markers[Mesenchymallike_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  Mesenchymal_Classic = to_upsetplot[to_upsetplot$Comparison %in% 'Mesenchymal-like_Classic-like_Mesenchymal-like_Orientation', ]$probeID,
  Mesenchymal_Codel = to_upsetplot[to_upsetplot$Comparison %in% 'Mesenchymal-like_Codel_Mesenchymal-like_Orientation', ]$probeID,
  Mesenchymal_GCIMPhigh = to_upsetplot[to_upsetplot$Comparison %in% 'Mesenchymal-like_G-CIMP-high_Mesenchymal-like_Orientation', ]$probeID,
  Mesenchymal_GCIMPlow = to_upsetplot[to_upsetplot$Comparison %in% 'Mesenchymal-like_G-CIMP-low_Mesenchymal-like_Orientation', ]$probeID,
  Mesenchymal_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'Mesenchymal-like_LGm6-GBM_Mesenchymal-like_Orientation', ]$probeID
)

upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

Mesenchymal_probeset_7p <- x1[ rowSums(x$New_data) == 5]

Extreme_upset_cutoff_Glioma_probes_vector <- c(
  Classic_probeset_70p,
  Codel_probeset_147p,
  CGIMPhigh_probeset_13p,
  CGIMPlow_probeset_129p,
  LGm6GBM_probeset_12p,
  Mesenchymal_probeset_7p)


# saveRDS(Extreme_upset_cutoff_Glioma_probes_vector, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/Extreme_upset_cutoff_Glioma_probes_vector.rds')


