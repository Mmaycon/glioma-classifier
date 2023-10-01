#**MVF visualization of gliomas DNA methylation data**#

# Load DNAmet with all probes (~ 380k probes, mask/chr/NAs probes solved)
DNAmtx_whoArray <- readRDS('/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DNAmtx_glioma_TCGA_proc_NAreplacedKNN.rds')
library(stringr)
colnames(DNAmtx_whoArray) <- str_split_fixed(as.character(colnames(DNAmtx_whoArray)), "[-][0-9][0-9][A-Z]", 2)[,1] #fixing sample names
dim(DNAmtx_whoArray) #383940    689

# Load TCGA gliomas clinical metadata
load("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/from_local_mc/metadata.rda")
dim(metadata) # 878  52

# Filter glioma samples accordandly to metadata and DNAmtx 
colnames(DNAmtx_whoArray) <- str_split_fixed(as.character(colnames(DNAmtx_whoArray)), "[-][0-9][0-9][A-Z]", 2)[,1]
DNAmtx_whoArray <- DNAmtx_whoArray[, colnames(DNAmtx_whoArray) %in% metadata$Case]
metadata <- metadata[metadata$Case %in% colnames(DNAmtx_whoArray), ]
DNAmtx_whoArray <- DNAmtx_whoArray[, !duplicated(colnames(DNAmtx_whoArray))]
dim(DNAmtx_whoArray) #383940    641
dim(metadata) #641  52

# Load Most Variable Features fine tunned 
tuned_probes <- readRDS('/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/Extreme_upset_cutoff_Glioma_probes_vector.rds')
DNAmtx_tProbes <- DNAmtx_whoArray[rownames(DNAmtx_whoArray) %in% tuned_probes, ]
head(DNAmtx_tProbes)
dim(DNAmtx_tProbes) #371 641


# Plot DNAmtx_whoArray PCA ----------
pca_data <- DNAmtx_whoArray
dim(pca_data) # 383940   641
table(is.na(pca_data))
# FALSE      
# 264534660

pca2 <- prcomp(t(pca_data)) # this is the part where you calculate the PCA
aux <- as.data.frame(pca2$x[, 1:3]) #get only the info you need
pca_metadata <- merge(metadata, aux, by.y=0, by.x="Case", all.x=T)
# saveRDS(pca_metadata, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/glioma_TCGA_wholearray_EPIC450k_PCA_out.rds')

library(viridis)
n_colors <- length(table(pca_metadata$Ceccarelli_six_subtypes))
pal <- viridis(n = n_colors, option = "K", direction = -1)


library(ggplot2); theme_set(theme_classic())
ggplot(pca_metadata, aes(x=PC1, y=PC2, colour=Ceccarelli_six_subtypes)) +
  geom_point() +
  scale_color_manual(values=c("orange","green", "purple", "red","darkgreen","darkblue"), name="Legend") +
  #scale_fill_manual(values=c(pal), name="Legend") + 
  #scale_color_manual(values=c(pal), name="Legend") +
  #scale_color_viridis(discrete=TRUE) +
  xlab(paste0("PC1 (",prettyNum(summary(pca2)$importance[2,1]*100, digits = 2),"%)")) +
  ylab(paste0("PC2 (",prettyNum(summary(pca2)$importance[2,2]*100, digits = 2 ),"%)")) +
  ggtitle("Gliomas TCGA - EPIC/450k genome wide") +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
  ) 

ggsave(file='/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/glioma_TCGA_wholearray_EPIC450k_PCA_out.pdf',width = 12, height = 6, dpi = 300, 
) 


# t-SNE 
library(Rtsne)
tsne_realData <- Rtsne(t(pca_data), perplexity=30, check_duplicates = FALSE) # #function to run t-sne
### Plotando o tnse 
# add as 2 colunas (V1 e V2) com os dados do tsne
pdata.teste.tsne <- metadata #pData inicial
pdata.teste.tsne$V1 <- tsne_realData$Y[,1]
pdata.teste.tsne$V2 <- tsne_realData$Y[,2]

library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne[, ], aes(x=V1, V2, colour=Ceccarelli_six_subtypes)) +
  geom_point() +
  scale_color_manual(values=c("orange","green", "purple", "red","darkgreen","darkblue"), name="Legend") +
  # scale_color_viridis(discrete=TRUE) +
  
  ggtitle(label = "t-SNE Genome wide (380k probes)",
          subtitle = "perplexity=30")

# ggsave("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/glioma_TCGA_wholearray_EPIC450k_tSNE_out.pdf",width = 6, height = 4)


library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne[, ], aes(x=V1, V2, colour=IDH.status)) +
  geom_point() +
  scale_color_manual(values=c("orange","green", "purple", "red","darkgreen","darkblue"), name="Legend") +
  # scale_color_viridis(discrete=TRUE) +
  
  ggtitle(label = "t-SNE Genome wide (380k probes)",
          subtitle = "perplexity=30")






# s3s3s3s3s3s3s3s3s3ss3s3s3s3s
# s3s3s3s3s3s3s3s3s3ss3s3s3s3s
# s3s3s3s3s3s3s3s3s3ss3s3s3s3s
# Plot PCA each have the MostVariableFeature fine tuned ---------

# PCA 
pca_data <- DNAmtx_tProbes
dim(pca_data) # 371 641
table(is.na(pca_data))
# FALSE      
# 237811

pca2 <- prcomp(t(pca_data)) # this is the part where you calculate the PCA
aux <- as.data.frame(pca2$x[, 1:3]) #get only the info you need
pca_metadata <- merge(metadata, aux, by.y=0, by.x="Case", all.x=T)
#saveRDS(pca_metadata, file = '/.rds')


library(ggplot2); theme_set(theme_classic())
ggplot(pca_metadata, aes(x=PC1, y=PC2, colour=Ceccarelli_six_subtypes)) +
  geom_point() +
  scale_color_manual(values=c("orange","green", "purple", "red","darkgreen","darkblue"), name="Legend") +
  #scale_fill_manual(values=c(pal), name="Legend") + 
  #scale_color_manual(values=c(pal), name="Legend") +
  #scale_color_viridis(discrete=TRUE) +
  xlab(paste0("PC1 (",prettyNum(summary(pca2)$importance[2,1]*100, digits = 2),"%)")) +
  ylab(paste0("PC2 (",prettyNum(summary(pca2)$importance[2,2]*100, digits = 2 ),"%)")) +
  ggtitle("Gliomas TCGA - EPIC/450k 371 fine tuned probes") +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
  ) 

# ggsave(file='/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/glioma_TCGA_wholearray_EPIC450k_371finetunedprobes_PCA_out.pdf',width = 12, height = 6, dpi = 300,
# ) 


# t-SNE 
library(Rtsne)
tsne_realData <- Rtsne(t(pca_data), perplexity=30, check_duplicates = FALSE) # #function to run t-sne
### Plotando o tnse 
# add as 2 colunas (V1 e V2) com os dados do tsne
pdata.teste.tsne <- metadata #pData inicial
pdata.teste.tsne$V1 <- tsne_realData$Y[,1]
pdata.teste.tsne$V2 <- tsne_realData$Y[,2]

library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne[, ], aes(x=V1, V2, colour=Ceccarelli_six_subtypes)) +
  geom_point() +
  scale_color_manual(values=c("orange","green", "purple", "red","darkgreen","darkblue"), name="Legend") +
  # scale_color_viridis(discrete=TRUE) +
  
  ggtitle(label = "t-SNE Gliomas TCGA - EPIC/450k 371 fine tuned probes",
          subtitle = "perplexity=30")

# ggsave("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/glioma_TCGA_wholearray_EPIC450k_371finetunedprobes_tSNE_out.pdf",width = 6, height = 4)




