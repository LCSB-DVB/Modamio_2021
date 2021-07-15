## JM adaptation from EG 

## gene level analysis of transcriptomic data starting from foldchage values 
library(dplyr)
library(ggplot2) 


## set path and load file 

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis logFoldChange")
getwd() 

library(readxl)
midbrain_sncatrip_wt_edger_ranking_genelevel_with_logfc <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis logFoldChange/snca_mutant_vs_wildtype_2021_rsubread_edger_ranking_complete.xls")
View(snca_mutant_vs_wildtype_2021_rsubread_edger_ranking_complete.xls)

data <- midbrain_sncatrip_wt_edger_ranking_genelevel_with_logfc
head(data)

## filter significant genes based on pvalue and FDR value 
library(dplyr)

# select pvalu3 significant  < 0.05
length(which(data$PValue < 0.05)) # number of gene with Pvalue lower than 0.05 
data_Pvaluelist = rownames(data)[which(data$PValue < 0.05)]# list genes  
data_PvalueFrame  <- data %>% filter(data$PValue < 0.05)
write.csv(data_PvalueFrame, file = 'Pvalue_significant.csv')


# select FDR value significant < 0.05
length(which(data$FDR < 0.05)) # number of gene with FDR lower than 0.05 
data_FDRlist = rownames(data)[which(data$FDR < 0.05)]# list genes  
data_FRDFrame <- data %>% filter(data$FDR < 0.05)
write.csv(data_FRDFrame, file = 'FDR_significant.csv')

################################################################################################
################################ neuroprotective thropic genes (E.glaab)

## load protective/neurothrophic genes 
library(readxl)
neurothrophic <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data/neuropro_v1.xls")
View(neurothrophic)
head(neurothrophic)

# cross with FDR significant
### obtain data again 

neuroprotective_genes <- data.frame(neurothrophic$`Gene Symbol`)
colnames(neuroprotective_genes) <- c("`Gene symbol`")
View(neuroprotective_genes)
head(neuroprotective_genes)
library(dplyr) 
library(ggplot2) 
## FDR 

target <- c(intersect(neurothrophic$`Gene Symbol`, data$`Gene symbol`))
neuroprotective_genes_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(neurothrophic$`Gene Symbol`, data_PvalueFrame$`Gene symbol`))
neuroprotective_genes_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(neurothrophic$`Gene Symbol`, data_FRDFrame$`Gene symbol`))
neuroprotective_genes_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis logFoldChange")
write.csv(neuroprotective_genes_TotalData, file = 'neuroprotective_genes_TotalData.csv')
write.csv(neuroprotective_genes_Pvalue, file = 'neuroprotective_genes_Pvalue.csv')
write.csv(neuroprotective_genes_FDR, file = 'neuroprotective_genes_FDR.csv')

# plotting  FDR
library(ggplot2) 

# add colour & position labels & levels
neuroprotective_genes_FDR2 <-neuroprotective_genes_FDR[with(neuroprotective_genes_FDR, order(neuroprotective_genes_FDR$logFC)),] 
neuroprotective_genes_FDR2$order <- 1:nrow(neuroprotective_genes_FDR2)
neuroprotective_genes_FDR2[,c("order",setdiff(names(neuroprotective_genes_FDR2),"order"))]


neuroprotective_genes_FDR2$colour <- ifelse(neuroprotective_genes_FDR2$logFC < 0, "firebrick1","steelblue")
neuroprotective_genes_FDR2$Gene_Expression <- ifelse(neuroprotective_genes_FDR2$logFC < 0, "down-regulated","up-regulated")
neuroprotective_genes_FDR2$hjust <- ifelse(neuroprotective_genes_FDR2$logFC > 0, 1.3, -0.3)

 ggplot(neuroprotective_genes_FDR2, aes(order,logFC, label = `Gene symbol`,
                                             hjust = hjust)) + geom_text(aes(y = 0,
                                                                             colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with FDR < 0.05", 
       title= "neuroprotective/thropic genes") +
  coord_flip()+
  theme_classic()

 

 ## load libraries
 library(tidyverse)
 
 neuroprotective_genes_FDR %>%
   mutate(upregulated = logFC > 0) -> neuroprotective_genes_FDR
 
 ggplot(data = neuroprotective_genes_FDR,
        aes(x = reorder(`Gene symbol`, logFC), y = logFC,
            fill = upregulated))+
   geom_bar(stat = "identity")+
   coord_flip()+
   labs(x = "Genes FDR < 0.05", y = "LogFC",
        title = "Neuroprotective/Neurothropic genes",
        subtitles = "Fold change value for significant genes with FDR < 0.05")+
   scale_fill_manual(values=c('red4','navyblue'))+
   #theme_minimal()+
   guides(fill = FALSE)
 
## 
# plotting  Pvalue
library(ggplot2) 

# add colour & position labels & levels
neuroprotective_genes_Pvalue2 <-neuroprotective_genes_Pvalue[with(neuroprotective_genes_Pvalue, order(neuroprotective_genes_Pvalue$logFC)),] 
neuroprotective_genes_Pvalue2$order <- 1:nrow(neuroprotective_genes_Pvalue2)
neuroprotective_genes_Pvalue2[,c("order",setdiff(names(neuroprotective_genes_Pvalue2),"order"))]


neuroprotective_genes_Pvalue2$colour <- ifelse(neuroprotective_genes_Pvalue2$logFC < 0, "firebrick1","steelblue")
neuroprotective_genes_Pvalue2$Gene_Expression <- ifelse(neuroprotective_genes_Pvalue2$logFC < 0, "down-regulated","up-regulated")
neuroprotective_genes_Pvalue2$hjust <- ifelse(neuroprotective_genes_Pvalue2$logFC > 0, 1.3, -0.3)

ggplot(neuroprotective_genes_Pvalue2, aes(order,logFC, label = `Gene symbol`,
                                       hjust = hjust)) + geom_text(aes(y = 0,
                                                                       colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with Pvalue < 0.05", 
       title= "neuroprotective/thropic genes") +
  coord_flip()+
  theme_classic()

################################################################################################
################################################################################################
# Data from post-mortem human samples 

PMBT_integration_results <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis logFoldChange/crossed info with pathways/PMBT_integration_results.xlsx")
View(PMBT_integration_results)

# filter by P value

# select pvalu3 significant  < 0.05
length(which(PMBT_integration_results$`Adjusted P value` < 0.05)) # number of gene with Pvalue lower than 0.05 
postmortem_Pvaluelist = rownames(PMBT_integration_results)[which(PMBT_integration_results$`Adjusted P value` < 0.05)]# list genes  
postmortem_Pvalueframe  <- PMBT_integration_results %>% filter(PMBT_integration_results$`Adjusted P value` < 0.05)


# 153 common genes deregulated in postmortem samples PD (Pvalue <0,05) and hMos (FDR <0.05)
target <- c(intersect(postmortem_Pvalueframe$SYMBOL, data_FRDFrame$`Gene symbol`))
postmortemPvalue_and_hMO <- filter(postmortem_Pvalueframe, SYMBOL %in% target)
hmoFDR_and_postmortem <- filter(data_FRDFrame, `Gene symbol` %in% target)
rm(target)

## load libraries
library(tidyverse)

hmoFDR_and_postmortem %>%
  mutate(upregulated = logFC > 0) -> hmoFDR_and_postmortem

ggplot(data = hmoFDR_and_postmortem,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Deregulated gene in hMOs and PD postmortem samples",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)


write.csv(postmortemPvalue_and_hMO, file = 'postmortem_hMO.csv')
write.csv(hmoFDR_and_postmortem, file = 'hmo_postmortem.csv')

################################################################################################
################################################################################################
# Data from astrocytes - human 
library(readxl)
geneontology_Astrocyte_homoSapiens <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis logFoldChange\crossed info with pathways/geneontology_Astrocyte_homoSapiens.xlsx")
View(geneontology_Astrocyte_homoSapiens)


# 153 common genes deregulated in postmortem samples PD (Pvalue <0,05) and hMos (FDR <0.05)
target <- c(intersect(geneontology_Astrocyte_homoSapiens$Gene_symbol, data_FRDFrame$`Gene symbol`))
astrocyte_genes <- filter(data_FRDFrame, `Gene symbol` %in% target)
rm(target)

# add colour & position labels & levels
astrocyte_genes2 <-astrocyte_genes[with(astrocyte_genes, order(astrocyte_genes$logFC)),] 
astrocyte_genes2$order <- 1:nrow(astrocyte_genes2)
astrocyte_genes2[,c("order",setdiff(names(astrocyte_genes2),"order"))]


astrocyte_genes2$colour <- ifelse(astrocyte_genes2$logFC < 0, "darkred","darkblue")
astrocyte_genes2$Gene_Expression <- ifelse(astrocyte_genes2$logFC < 0, "down-regulated","up-regulated")
astrocyte_genes2$hjust <- ifelse(astrocyte_genes2$logFC > 0, 1.3, -0.3)

ggplot(astrocyte_genes2, aes(order,logFC, label = `Gene symbol`,
                            hjust = hjust)) + geom_text(aes(y = 0,
                                                            colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with FDR < 0.05", 
       title= "Astrocyte asociated genes") +
  coord_flip()+
  theme_classic()

astrocyte_genes %>%
  mutate(upregulated = logFC > 0) -> astrocyte_genes

ggplot(data = astrocyte_genes,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Astocyte associated genes",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)


################################################################################################
################################ pathcards data 

## load genes PD related 

PDgenes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/PD asociated genes eDGAR database simplified.xlsx")
View(PDgenes)
head(PDgenes)

# obtain cross PD genes with total data 
PD_data <- data.frame(intersect(PDgenes$`Gene name`, data$`Gene symbol`))
PD_data_sum <- data.frame(matrix(ncol = 6, nrow = 27))
PD_data_sum[] = NA
colnames(PD_data_sum) <- colnames(data) 

for (i in 1:length(PD_data$intersect.PDgenes..Gene.name...data..Gene.symbol..)){
  PD_data_sum[i,] = filter(data, `Gene symbol` == PD_data$intersect.PDgenes..Gene.name...data..Gene.symbol..[i])
}

# obtain cross PD genes with pValue significant data 

PD_dataPvalue <- data.frame(intersect(PDgenes$`Gene name`, data_PvalueFrame$`Gene symbol`))
PD_dataPvalue_sum <- data.frame(matrix(ncol = 6, nrow = 9))
PD_dataPvalue_sum[] = NA
colnames(PD_dataPvalue_sum) <- colnames(data) 

for (i in 1:length(PD_dataPvalue$intersect.PDgenes..Gene.name...data_PvalueFrame..Gene.symbol..)){
  PD_dataPvalue_sum[i,] = filter(data, `Gene symbol` == PD_dataPvalue$intersect.PDgenes..Gene.name...data_PvalueFrame..Gene.symbol..[i])
}

# obtain cross PD genes with FDR significant data 

PD_dataFDR <- data.frame(intersect(PDgenes$`Gene name`, data_FRDFrame$`Gene symbol`))
PD_dataFDR_sum = filter(data_FRDFrame,`Gene symbol` == PD_dataFDR$intersect.PDgenes..Gene.name...data_FRDFrame..Gene.symbol..)

for (i in 1:length(PD_dataFDR$intersect.PDgenes..Gene.name...data_FRDFrame..Gene.symbol..)){
  PD_dataFDR_sum[i,] = filter(data, `Gene symbol` == PD_dataFDR$intersect.PDgenes..Gene.name...data_FRDFrame..Gene.symbol..[i])
}

# saving 

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(PD_data_sum, file = 'PD_data_sum.csv')
write.csv(PD_dataPvalue_sum, file = 'PD_dataPvalue_sum.csv')
write.csv(PD_dataFDR_sum, file = 'PD_dataFDR_sum.csv')


# plotting  FDR
library(ggplot2) 

# add colour & position labels & levels
PD_dataFDR_sum2 <-PD_dataFDR_sum[with(PD_dataFDR_sum, order(PD_dataFDR_sum$logFC)),] 
PD_dataFDR_sum2$order <- 1:nrow(PD_dataFDR_sum2)
PD_dataFDR_sum2[,c("order",setdiff(names(PD_dataFDR_sum2),"order"))]


PD_dataFDR_sum2$colour <- ifelse(PD_dataFDR_sum2$logFC < 0, "firebrick1","steelblue")
PD_dataFDR_sum2$Gene_Expression <- ifelse(PD_dataFDR_sum2$logFC < 0, "down-regulated","up-regulated")
PD_dataFDR_sum2$hjust <- ifelse(PD_dataFDR_sum2$logFC > 0, 1.3, -0.3)

ggplot(PD_dataFDR_sum2, aes(order,logFC, label = `Gene symbol`,
                                          hjust = hjust)) + geom_text(aes(y = 0,
                                                                          colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with Pvalue < 0.05", 
       title= "PD asociated genes") +
  coord_flip()+
  theme_classic()

PD_dataFDR_sum %>%
  mutate(upregulated = logFC > 0) -> PD_dataFDR_sum

ggplot(data = PD_dataFDR_sum,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "PD associated genes",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)



##########################################
## load genes dopaminergic diferentiation 

library(readxl)
pathcards_Dopaminergic_neurogenesis_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards Dopaminergic neurogenesis genes.xlsx")
View(pathcards_Dopaminergic_neurogenesis_genes)

dopaminergic_diferentiation <- pathcards_Dopaminergic_neurogenesis_genes
View(dopaminergic_diferentiation)



target <- c(intersect(dopaminergic_diferentiation$Gene_symbol, data$`Gene symbol`))
dopaminergic_diferentiation_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(dopaminergic_diferentiation$Gene_symbol, data_PvalueFrame$`Gene symbol`))
dopaminergic_diferentiation_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(dopaminergic_diferentiation$Gene_symbol, data_FRDFrame$`Gene symbol`))
dopaminergic_diferentiation_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(dopaminergic_diferentiation_TotalData, file = 'dopaminergic_diferentiation_TotalData.csv')
write.csv(dopaminergic_diferentiation_Pvalue, file = 'dopaminergic_diferentiation_Pvalue.csv')
write.csv(dopaminergic_diferentiation_FDR, file = 'dopaminergic_diferentiation_FDR.csv')

# plotting 
library(ggplot2) 

# add colour & position labels & levels
dopaminergic_diferentiation_FDR2 <-dopaminergic_diferentiation_FDR[with(dopaminergic_diferentiation_FDR, order(dopaminergic_diferentiation_FDR$logFC)),] 
dopaminergic_diferentiation_FDR2$order <- 1:nrow(dopaminergic_diferentiation_FDR2)
dopaminergic_diferentiation_FDR2[,c("order",setdiff(names(dopaminergic_diferentiation_FDR2),"order"))]


dopaminergic_diferentiation_FDR2$colour <- ifelse(dopaminergic_diferentiation_FDR2$logFC < 0, "firebrick1","steelblue")
dopaminergic_diferentiation_FDR2$Gene_Expression <- ifelse(dopaminergic_diferentiation_FDR2$logFC < 0, "down-regulated","up-regulated")
dopaminergic_diferentiation_FDR2$hjust <- ifelse(dopaminergic_diferentiation_FDR2$logFC > 0, 1.3, -0.3)

ggplot(dopaminergic_diferentiation_FDR2, aes(order,logFC, label = `Gene symbol`,
                                          hjust = hjust)) + geom_text(aes(y = 0,
                                                                          colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with FDR < 0.05", 
       title= "Dopaminergic Diferentiation") +
  coord_flip()+
  theme_classic()


dopaminergic_diferentiation_FDR %>%
  mutate(upregulated = logFC > 0) -> dopaminergic_diferentiation_FDR

ggplot(data = dopaminergic_diferentiation_FDR,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Genes associated to Dopaminergic Diferentiation",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)



##########################################
## load genes Tyrosine metabolism

pathcards_Tyrosine_metabolism_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards Tyrosine metabolism genes.xlsx")
 View(pathcards_Tyrosine_metabolism_genes)


Tyrosine_metabolism <- pathcards_Tyrosine_metabolism_genes
View(Tyrosine_metabolism)

target <- c(intersect(Tyrosine_metabolism$Gene_symbol, data$`Gene symbol`))
Tyrosine_metabolism_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Tyrosine_metabolism$Gene_symbol, data_PvalueFrame$`Gene symbol`))
Tyrosine_metabolism_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Tyrosine_metabolism$Gene_symbol, data_FRDFrame$`Gene symbol`))
Tyrosine_metabolism_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(Tyrosine_metabolism_TotalData, file = 'Tyrosine_metabolism_TotalData.csv')
write.csv(Tyrosine_metabolism_Pvalue, file = 'Tyrosine_metabolism_Pvalue.csv')
write.csv(Tyrosine_metabolism_FDR, file = 'Tyrosine_metabolism_FDR.csv')

# add colour & position labels & levels
Tyrosine_metabolism_Pvalue2 <-Tyrosine_metabolism_Pvalue[with(Tyrosine_metabolism_Pvalue, order(Tyrosine_metabolism_Pvalue$logFC)),] 
Tyrosine_metabolism_Pvalue2$order <- 1:nrow(Tyrosine_metabolism_Pvalue2)
Tyrosine_metabolism_Pvalue2[,c("order",setdiff(names(Tyrosine_metabolism_Pvalue2),"order"))]


Tyrosine_metabolism_Pvalue2$colour <- ifelse(Tyrosine_metabolism_Pvalue2$logFC < 0, "firebrick1","steelblue")
Tyrosine_metabolism_Pvalue2$Gene_Expression <- ifelse(Tyrosine_metabolism_Pvalue2$logFC < 0, "down-regulated","up-regulated")
Tyrosine_metabolism_Pvalue2$hjust <- ifelse(Tyrosine_metabolism_Pvalue2$logFC > 0, 1.3, -0.3)

ggplot(Tyrosine_metabolism_Pvalue2, aes(order,logFC, label = `Gene symbol`,
                                             hjust = hjust)) + geom_text(aes(y = 0,
                                                                             colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with Pvalue < 0.05", 
       title= "Tyrosine metabolism") +
  coord_flip() +
  theme_classic()


Tyrosine_metabolism_FDR %>%
  mutate(upregulated = logFC > 0) -> Tyrosine_metabolism_FDR

ggplot(data = Tyrosine_metabolism_FDR,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Genes associated to Tyrosine metabolism",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)



Tyrosine_metabolism_Pvalue2
ggplot(data = Tyrosine_metabolism_Pvalue2,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = Gene_Expression))+
  geom_bar(stat = "identity")+
  #coord_flip()+
  labs(x = "Genes Pvalue < 0.05", y = "LogFC",
       title = "Tyrosine metabolism",
       subtitles = "Fold change value for significant genes with Pvalue < 0.05")+
  scale_fill_manual(values=c('gray57','goldenrod2'))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  
  theme_classic() +
  theme(
    axis.line = element_line(colour = 'grey', size = 0.1) ,
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=9, angle = 90),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size=9),  #summaryPlot$colourPath
    axis.ticks.y = element_line(),
    legend.position="right",
    legend.text = element_text(size = 11, face = "bold") ,
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.key.width = unit(0.6,"cm"),
    plot.title = element_text(size = 15, hjust=0.5, vjust= 1, face = "bold"),
    plot.subtitle = element_blank(),#element_text(size = 2, hjust=0.5)
    strip.text = element_text(size=12, vjust=0.5),
    strip.background = element_rect(fill="lightgray"),
    # panel.border = element_rect(fill = NA, color = "black"),
    panel.spacing.y = unit(0.8, "lines"),
    strip.switch.pad.wrap=unit(20, "lines")
  ) 

##########################################
## load genes Parkinsons Disease Pathway

pathcards_Parkinson_s_disease_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards Parkinson's disease  genes.xlsx")
View(pathcards_Parkinson_s_disease_genes)

Parkinsons_Disease_Pathway <- pathcards_Parkinson_s_disease_genes
View(Parkinsons_Disease_Pathway)

target <- c(intersect(Parkinsons_Disease_Pathway$Gene_symbol, data$`Gene symbol`))
Parkinsons_Disease_Pathway_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Parkinsons_Disease_Pathway$Gene_symbol, data_PvalueFrame$`Gene symbol`))
Parkinsons_Disease_Pathway_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Parkinsons_Disease_Pathway$Gene_symbol, data_FRDFrame$`Gene symbol`))
Parkinsons_Disease_Pathway_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(Parkinsons_Disease_Pathway_TotalData, file = 'Parkinsons_Disease_Pathway_TotalData.csv')
write.csv(Parkinsons_Disease_Pathway_Pvalue, file = 'Parkinsons_Disease_Pathway_Pvalue.csv')
write.csv(Parkinsons_Disease_Pathway_FDR, file = 'Parkinsons_Disease_Pathway_FDR.csv')

# add colour & position labels & levels
Parkinsons_Disease_Pathway_Pvalue2 <-Parkinsons_Disease_Pathway_Pvalue[with(Parkinsons_Disease_Pathway_Pvalue, order(Parkinsons_Disease_Pathway_Pvalue$logFC)),] 
Parkinsons_Disease_Pathway_Pvalue2$order <- 1:nrow(Parkinsons_Disease_Pathway_Pvalue2)
Parkinsons_Disease_Pathway_Pvalue2[,c("order",setdiff(names(Parkinsons_Disease_Pathway_Pvalue2),"order"))]


Parkinsons_Disease_Pathway_Pvalue2$colour <- ifelse(Parkinsons_Disease_Pathway_Pvalue2$logFC < 0, "firebrick1","steelblue")
Parkinsons_Disease_Pathway_Pvalue2$Gene_Expression <- ifelse(Parkinsons_Disease_Pathway_Pvalue2$logFC < 0, "down-regulated","up-regulated")
Parkinsons_Disease_Pathway_Pvalue2$hjust <- ifelse(Parkinsons_Disease_Pathway_Pvalue2$logFC > 0, 1.3, -0.3)

ggplot(Parkinsons_Disease_Pathway_Pvalue2, aes(order,logFC, label = `Gene symbol`,
                                             hjust = hjust)) + geom_text(aes(y = 0,
                                                                             colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with Pvalue < 0.05", 
       title= "Parkinson's disease associated genes") +
  coord_flip() +
  theme_classic()



Parkinsons_Disease_Pathway_FDR %>%
  mutate(upregulated = logFC > 0) -> Parkinsons_Disease_Pathway_FDR

ggplot(data = Parkinsons_Disease_Pathway_FDR,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Genes associated to Parkinson's Disease",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)



##########################################
## load genes Dopamine metabolism

pathcards_Dopamine_metabolism_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards Dopamine metabolism genes.xlsx")
View(pathcards_Dopamine_metabolism_genes)

Dopamine_metabolism <- pathcards_Dopamine_metabolism_genes
View(Dopamine_metabolism)

target <- c(intersect(Dopamine_metabolism$Gene_symbol, data$`Gene symbol`))
Dopamine_metabolism_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Dopamine_metabolism$Gene_symbol, data_PvalueFrame$`Gene symbol`))
Dopamine_metabolism_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Dopamine_metabolism$Gene_symbol, data_FRDFrame$`Gene symbol`))
Dopamine_metabolism_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(Dopamine_metabolism_TotalData, file = 'Dopamine_metabolism_TotalData.csv')
write.csv(Dopamine_metabolism_Pvalue, file = 'Dopamine_metabolism_Pvalue.csv')
write.csv(Dopamine_metabolism_FDR, file = 'Dopamine_metabolism_FDR.csv')


# add colour & position labels & levels
Dopamine_metabolism_TotalData2 <-Dopamine_metabolism_TotalData[with(Dopamine_metabolism_TotalData, order(Dopamine_metabolism_TotalData$logFC)),] 
Dopamine_metabolism_TotalData2$order <- 1:nrow(Dopamine_metabolism_TotalData2)
Dopamine_metabolism_TotalData2 <- Dopamine_metabolism_TotalData2[,c("order",setdiff(names(Dopamine_metabolism_TotalData2),"order"))]


Dopamine_metabolism_TotalData2$colour <- ifelse(Dopamine_metabolism_TotalData2$logFC < 0, "firebrick1","steelblue")
Dopamine_metabolism_TotalData2$Gene_Expression <- ifelse(Dopamine_metabolism_TotalData2$logFC < 0, "down-regulated","up-regulated")
Dopamine_metabolism_TotalData2$hjust <- ifelse(Dopamine_metabolism_TotalData2$logFC > 0, 1.3, -0.3)

ggplot(Dopamine_metabolism_TotalData2, aes(order,logFC, label = `Gene symbol`,
                                               hjust = hjust)) + geom_text(aes(y = 0,
                                                                               colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for all genes (MAOB Pvalue < 0.05)", 
       title= "Dopamine metabolism associated genes") +
  coord_flip() +
  theme_classic()



Dopamine_metabolism_Pvalue %>%
  mutate(upregulated = logFC > 0) -> Dopamine_metabolism_Pvalue

ggplot(data = Dopamine_metabolism_Pvalue,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes Pvalue < 0.05", y = "LogFC",
       title = "Genes associated to Dopamine metabolism",
       subtitles = "Fold change value for significant genes with Pvalue < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)



Dopamine_metabolism_TotalData2
ggplot(data = Dopamine_metabolism_TotalData2,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = Gene_Expression))+
  geom_bar(stat = "identity")+
  #coord_flip()+
  labs(x = "Genes Pvalue < 0.05", y = "LogFC",
       title = "Dopamine metabolism",
       subtitles = "Fold change value for significant genes with Pvalue < 0.05")+
  scale_fill_manual(values=c('gray57','goldenrod2'))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  
  theme_classic() +
  theme(
    axis.line = element_line(colour = 'grey', size = 0.1) ,
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=9, angle = 90),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size=9),  #summaryPlot$colourPath
    axis.ticks.y = element_line(),
    legend.position="right",
    legend.text = element_text(size = 11, face = "bold") ,
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.key.width = unit(0.6,"cm"),
    plot.title = element_text(size = 15, hjust=0.5, vjust= 1, face = "bold"),
    plot.subtitle = element_blank(),#element_text(size = 2, hjust=0.5)
    strip.text = element_text(size=12, vjust=0.5),
    strip.background = element_rect(fill="lightgray"),
    # panel.border = element_rect(fill = NA, color = "black"),
    panel.spacing.y = unit(0.8, "lines"),
    strip.switch.pad.wrap=unit(20, "lines")
  ) 



##########################################
## load genes axon guiadance 

pathcards_Axonal_guiadance_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards Axonal guiadance genes.xlsx")
View(pathcards_Axonal_guiadance_genes)

axon_guiadance  <- pathcards_Axonal_guiadance_genes
View(axon_guiadance)

target <- c(intersect(axon_guiadance$Gene_symbol, data$`Gene symbol`))
axon_guiadance_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(axon_guiadance$Gene_symbol, data_PvalueFrame$`Gene symbol`))
axon_guiadance_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(axon_guiadance$Gene_symbol, data_FRDFrame$`Gene symbol`))
axon_guiadance_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(axon_guiadance_TotalData, file = 'axon_guiadance_TotalData.csv')
write.csv(axon_guiadance_Pvalue, file = 'axon_guiadance_Pvalue.csv')
write.csv(axon_guiadance_FDR, file = 'axon_guiadance_FDR.csv')


# add colour & position labels & levels
axon_guiadance_FDR2 <-axon_guiadance_FDR[with(axon_guiadance_FDR, order(axon_guiadance_FDR$logFC)),] 
axon_guiadance_FDR2$order <- 1:nrow(axon_guiadance_FDR2)
axon_guiadance_FDR2[,c("order",setdiff(names(axon_guiadance_FDR2),"order"))]


axon_guiadance_FDR2$colour <- ifelse(axon_guiadance_FDR2$logFC < 0, "firebrick1","steelblue")
axon_guiadance_FDR2$Gene_Expression <- ifelse(axon_guiadance_FDR2$logFC < 0, "down-regulated","up-regulated")
axon_guiadance_FDR2$hjust <- ifelse(axon_guiadance_FDR2$logFC > 0, 1.3, -0.3)

ggplot(axon_guiadance_FDR2, aes(order,logFC, label = `Gene symbol`,
                                             hjust = hjust)) + geom_text(aes(y = 0,
                                                                             colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with FDR < 0.05", 
       title= "Axonal projection/guiadance") +
  coord_flip()+
  theme_classic()



axon_guiadance_FDR %>%
  mutate(upregulated = logFC > 0) -> axon_guiadance_FDR

ggplot(data = axon_guiadance_FDR,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes Pvalue < 0.05", y = "LogFC",
       title = "Genes associated to axon guiadance",
       subtitles = "Fold change value for significant genes with Pvalue < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)


##########################################
## Neural Stem Cell 

pathcards_Neuronal_stemm_cells_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards Neuronal stemm cells genes.xlsx")
View(pathcards_Neuronal_stemm_cells_genes)

Neural_Stem_Cell_diff  <- pathcards_Neuronal_stemm_cells_genes
View(Neural_Stem_Cell_diff)

target <- c(intersect(Neural_Stem_Cell_diff$Gene_symbol, data$`Gene symbol`))
Neural_Stem_Cell_diff_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Neural_Stem_Cell_diff$Gene_symbol, data_PvalueFrame$`Gene symbol`))
Neural_Stem_Cell_diff_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Neural_Stem_Cell_diff$Gene_symbol, data_FRDFrame$`Gene symbol`))
Neural_Stem_Cell_diff_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(Neural_Stem_Cell_diff_TotalData, file = 'Neural_Stem_Cell_diff_TotalData.csv')
write.csv(Neural_Stem_Cell_diff_Pvalue, file = 'Neural_Stem_Cell_diff_Pvalue.csv')
write.csv(Neural_Stem_Cell_diff_FDR, file = 'Neural_Stem_Cell_diff_FDR.csv')

# add colour & position labels & levels
Neural_Stem_Cell_diff_FDR2 <-Neural_Stem_Cell_diff_FDR[with(Neural_Stem_Cell_diff_FDR, order(Neural_Stem_Cell_diff_FDR$logFC)),] 
Neural_Stem_Cell_diff_FDR2$order <- 1:nrow(Neural_Stem_Cell_diff_FDR2)
Neural_Stem_Cell_diff_FDR2[,c("order",setdiff(names(Neural_Stem_Cell_diff_FDR2),"order"))]


Neural_Stem_Cell_diff_FDR2$colour <- ifelse(Neural_Stem_Cell_diff_FDR2$logFC < 0, "firebrick1","steelblue")
Neural_Stem_Cell_diff_FDR2$Gene_Expression <- ifelse(Neural_Stem_Cell_diff_FDR2$logFC < 0, "down-regulated","up-regulated")
Neural_Stem_Cell_diff_FDR2$hjust <- ifelse(Neural_Stem_Cell_diff_FDR2$logFC > 0, 1.3, -0.3)

ggplot(Neural_Stem_Cell_diff_FDR2, aes(order,logFC, label = `Gene symbol`,
                                hjust = hjust)) + geom_text(aes(y = 0,
                                                                colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with FDR < 0.05", 
       title= "Neuronal stemm cell associated genes") +
  coord_flip()+
  ylim(-5, 4)+
  theme_classic()


Neural_Stem_Cell_diff_FDR %>%
  mutate(upregulated = logFC > 0) -> Neural_Stem_Cell_diff_FDR

ggplot(data = Neural_Stem_Cell_diff_FDR,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes Pvalue < 0.05", y = "LogFC",
       title = "Genes associated to Neural stemm cell differentiation",
       subtitles = "Fold change value for significant genes with Pvalue < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)


##########################################
## Synaptic vesicle cycle

pathcards_Synaptic_vesicle_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards Synaptic vesicle genes.xlsx")
 View(pathcards_Synaptic_vesicle_genes)

Synaptic_vesicle_cycle <- pathcards_Synaptic_vesicle_genes
View(Synaptic_vesicle_cycle)

target <- c(intersect(Synaptic_vesicle_cycle$Gene_symbol, data$`Gene symbol`))
Synaptic_vesicle_cycle_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Synaptic_vesicle_cycle$Gene_symbol, data_PvalueFrame$`Gene symbol`))
Synaptic_vesicle_cycle_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Synaptic_vesicle_cycle$Gene_symbol, data_FRDFrame$`Gene symbol`))
Synaptic_vesicle_cycle_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(Synaptic_vesicle_cycle_TotalData, file = 'Synaptic_vesicle_cycle_TotalData.csv')
write.csv(Synaptic_vesicle_cycle_Pvalue, file = 'Synaptic_vesicle_cycle_Pvalue.csv')
write.csv(Synaptic_vesicle_cycle_FDR, file = 'Synaptic_vesicle_cycle_FDR.csv')

# add colour & position labels & levels
Synaptic_vesicle_cycle_FDR2 <-Synaptic_vesicle_cycle_FDR[with(Synaptic_vesicle_cycle_FDR, order(Synaptic_vesicle_cycle_FDR$logFC)),] 
Synaptic_vesicle_cycle_FDR2$order <- 1:nrow(Synaptic_vesicle_cycle_FDR2)
Synaptic_vesicle_cycle_FDR2[,c("order",setdiff(names(Synaptic_vesicle_cycle_FDR2),"order"))]


Synaptic_vesicle_cycle_FDR2$colour <- ifelse(Synaptic_vesicle_cycle_FDR2$logFC < 0, "firebrick1","steelblue")
Synaptic_vesicle_cycle_FDR2$Gene_Expression <- ifelse(Synaptic_vesicle_cycle_FDR2$logFC < 0, "down-regulated","up-regulated")
Synaptic_vesicle_cycle_FDR2$hjust <- ifelse(Synaptic_vesicle_cycle_FDR2$logFC > 0, 1.3, -0.3)

ggplot(Synaptic_vesicle_cycle_FDR2, aes(order,logFC, label = `Gene symbol`,
                                       hjust = hjust)) + geom_text(aes(y = 0,
                                                                       colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with FDR < 0.05", 
       title= "Synaptic vesicle associated genes") +
  coord_flip()+
  ylim(-2, 6)+
  theme_classic()

## 
ggplot(Synaptic_vesicle_cycle_FDR2, aes(order,logFC, label = `Gene symbol`,
       hjust = hjust)) + geom_text(aes(y = 0, colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with FDR < 0.05", 
       title= "Synaptic vesicle associated genes") +
  #coord_flip()+
  ylim(-2, 6)+
  #theme_classic()

## load libraries
library(tidyverse)

Synaptic_vesicle_cycle_FDR %>%
  mutate(upregulated = logFC > 0) -> Synaptic_vesicle_cycle_FDR

ggplot(data = Synaptic_vesicle_cycle_FDR,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Genes associated to synaptic vesicles ",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)


##########################################
## load genes Neurotransmitter Release Cycle 


pathcards_Neurotransmiter_release_cycle_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards Neurotransmiter release cycle genes.xlsx")
View(pathcards_Neurotransmiter_release_cycle_genes)

Neurotransmitter_Release_Cycle <- pathcards_Neurotransmiter_release_cycle_genes
View(Neurotransmitter_Release_Cycle)

target <- c(intersect(Neurotransmitter_Release_Cycle$Gene_symbol, data$`Gene symbol`))
Neurotransmitter_Release_Cycle_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Neurotransmitter_Release_Cycle$Gene_symbol, data_PvalueFrame$`Gene symbol`))
Neurotransmitter_Release_Cycle_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Neurotransmitter_Release_Cycle$Gene_symbol, data_FRDFrame$`Gene symbol`))
Neurotransmitter_Release_Cycle_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(Neurotransmitter_Release_Cycle_TotalData, file = 'Neurotransmitter_Release_Cycle_TotalData.csv')
write.csv(Neurotransmitter_Release_Cycle_Pvalue, file = 'Neurotransmitter_Release_Cycle_Pvalue.csv')
write.csv(Neurotransmitter_Release_Cycle_FDR, file = 'Neurotransmitter_Release_Cycle_FDR.csv')

# add colour & position labels & levels
Neurotransmitter_Release_Cycle_Pvalue2 <-Neurotransmitter_Release_Cycle_Pvalue[with(Neurotransmitter_Release_Cycle_Pvalue, order(Neurotransmitter_Release_Cycle_Pvalue$logFC)),] 
Neurotransmitter_Release_Cycle_Pvalue2$order <- 1:nrow(Neurotransmitter_Release_Cycle_Pvalue2)
Neurotransmitter_Release_Cycle_Pvalue2[,c("order",setdiff(names(Neurotransmitter_Release_Cycle_Pvalue2),"order"))]


Neurotransmitter_Release_Cycle_Pvalue2$colour <- ifelse(Neurotransmitter_Release_Cycle_Pvalue2$logFC < 0, "firebrick1","steelblue")
Neurotransmitter_Release_Cycle_Pvalue2$Gene_Expression <- ifelse(Neurotransmitter_Release_Cycle_Pvalue2$logFC < 0, "down-regulated","up-regulated")
Neurotransmitter_Release_Cycle_Pvalue2$hjust <- ifelse(Neurotransmitter_Release_Cycle_Pvalue2$logFC > 0, 1.3, -0.3)

ggplot(Neurotransmitter_Release_Cycle_Pvalue2, aes(order,logFC, label = `Gene symbol`,
                                        hjust = hjust)) + geom_text(aes(y = 0,
                                                                        colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with Pvalue < 0.05", 
       title= "Neurotransmitter Release Cycle associated genes") +
  coord_flip()+
  # ylim(-4, 2)+
  theme_classic()

## load libraries
library(tidyverse)

Neurotransmitter_Release_Cycle_FDR %>%
  mutate(upregulated = logFC > 0) -> Neurotransmitter_Release_Cycle_FDR

ggplot(data = Neurotransmitter_Release_Cycle_FDR,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Genes associated to Neurotransmiter release",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)


##########################################
## load genes Dopamine Neurotransmitter

pathcards_Dopamine_Neurotrasnsmitter_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards Dopamine Neurotrasnsmitter genes.xlsx")
View(pathcards_Dopamine_Neurotrasnsmitter_genes)

Dopamine_Neurotransmitter <- pathcards_Dopamine_Neurotrasnsmitter_genes
View(Dopamine_Neurotransmitter)

target <- c(intersect(Dopamine_Neurotransmitter$Gene_symbol, data$`Gene symbol`))
Dopamine_Neurotransmitter_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Dopamine_Neurotransmitter$Gene_symbol, data_PvalueFrame$`Gene symbol`))
Dopamine_Neurotransmitter_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Dopamine_Neurotransmitter$Gene_symbol, data_FRDFrame$`Gene symbol`))
Dopamine_Neurotransmitter_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(Dopamine_Neurotransmitter_TotalData, file = 'Dopamine_Neurotransmitter_TotalData.csv')
write.csv(Dopamine_Neurotransmitter_Pvalue, file = 'Dopamine_Neurotransmitter_Pvalue.csv')
write.csv(Dopamine_Neurotransmitter_FDR, file = 'Dopamine_Neurotransmitter_FDR.csv')

# add colour & position labels & levels
Dopamine_Neurotransmitter_Pvalue2 <-Dopamine_Neurotransmitter_Pvalue[with(Dopamine_Neurotransmitter_Pvalue, order(Dopamine_Neurotransmitter_Pvalue$logFC)),] 
Dopamine_Neurotransmitter_Pvalue2$order <- 1:nrow(Dopamine_Neurotransmitter_Pvalue2)
Dopamine_Neurotransmitter_Pvalue2[,c("order",setdiff(names(Dopamine_Neurotransmitter_Pvalue2),"order"))]


Dopamine_Neurotransmitter_Pvalue2$colour <- ifelse(Dopamine_Neurotransmitter_Pvalue2$logFC < 0, "firebrick1","steelblue")
Dopamine_Neurotransmitter_Pvalue2$Gene_Expression <- ifelse(Dopamine_Neurotransmitter_Pvalue2$logFC < 0, "down-regulated","up-regulated")
Dopamine_Neurotransmitter_Pvalue2$hjust <- ifelse(Dopamine_Neurotransmitter_Pvalue2$logFC > 0, 1.3, -0.3)

ggplot(Dopamine_Neurotransmitter_Pvalue2, aes(order,logFC, label = `Gene symbol`,
                                                   hjust = hjust)) + geom_text(aes(y = 0,
                                                                                   colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with Pvalue < 0.05", 
       title= "Dopamine Neurotransmitter associated genes") +
  coord_flip()+
  # ylim(-4, 2)+
  theme_classic()

Dopamine_Neurotransmitter_FDR %>%
  mutate(upregulated = logFC > 0) -> Dopamine_Neurotransmitter_FDR

ggplot(data = Dopamine_Neurotransmitter_FDR,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Genes associated to Dopamine neurotransmiter",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)

Dopamine_Neurotransmitter_Pvalue2

ggplot(data = Dopamine_Neurotransmitter_Pvalue2,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = Gene_Expression))+
  geom_bar(stat = "identity")+
  #coord_flip()+
  labs(x = "Genes Pvalue < 0.05", y = "LogFC",
       title = "Dopamine neurotransmiter",
       subtitles = "Fold change value for significant genes with Pvalue < 0.05")+
  scale_fill_manual(values=c('gray57','goldenrod2'))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  
  theme_classic() +
  theme(
    axis.line = element_line(colour = 'grey', size = 0.1) ,
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=9, angle = 90),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size=9),  #summaryPlot$colourPath
    axis.ticks.y = element_line(),
    legend.position="right",
    legend.text = element_text(size = 11, face = "bold") ,
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.key.width = unit(0.6,"cm"),
    plot.title = element_text(size = 15, hjust=0.5, vjust= 1, face = "bold"),
    plot.subtitle = element_blank(),#element_text(size = 2, hjust=0.5)
    strip.text = element_text(size=12, vjust=0.5),
    strip.background = element_rect(fill="lightgray"),
    # panel.border = element_rect(fill = NA, color = "black"),
    panel.spacing.y = unit(0.8, "lines"),
    strip.switch.pad.wrap=unit(20, "lines")
  ) 



##########################################
## load genes Alpha-synuclein

pathcards_Alpha_synuclein_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards Alpha synuclein genes.xlsx")
View(pathcards_Alpha_synuclein_genes)

Alpha_synuclein <- pathcards_Alpha_synuclein_genes
View(Alpha_synuclein)

target <- c(intersect(Alpha_synuclein$Gene_symbol, data$`Gene symbol`))
Alpha_synuclein_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Alpha_synuclein$Gene_symbol, data_PvalueFrame$`Gene symbol`))
Alpha_synuclein_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Alpha_synuclein$Gene_symbol, data_FRDFrame$`Gene symbol`))
Alpha_synuclein_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(Alpha_synuclein_TotalData, file = 'Alpha_synuclein_TotalData.csv')
write.csv(Alpha_synuclein_Pvalue, file = 'Alpha_synuclein_Pvalue.csv')
write.csv(Alpha_synuclein_FDR, file = 'Alpha_synuclein_FDR.csv')

# add colour & position labels & levels
Alpha_synuclein_Pvalue2 <-Alpha_synuclein_Pvalue[with(Alpha_synuclein_Pvalue, order(Alpha_synuclein_Pvalue$logFC)),] 
Alpha_synuclein_Pvalue2$order <- 1:nrow(Alpha_synuclein_Pvalue2)
Alpha_synuclein_Pvalue2[,c("order",setdiff(names(Alpha_synuclein_Pvalue2),"order"))]


Alpha_synuclein_Pvalue2$colour <- ifelse(Alpha_synuclein_Pvalue2$logFC < 0, "firebrick1","steelblue")
Alpha_synuclein_Pvalue2$Gene_Expression <- ifelse(Alpha_synuclein_Pvalue2$logFC < 0, "down-regulated","up-regulated")
Alpha_synuclein_Pvalue2$hjust <- ifelse(Alpha_synuclein_Pvalue2$logFC > 0, 1.3, -0.3)

ggplot(Alpha_synuclein_Pvalue2, aes(order,logFC, label = `Gene symbol`,
                                              hjust = hjust)) + geom_text(aes(y = 0,
                                                                              colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with Pvalue < 0.05", 
       title= "Alpha-synuclein associated genes") +
  coord_flip()+
  # ylim(-4, 2)+
  theme_classic()

Alpha_synuclein_FDR %>%
  mutate(upregulated = logFC > 0) -> Alpha_synuclein_FDR

ggplot(data = Alpha_synuclein_FDR,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Genes associated to Alpha-synuclein",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)


Alpha_synuclein_Pvalue %>%
  mutate(upregulated = logFC > 0) -> Alpha_synuclein_Pvalue

ggplot(data = Alpha_synuclein_Pvalue,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes Pvalue < 0.05", y = "LogFC",
       title = "Genes associated to Alpha-synuclein",
       subtitles = "Fold change value for significant genes with Pvalue < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)

##########################################
## load genes Autophagy Pathway

pathcards_Autophagy_pathway_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards Autophagy pathway genes.xlsx")
View(pathcards_Autophagy_pathway_genes)

Autophagy_Pathway<- pathcards_Autophagy_pathway_genes
View(Autophagy_Pathway)

target <- c(intersect(Autophagy_Pathway$Gene_symbol, data$`Gene symbol`))
Autophagy_Pathway_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Autophagy_Pathway$Gene_symbol, data_PvalueFrame$`Gene symbol`))
Autophagy_Pathway_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(Autophagy_Pathway$Gene_symbol, data_FRDFrame$`Gene symbol`))
Autophagy_Pathway_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(Autophagy_Pathway_TotalData, file = 'Autophagy_Pathway_TotalData.csv')
write.csv(Autophagy_Pathway_Pvalue, file = 'Autophagy_Pathway_Pvalue.csv')
write.csv(Autophagy_Pathway_FDR, file = 'Autophagy_Pathway_FDR.csv')

# add colour & position labels & levels
Autophagy_Pathway_Pvalue2 <-Autophagy_Pathway_Pvalue[with(Autophagy_Pathway_Pvalue, order(Autophagy_Pathway_Pvalue$logFC)),] 
Autophagy_Pathway_Pvalue2$order <- 1:nrow(Autophagy_Pathway_Pvalue2)
Autophagy_Pathway_Pvalue2[,c("order",setdiff(names(Autophagy_Pathway_Pvalue2),"order"))]


Autophagy_Pathway_Pvalue2$colour <- ifelse(Autophagy_Pathway_Pvalue2$logFC < 0, "firebrick1","steelblue")
Autophagy_Pathway_Pvalue2$Gene_Expression <- ifelse(Autophagy_Pathway_Pvalue2$logFC < 0, "down-regulated","up-regulated")
Autophagy_Pathway_Pvalue2$hjust <- ifelse(Autophagy_Pathway_Pvalue2$logFC > 0, 1.3, -0.3)

ggplot(Autophagy_Pathway_Pvalue2, aes(order,logFC, label = `Gene symbol`,
                                    hjust = hjust)) + geom_text(aes(y = 0,
                                                                    colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with Pvalue < 0.05", 
       title= "Autophagy pathway associated genes") +
  coord_flip()+
  ylim(-1, 3)+
  theme_classic()


##########################################
## load genes Autophagy Pathway

library(readxl)
GO_term_Neuronal_Maturation_20200407_024556 <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/GO_term_Neuronal_Maturation_20200407_024556.xlsx", 
                                                          sheet = "modif raw data")
View(GO_term_Neuronal_Maturation_20200407_024556)
Neuronal_maturation <- GO_term_Neuronal_Maturation_20200407_024556

target <- c(intersect(Neuronal_maturation$Gene_symbol, data_FRDFrame$`Gene symbol`))
Neuronal_maturation_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

# add colour & position labels & levels
Neuronal_maturation_FDR <-Neuronal_maturation_FDR[with(Neuronal_maturation_FDR, order(Neuronal_maturation_FDR$logFC)),] 
Neuronal_maturation_FDR$order <- 1:nrow(Neuronal_maturation_FDR)
Neuronal_maturation_FDR[,c("order",setdiff(names(Neuronal_maturation_FDR),"order"))]


Neuronal_maturation_FDR$colour <- ifelse(Neuronal_maturation_FDR$logFC < 0, "firebrick1","steelblue")
Neuronal_maturation_FDR$Gene_Expression <- ifelse(Neuronal_maturation_FDR$logFC < 0, "down-regulated","up-regulated")
Neuronal_maturation_FDR$hjust <- ifelse(Neuronal_maturation_FDR$logFC > 0, 1.3, -0.3)

ggplot(Neuronal_maturation_FDR, aes(order,logFC, label = `Gene symbol`,
                                      hjust = hjust)) + geom_text(aes(y = 0,
                                                                      colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with Pvalue < 0.05", 
       title= "Autophagy pathway associated genes") +
  coord_flip()+
  ylim(-1, 3)+
  theme_classic()

Neuronal_maturation_FDR %>%
  mutate(upregulated = logFC > 0) -> Neuronal_maturation_FDR

ggplot(data = Neuronal_maturation_FDR,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Genes associated to Synapse",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)




##########################################
## load genes protein-protein interactiosn at synapse


##########################################
## load genes protein-protein interactiosn at synapse

pathcards_PPI_at_synapse_genes <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/analysis logFC_GO 03042020/data obtained online/pathcards PPI at synapse genes.xlsx")
View(pathcards_PPI_at_synapse_genes)

PPI_at_synapse<- pathcards_PPI_at_synapse_genes
View(PPI_at_synapse)

target <- c(intersect(PPI_at_synapse$Gene_symbol, data$`Gene symbol`))
PPI_at_synapse_TotalData <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(PPI_at_synapse$Gene_symbol, data_PvalueFrame$`Gene symbol`))
PPI_at_synapse_Pvalue <- filter(data, `Gene symbol` %in% target)
rm(target)
target <- c(intersect(PPI_at_synapse$Gene_symbol, data_FRDFrame$`Gene symbol`))
PPI_at_synapse_FDR <- filter(data, `Gene symbol` %in% target)
rm(target)

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2020/Transcriptomic data/Pathway and network analysis based on gene fold change data")
write.csv(PPI_at_synapse_TotalData, file = 'PPI_at_synapse_TotalData.csv')
write.csv(PPI_at_synapse_Pvalue, file = 'PPI_at_synapse_Pvalue.csv')
write.csv(PPI_at_synapse_FDR, file = 'PPI_at_synapse_FDR.csv')

# add colour & position labels & levels
PPI_at_synapse_Pvalue2 <-PPI_at_synapse_Pvalue[with(PPI_at_synapse_Pvalue, order(PPI_at_synapse_Pvalue$logFC)),] 
PPI_at_synapse_Pvalue2$order <- 1:nrow(PPI_at_synapse_Pvalue2)
PPI_at_synapse_Pvalue2[,c("order",setdiff(names(PPI_at_synapse_Pvalue2),"order"))]


PPI_at_synapse_Pvalue2$colour <- ifelse(PPI_at_synapse_Pvalue2$logFC < 0, "firebrick1","steelblue")
PPI_at_synapse_Pvalue2$Gene_Expression <- ifelse(PPI_at_synapse_Pvalue2$logFC < 0, "down-regulated","up-regulated")
PPI_at_synapse_Pvalue2$hjust <- ifelse(PPI_at_synapse_Pvalue2$logFC > 0, 1.3, -0.3)

ggplot(PPI_at_synapse_Pvalue2, aes(order,logFC, label = `Gene symbol`,
                                      hjust = hjust)) + geom_text(aes(y = 0,
                                                                      colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with Pvalue < 0.05", 
       title= "protein-protein interaction at synapse associated genes") +
  coord_flip()+
  #ylim(-1, 2)+
  theme_classic()

PPI_at_synapse_FDR2 <-PPI_at_synapse_FDR[with(PPI_at_synapse_FDR, order(PPI_at_synapse_FDR$logFC)),] 
PPI_at_synapse_FDR2$order <- 1:nrow(PPI_at_synapse_FDR2)
PPI_at_synapse_FDR2[,c("order",setdiff(names(PPI_at_synapse_FDR2),"order"))]


PPI_at_synapse_FDR2$colour <- ifelse(PPI_at_synapse_FDR2$logFC < 0, "firebrick1","steelblue")
PPI_at_synapse_FDR2$Gene_Expression <- ifelse(PPI_at_synapse_FDR2$logFC < 0, "down-regulated","up-regulated")
PPI_at_synapse_FDR2$hjust <- ifelse(PPI_at_synapse_FDR2$logFC > 0, 1.3, -0.3)

ggplot(PPI_at_synapse_FDR2, aes(order,logFC, label = `Gene symbol`,
                                   hjust = hjust)) + geom_text(aes(y = 0,
                                                                   colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
  aes(show.legend = FALSE)+
  labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
  labs(subtitle="Fold change value for significant genes with FDR < 0.05", 
       title= "protein-protein interaction at synapse associated genes") +
  coord_flip()+
  ylim(-1, 3)+
  theme_classic()

PPI_at_synapse_FDR %>%
  mutate(upregulated = logFC > 0) -> PPI_at_synapse_FDR

ggplot(data = PPI_at_synapse_FDR,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Genes associated to Synapse",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('navyblue'))+
  #theme_minimal()+
  guides(fill = FALSE)



## genes asociated to synapse (combining PPI_at_synapse_FDR and pathcards synaptic vesicle)

synapseTotal <- rbind(PPI_at_synapse_FDR,Synaptic_vesicle_cycle_FDR)
synapseTotal <- unique(synapseTotal)


synapseTotal %>%
  mutate(upregulated = logFC > 0) -> synapseTotal

ggplot(data = synapseTotal,
       aes(x = reorder(`Gene symbol`, logFC), y = logFC,
           fill = upregulated))+
  geom_bar(stat = "identity")+
  #coord_flip()+
  labs(x = "Genes FDR < 0.05", y = "LogFC",
       title = "Genes associated to Synapse and synaptic vesicle",
       subtitles = "Fold change value for significant genes with FDR < 0.05")+
  scale_fill_manual(values=c('red4','navyblue'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #theme_minimal()+
  guides(fill = FALSE)

  
  synapseTotal2 <-synapseTotal[with(synapseTotal, order(synapseTotal$logFC)),] 
  synapseTotal2$order <- 1:nrow(synapseTotal2)
  synapseTotal2[,c("order",setdiff(names(synapseTotal2),"order"))]
  synapseTotal2$Gene_Expression <- ifelse(synapseTotal2$logFC < 0, "down-regulated","up-regulated")
 
  
  synapseTotal2
  ggplot(data = synapseTotal2,
         aes(x = reorder(`Gene symbol`, logFC), y = logFC,
             fill = Gene_Expression))+
    geom_bar(stat = "identity")+
    #coord_flip()+
    labs(x = "Genes FDR < 0.05", y = "LogFC",
         title = "Genes associated to Synapse and synaptic vesicle",
         subtitles = "Fold change value for significant genes with FDR < 0.05")+
    scale_fill_manual(values=c('gray57','goldenrod2'))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
 # guides(fill = FALSE)
  
  
  
  synapseTotal2
  ggplot(data = synapseTotal2,
         aes(x = reorder(`Gene symbol`, logFC), y = logFC,
             fill = Gene_Expression))+
    geom_bar(stat = "identity")+
    #coord_flip()+
    labs(x = "Genes FDR < 0.05", y = "LogFC",
         title = "Genes associated to Synapse and synaptic vesicle",
         subtitles = "Fold change value for significant genes with FDR < 0.05")+
    scale_fill_manual(values=c('gray57','goldenrod2'))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  
  theme_classic() +
    theme(
      axis.line = element_line(colour = 'grey', size = 0.1) ,
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=9, angle = 90),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(size=9),  #summaryPlot$colourPath
      axis.ticks.y = element_line(),
      legend.position="right",
      legend.text = element_text(size = 11, face = "bold") ,
      legend.title = element_blank(),
      legend.key.size = unit(0.7, "cm"),
      legend.key.width = unit(0.6,"cm"),
      plot.title = element_text(size = 15, hjust=0.5, vjust= 1, face = "bold"),
      plot.subtitle = element_blank(),#element_text(size = 2, hjust=0.5)
      strip.text = element_text(size=12, vjust=0.5),
      strip.background = element_rect(fill="lightgray"),
      # panel.border = element_rect(fill = NA, color = "black"),
      panel.spacing.y = unit(0.8, "lines"),
      strip.switch.pad.wrap=unit(20, "lines")
    ) 
  
  
  
################################################################################################
################################################################################################
                                # PATHWAY ANALYSIS (networks)
################################################################################################
################################################################################################
################################################################################################



################################################################################################
                            #    genego_processes_network  #
################################################################################################


library(tidyverse)

#library(readxl)
#snca_mutant_vs_wildtype_rsubread_edger_genego_process_networks_2021 <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/process_networks/snca_mutant_vs_wildtype_rsubread_edger_genego_process_networks_2021.xls")
#View(snca_mutant_vs_wildtype_rsubread_edger_genego_process_networks_2021)
# manually add in order to skipt first two rows 

library(readxl)
snca_mutant_vs_wildtype_rsubread_edger_genego_process_networks_2021 <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/process_networks/snca_mutant_vs_wildtype_rsubread_edger_genego_process_networks_2021.xls", 
                                                                                  skip = 2)
View(snca_mutant_vs_wildtype_rsubread_edger_genego_process_networks_2021)


networks_genego <- snca_mutant_vs_wildtype_rsubread_edger_genego_process_networks_2021
head(networks_genego)

# create input files 

length(networks_genego$Networks)
genebyNetwork <- data.frame(matrix(ncol = 160, nrow = 0))
genebyNetwork <- data.frame(matrix(ncol = 160, nrow = 2000))
genebyNetwork[] = NA
colnames(genebyNetwork) <- networks_genego$Networks

genebyNetwork_intersectData <- genebyNetwork
genebyNetwork_intersectFDRsignificant <- genebyNetwork
genebyNetwork_intersectPvaluesignificant <- genebyNetwork

# obtain all information with filtering 

for (i in 1:length(genebyNetwork)){
  data_a = networks_genego[i,9]
  data_a_a <- separate_rows(data_a,`Network Objects from Active Data`, sep = ", ")
  data_a_a[(length(data_a_a$`Network Objects from Active Data`)+1):2000,]= "NA"   #+1:2000,
  genebyNetwork[,i] <- data_a_a # create general file genes by pathway 
  
  data_a_a_intersectTotal <- data.frame(intersect(data_a_a$`Network Objects from Active Data`, data$`Gene symbol`))
  colnames(data_a_a_intersectTotal)[1] <- "Gene_symbol"
  data_a_a_intersectTotal[(length(data_a_a_intersectTotal$`Gene_symbol`)+1):2000,] = "NA"    #+1:2000,
  genebyNetwork_intersectData[,i] <- data_a_a_intersectTotal
  
  data_a_a_intersectFDRsig <- data.frame(intersect(data_a_a$`Network Objects from Active Data`, data_FRDFrame$`Gene symbol`))
  colnames(data_a_a_intersectFDRsig)[1] <- "Gene_symbol"
  data_a_a_intersectFDRsig[(length(data_a_a_intersectFDRsig$`Gene_symbol`)+1):2000,] = "NA"   #+1:2000,
  genebyNetwork_intersectFDRsignificant[,i] <- data_a_a_intersectFDRsig
  
  data_a_a_intersectPvaluesig <- data.frame(intersect(data_a_a$`Network Objects from Active Data`, data_PvalueFrame$`Gene symbol`))
  colnames(data_a_a_intersectPvaluesig)[1] <- "Gene_symbol"  
  data_a_a_intersectPvaluesig[(length(data_a_a_intersectPvaluesig$`Gene_symbol`)+1):2000,] = "NA"
  genebyNetwork_intersectPvaluesignificant[,i] <- data_a_a_intersectPvaluesig
  
  rm(data_a_a_intersectTotal, data_a_a_intersectFDRsig , data_a_a_intersectPvaluesig )
  rm(data_a, data_a_a )
}


setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/process_networks/")
write.csv(genebyNetwork, file = 'genebyNetwork.csv')
write.csv(genebyNetwork_intersectData, file = 'genebyNetwork_intersectData.csv')
write.csv(genebyNetwork_intersectFDRsignificant, file = 'genebyNetwork_intersectFDRsignificant.csv')
write.csv(genebyNetwork_intersectPvaluesignificant, file = 'genebyNetwork_intersectPvaluesignificant.csv')

######################### plot individual manually selected pathway 
# only one network present a FDR significant 


head(genebyNetwork_intersectFDRsignificant[1,])
#Neurophysiological process_Transmission of nerve impulse
  
  ##########################Neurophysiological process_Transmission of nerve impulse by FFR #########################
  head(genebyNetwork_intersectFDRsignificant)
  nerve_impulse <- data.frame(genebyNetwork_intersectFDRsignificant[,1])
  nerve_impulse <- data.frame(intersect(nerve_impulse$genebyNetwork_intersectFDRsignificant...1., data$`Gene symbol`))
  
  nerve_impulse_sum <- data.frame(matrix(ncol = 6, nrow = 2))
  colnames(nerve_impulse_sum) <- colnames(data)
  
  for (i in 1:length(nerve_impulse$intersect.nerve_impulse.genebyNetwork_intersectFDRsignificant...1...)){
    nerve_impulse_sum[i,] = filter(data,`Gene symbol` == nerve_impulse$intersect.nerve_impulse.genebyNetwork_intersectFDRsignificant...1...[i])}
  
  #colnames(nerve_impulse_sum) <- colnames(data)
  head(nerve_impulse_sum)
  
  # order logFC
  nerve_impulse_sum2 <-nerve_impulse_sum[with(nerve_impulse_sum, order(nerve_impulse_sum$logFC)),] 
  nerve_impulse_sum2$order <- 1:nrow(nerve_impulse_sum2)
  nerve_impulse_sum2[,c("order",setdiff(names(nerve_impulse_sum),"order"))]
  
  # plotting 
  library(ggplot2) 
  
  # add colour & position labels & levels
  nerve_impulse_sum2$colour <- ifelse(nerve_impulse_sum2$logFC < 0, "firebrick1","steelblue")
  nerve_impulse_sum2$Gene_Expression <- ifelse(nerve_impulse_sum2$logFC < 0, "down-regulated","up-regulated")
  nerve_impulse_sum2$hjust <- ifelse(nerve_impulse_sum2$logFC > 0, 1.3, -0.3)
  
  ggplot(nerve_impulse_sum2, aes(order,logFC, label = `Gene symbol`,
                                            hjust = hjust)) + geom_text(aes(y = 0,
                                                                            colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
    aes(show.legend = FALSE)+
    labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
    labs(subtitle="Fold change value for significant genes with FDR  < 0.05", 
         title= "Neurophysiological process_Transmission of nerve impulse") +
    ylim(-2, 5)+
    coord_flip() 
  
  rm(nerve_impulse,nerve_impulse_sum )
  
  
  ##########################Neurophysiological process_Transmission of nerve impulse by Pvalue #########################
  
  nerve_impulse <- data.frame(genebyNetwork_intersectPvaluesignificant[,1])
  nerve_impulse <- data.frame(intersect(nerve_impulse$genebyNetwork_intersectPvaluesignificant...1., data$`Gene symbol`))
  
  nerve_impulse_sum <- data.frame(matrix(ncol = 6, nrow = 11))
  colnames(nerve_impulse_sum) <- colnames(data)
  
  for (i in 1:length(nerve_impulse$intersect.nerve_impulse.genebyNetwork_intersectPvaluesignificant...1...)){
    nerve_impulse_sum[i,] = filter(data,`Gene symbol`  == nerve_impulse$intersect.nerve_impulse.genebyNetwork_intersectPvaluesignificant...1...[i])}
  
  #colnames(nerve_impulse_sum) <- colnames(data)
  head(nerve_impulse_sum)
  
  # order logFC
  nerve_impulse_sum2 <-nerve_impulse_sum[with(nerve_impulse_sum, order(nerve_impulse_sum$logFC)),] 
  nerve_impulse_sum2$order <- 1:nrow(nerve_impulse_sum2)
  nerve_impulse_sum2[,c("order",setdiff(names(nerve_impulse_sum),"order"))]
  
  # plotting 
  library(ggplot2) 
  
  # add colour & position labels & levels
  nerve_impulse_sum2$colour <- ifelse(nerve_impulse_sum2$logFC < 0, "firebrick1","steelblue")
  nerve_impulse_sum2$Gene_Expression <- ifelse(nerve_impulse_sum2$logFC < 0, "down-regulated","up-regulated")
  nerve_impulse_sum2$hjust <- ifelse(nerve_impulse_sum2$logFC > 0, 1.3, -0.3)
  
  ggplot(nerve_impulse_sum2, aes(order,logFC, label = `Gene symbol`,
                                 hjust = hjust)) + geom_text(aes(y = 0,
                                                                 colour = Gene_Expression)) +geom_bar(stat = "identity", aes(fill = Gene_Expression))+
    aes(show.legend = FALSE)+
    labs( x = "Ranked genes based on Fold Change", y="log Fold Change") +
    labs(subtitle="Fold change value for significant genes with Pvalue < 0.05", 
         title= "Neurophysiological process_Transmission of nerve impulse") +
    ylim(-2, 5)+
    coord_flip() 
  
  nerve_impulse_sum2 <- nerve_impulse_sum2[1:10,]
 #snca_mutant_vs_wildtype_rsubread_edger_genego_process_networks_2021  -> transmission of neerv impulse
  ggplot(data = nerve_impulse_sum2,
         aes(x = reorder(`Gene symbol`, logFC), y = logFC,
             fill = Gene_Expression))+
    geom_bar(stat = "identity")+
    #coord_flip()+
    labs(x = "Genes FDR < 0.05", y = "LogFC",
         title = "Genes associated to transmission of nerve impulse",
         subtitles = "Fold change value for significant genes with FDR < 0.05")+
    scale_fill_manual(values=c('navyblue')) +
    #theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ############################################################################################################
  ############################################################################################################
  
  ################################ PATHWAY ANALYSIS (networks) summay per pathway ############################
  
  ############################################################################################################
  ############################################################################################################
  # 
  
  
  # pathways not sinificant ( information fom summary not needed)
  summary_network <- data.frame(matrix(ncol = 10, nrow = 160)) #length networks 
  colnames(summary_network) <-c("Networks" , "total_genes", "total_genes_mean" ,"total_genes_media", "Pvalue_genes", "Pvalue_genes_mean" ,"Pvalue_genes_media", "FDR_genes", "FDR_genes_mean" ,"FDR_genes_media")
  summary_network[,1] <- networks_genego$Networks
  
  
  for (i in 1:length(summary_network$Networks)){
    # for all data 
    target <- c(intersect(genebyNetwork_intersectData[,i], data$`Gene symbol`))
    temp <- filter(data, `Gene symbol` %in% target)
    summary_network[i,2] = length(temp$`Gene symbol`)
    summary_network[i,3] = mean(temp$logFC)
    summary_network[i,4] = median(temp$logFC)
    rm(target,temp)
    # for pvalue
    target <- c(intersect(genebyNetwork_intersectPvaluesignificant[,i], data$`Gene symbol`))
    temp <- filter(data, `Gene symbol` %in% target)
    summary_network[i,5] = length(temp$`Gene symbol`)
    summary_network[i,6] = mean(temp$logFC)
    summary_network[i,7] = median(temp$logFC)
    # for FDR
    target <- c(intersect(genebyNetwork_intersectFDRsignificant[,i], data$`Gene symbol`))
    temp <- filter(data, `Gene symbol` %in% target)
    summary_network[i,8] = length(temp$`Gene symbol`)
    summary_network[i,9] = mean(temp$logFC)
    summary_network[i,10] = median(temp$logFC)
  } 
  
  setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/process_networks/")
  write.csv(summary_network, file = 'summary_network.csv')
  
  
################################################################################################
################################################################################################
                                     # GO maps

library(readxl)
snca_mutant_vs_wildtype_rsubread_edger_genego_pathway_maps_2021 <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/pathway_maps/snca_mutant_vs_wildtype_rsubread_edger_genego_pathway_maps_2021.xls", 
                                                                                skip = 2)
genego_pathway_maps <-snca_mutant_vs_wildtype_rsubread_edger_genego_pathway_maps_2021


# extract patwhays with FDR <0.05 
# only first 16 have significant FDR 
#length(which(genego_pathway_maps$FDR < 0.05)) # number of gene with FDR lower than 0.05 
#index = rownames(genego_pathway_maps)[which(genego_pathway_maps$FDR < 0.05)]# list pathways 
genego_pathway_maps_FDR  <- genego_pathway_maps %>% filter(genego_pathway_maps$FDR < 0.05)
setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/pathway_maps/")
write.csv(data_PvalueFrame, file = 'genego_pathway_maps_FDR.csv')

######################
# obtain summary genes sigf FDR<0.05 - media and average 

length(genego_pathway_maps_FDR$Maps)
GO_genego_pathway_maps_FDR <- data.frame(matrix(ncol = 16, nrow = 0)) # 601 to extract only  patwhays with FDR <0.05 (significance)
GO_genego_pathway_maps_FDR <- data.frame(matrix(ncol = 16, nrow = 20))
GO_genego_pathway_maps_FDR[] = NA
colnames(GO_genego_pathway_maps_FDR) <- genego_pathway_maps_FDR$Maps[1:16]

GO_genego_pathway_maps_intersectFDRsignificant <- GO_genego_pathway_maps_FDR
GO_genego_pathway_maps__intersectPvaluesignificant <- GO_genego_pathway_maps_FDR



######################
# obtain all information with filtering 

for (i in 1:length(GO_genego_pathway_maps_FDR)){
  data_a = genego_pathway_maps_FDR[i,9]
  data_a_a <- separate_rows(data_a,`Network Objects from Active Data`, sep = ", ")
  data_a_a[(length(data_a_a$`Network Objects from Active Data`)+1):20,] = "empty"
  GO_genego_pathway_maps_FDR[,i] <- data_a_a # create general file genes by pathway 
  
  data_a_a_intersectFDRsig <- data.frame(intersect(data_a_a$`Network Objects from Active Data`, data_FRDFrame$`Gene symbol`))
  colnames(data_a_a_intersectFDRsig)[1] <- "Gene_symbol"
  data_a_a_intersectFDRsig[(length(data_a_a_intersectFDRsig$Gene_symbol)+1):20,] = "empty"
  GO_genego_pathway_maps_intersectFDRsignificant[,i] <- data_a_a_intersectFDRsig
  
  data_a_a_intersectPvaluesig <- data.frame(intersect(data_a_a$`Network Objects from Active Data`, data_PvalueFrame$`Gene symbol`))
  colnames(data_a_a_intersectPvaluesig)[1] <- "Gene_symbol"
  data_a_a_intersectPvaluesig[(length(data_a_a_intersectPvaluesig$Gene_symbol)+1):20,] = "empty"
  GO_genego_pathway_maps__intersectPvaluesignificant[,i] <- data_a_a_intersectPvaluesig
  
  rm(data_a_a_intersectTotal, data_a_a_intersectFDRsig , data_a_a_intersectPvaluesig )
  rm(data_a, data_a_a )
}

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/pathway_maps/")
write.csv(data_PvalueFrame, file = 'GO_genego_pathway_maps_FDR.csv')

######################
# obtain mean and media for pathways deregulated FDR <0.05 using genes with FDR also 0>05 

pathway_maps_FDR <- data.frame(matrix(ncol = 10, nrow = 16)) #length networks 
colnames(pathway_maps_FDR) <-c("Networks" , "total_genes", "total_genes_mean" ,"total_genes_media", "Pvalue_genes", "Pvalue_genes_mean" ,"Pvalue_genes_media", "FDR_genes", "FDR_genes_mean" ,"FDR_genes_media")
pathway_maps_FDR[,1] <- genego_pathway_maps_FDR$Maps[1:16] # obtaining onlt those with FDR <0.05

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/pathway_maps/")
pathway_maps_FDR <- read.csv("genego_pathway_maps_FDR.csv") 

for (i in 1:length(pathway_maps_FDR$Networks)){
  # for all data 
  target <- c(intersect(GO_genego_pathway_maps_FDR[,i], data$`Gene symbol`))
  temp <- filter(data, `Gene symbol` %in% target)
  pathway_maps_FDR[i,2] = length(temp$`Gene symbol`) # calculate number genes 
  pathway_maps_FDR[i,3] = mean(temp$logFC) # calculat mean logFC 
  pathway_maps_FDR[i,4] = median(temp$logFC) # calculate media logFC 
  rm(target,temp)
  # for pvalue
  target <- c(intersect(GO_genego_pathway_maps__intersectPvaluesignificant[,i], data$`Gene symbol`))
  temp <- filter(data, `Gene symbol` %in% target)
  pathway_maps_FDR[i,5] = length(temp$`Gene symbol`)
  pathway_maps_FDR[i,6] = mean(temp$logFC)
  pathway_maps_FDR[i,7] = median(temp$logFC)
  # for FDR
  target <- c(intersect(GO_genego_pathway_maps_intersectFDRsignificant[,i], data$`Gene symbol`))
  temp <- filter(data, `Gene symbol` %in% target)
  pathway_maps_FDR[i,8] = length(temp$`Gene symbol`)
  pathway_maps_FDR[i,9] = mean(temp$logFC)
  pathway_maps_FDR[i,10] = median(temp$logFC)
} 

## plotting all patwhays 


pathway_maps_FDR2<-pathway_maps_FDR[-c(10),] %>% # row number 10 is empty because it contains no sing genes 
  arrange(desc(pathway_maps_FDR[-c(10),])) 


# order logFC
pathway_maps_FDR3 <-pathway_maps_FDR2[with(pathway_maps_FDR2, order(pathway_maps_FDR2$FDR_genes_media)),] 
pathway_maps_FDR3$order <- 1:nrow(pathway_maps_FDR3)
pathway_maps_FDR3[,c("order",setdiff(names(pathway_maps_FDR3),"order"))]

# plotting 
library(ggplot2) 

# add colour & position labels & levels
pathway_maps_FDR3$colour <- ifelse(pathway_maps_FDR3$FDR_genes_media < 0, "firebrick1","steelblue")
pathway_maps_FDR3$Gene_Expression <- ifelse(pathway_maps_FDR3$FDR_genes_media < 0, "down-regulated","up-regulated")
pathway_maps_FDR3$hjust <- ifelse(pathway_maps_FDR3$FDR_genes_media > 0, 1.3, -0.3)

ggplot(pathway_maps_FDR3, aes(order,FDR_genes_media, label = Networks,
                               hjust = hjust)) + 
  geom_text(aes(y = -1,colour = FDR_genes_media)) +
  geom_bar(stat = "identity", aes(fill = FDR_genes_media))+
  geom_text(aes(label=FDR_genes_media))+
  aes(show.legend = FALSE)+
  labs( x = "genego pathway maps FDR<0.05", y="log Fold Change") +
  labs(subtitle="pathways with FDR<0.05 & media included genes with FDR < 0.05", 
       title= "Significantly differentiated genego pathways maps") +
  ylim(-2, 5)+
  coord_flip() 

################################################################################################
################################################################################################
                                    # geneonthology preocesses 


library(readxl)
snca_mutant_vs_wildtype_rsubread_edger_genego_geneontology_processes_2021 <- read_excel("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/geneontology_processes/snca_mutant_vs_wildtype_rsubread_edger_genego_geneontology_processes_2021.xls", 
                                                                                        skip = 2)
View(snca_mutant_vs_wildtype_rsubread_edger_genego_geneontology_processes_2021)

GO <- snca_mutant_vs_wildtype_rsubread_edger_genego_geneontology_processes_2021
# manually inserted -> skip two first 
# only first 601 have significant FDR 

# GO structure, GO pathways as column name, and genes as row 
library(tidyverse)

# create input files 

length(GO$Processes)
GO_genebyProcess <- data.frame(matrix(ncol = 601, nrow = 0)) # 601 to extract only  patwhays with FDR <0.05 (significance)
GO_genebyProcess <- data.frame(matrix(ncol = 601, nrow = 700))
GO_genebyProcess[] = NA
colnames(GO_genebyProcess) <- GO$Processes[1:601]

GO_genebyProcess_intersectData <- GO_genebyProcess
GO_genebyProcess_intersectFDRsignificant <- GO_genebyProcess
GO_genebyProcess_intersectPvaluesignificant <- GO_genebyProcess

# obtain all information with filtering 

for (i in 1:length(GO_genebyProcess)){
  data_a = GO[i,9]
  data_a_a <- separate_rows(data_a,`Network Objects from Active Data`, sep = ", ")
  data_a_a[(length(data_a_a$`Network Objects from Active Data`)+1):700,] = "empty"
  GO_genebyProcess[,i] <- data_a_a # create general file genes by pathway 
  
  #data_a_a_intersectTotal<- data.frame(matrix(ncol = 1, nrow = 2200))
  data_a_a_intersectTotal <- data.frame(intersect(data_a_a$`Network Objects from Active Data`, data$`Gene symbol`))
  colnames(data_a_a_intersectTotal)[1] <- "Gene_symbol"
  data_a_a_intersectTotal[(length(data_a_a_intersectTotal$Gene_symbol)+1):700,] = "empty"
  GO_genebyProcess_intersectData[,i] <- data_a_a_intersectTotal
  
  data_a_a_intersectFDRsig <- data.frame(intersect(data_a_a$`Network Objects from Active Data`, data_FRDFrame$`Gene symbol`))
  colnames(data_a_a_intersectFDRsig)[1] <- "Gene_symbol"
  data_a_a_intersectFDRsig[(length(data_a_a_intersectFDRsig$Gene_symbol)+1):700,] = "empty"
  GO_genebyProcess_intersectFDRsignificant[,i] <- data_a_a_intersectFDRsig
  
  data_a_a_intersectPvaluesig <- data.frame(intersect(data_a_a$`Network Objects from Active Data`, data_PvalueFrame$`Gene symbol`))
  colnames(data_a_a_intersectPvaluesig)[1] <- "Gene_symbol"
  data_a_a_intersectPvaluesig[(length(data_a_a_intersectPvaluesig$Gene_symbol)+1):700,] = "empty"
  GO_genebyProcess_intersectPvaluesignificant[,i] <- data_a_a_intersectPvaluesig
  
  rm(data_a_a_intersectTotal, data_a_a_intersectFDRsig , data_a_a_intersectPvaluesig )
  rm(data_a, data_a_a )
}


setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/geneontology_processes/")
write.csv(GO_genebyProcess, file = 'GO_genebyProcess.csv')
write.csv(GO_genebyProcess_intersectData, file = 'GO_genebyProcess_intersectData.csv')
write.csv(GO_genebyProcess_intersectFDRsignificant, file = 'GO_genebyProcess_intersectFDRsignificant.csv')
write.csv(GO_genebyProcess_intersectPvaluesignificant, file = 'GO_genebyProcess_intersectPvaluesignificant.csv')


################################################################################################
################################ PATHWAY ANALYSIS (networks) summay per pathway
# obtain all total genes crossed with logFC file (obtain total number and mean and media of logFC for those 
# that are significantly exprrssed with FDR <0.05)

summary_GO_processes <- data.frame(matrix(ncol = 10, nrow = 601)) #length networks 
colnames(summary_GO_processes) <-c("Networks" , "total_genes", "total_genes_mean" ,"total_genes_media", "Pvalue_genes", "Pvalue_genes_mean" ,"Pvalue_genes_media", "FDR_genes", "FDR_genes_mean" ,"FDR_genes_media")
summary_GO_processes[,1] <- GO$Processes[1:601] # obtaining onlt those with FDR <0.05


for (i in 1:length(summary_GO_processes$Networks)){
  # for all data 
  target <- c(intersect(GO_genebyProcess_intersectData[,i], data$`Gene symbol`))
  temp <- filter(data, `Gene symbol` %in% target)
  summary_GO_processes[i,2] = length(temp$`Gene symbol`) # calculate number genes 
  summary_GO_processes[i,3] = mean(temp$logFC) # calculat mean logFC 
  summary_GO_processes[i,4] = median(temp$logFC) # calculate media logFC 
  rm(target,temp)
  # for pvalue
  target <- c(intersect(GO_genebyProcess_intersectPvaluesignificant[,i], data$`Gene symbol`))
  temp <- filter(data, `Gene symbol` %in% target)
  summary_GO_processes[i,5] = length(temp$`Gene symbol`)
  summary_GO_processes[i,6] = mean(temp$logFC)
  summary_GO_processes[i,7] = median(temp$logFC)
  # for FDR
  target <- c(intersect(GO_genebyProcess_intersectFDRsignificant[,i], data$`Gene symbol`))
  temp <- filter(data, `Gene symbol` %in% target)
  summary_GO_processes[i,8] = length(temp$`Gene symbol`)
  summary_GO_processes[i,9] = mean(temp$logFC)
  summary_GO_processes[i,10] = median(temp$logFC)
} 

setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/geneontology_processes/")
write.csv(summary_GO_processes, file = 'summary_GO_processes.csv')

################################################################################################
################################ PATHWAY ANALYSIS (networks) summay per pathway - 
# obtain number of genes per pathway and if they are positive or negatively regulated 

summary_GO <- data.frame(matrix(ncol = 10, nrow = 601)) #length networks 
colnames(summary_GO) <-c("Networks" , "total_genes", "logFC_total_pos" ,"logFC_total_neg", "Pvalue_genes", "logFC_pvalue_pos" ,"logFC_pvalue_neg", "FDR_genes", "logFC_FDR_pos" ,"logFC_FDR_neg")
summary_GO[,1] <- GO$Processes[1:601] # obtaining onlt those with FDR <0.05


for (i in 1:length(summary_GO$Networks)){
  # for all data 
  target <- c(intersect(GO_genebyProcess_intersectData[,i], data$`Gene symbol`))
  temp <- filter(data, `Gene symbol` %in% target)
  summary_GO[i,2] = length(temp$`Gene symbol`) # calculate number genes 
  pos <- temp %>% filter(temp$logFC > 0)
  neg <- temp %>% filter(temp$logFC < 0)
  summary_GO[i,3] = length(pos$`Gene symbol`) # calculat mean logFC 
  summary_GO[i,4] = length(neg$`Gene symbol`) # calculate media logFC 
  rm(target,temp, pos, neg)
  # for pvalue
  target <- c(intersect(GO_genebyProcess_intersectPvaluesignificant[,i], data$`Gene symbol`))
  temp <- filter(data, `Gene symbol` %in% target)
  summary_GO[i,5] = length(temp$`Gene symbol`) # calculate number genes 
  pos <- temp %>% filter(temp$logFC > 0)
  neg <- temp %>% filter(temp$logFC < 0)
  summary_GO[i,6] = length(pos$`Gene symbol`) # calculat mean logFC 
  summary_GO[i,7] = length(neg$`Gene symbol`) # calculate media logFC 
  rm(target,temp, pos, neg)
  # for FDR
  target <- c(intersect(GO_genebyProcess_intersectFDRsignificant[,i], data$`Gene symbol`))
  temp <- filter(data, `Gene symbol` %in% target)
  summary_GO[i,8] = length(temp$`Gene symbol`) # calculate number genes 
  pos <- temp %>% filter(temp$logFC > 0)
  neg <- temp %>% filter(temp$logFC < 0)
  summary_GO[i,9] = length(pos$`Gene symbol`) # calculat mean logFC 
  summary_GO[i,10] = length(neg$`Gene symbol`) # calculate media logFC 
  rm(target,temp, pos, neg)
} 


## manually selecting pathways of interest 

# from already selected file (summary_GO_processes_manuallyAdapted)

target <- c(intersect(summary_GO_processes_manuallyAdapted$Networks, summary_GO$Networks))
summary_GO_selected <- filter(summary_GO, Networks %in% target)
# save and manually add the clasification 
setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/geneontology_processes/")
write.csv(summary_GO_selected, file = 'summary_GO_selected.csv')
#load modify version 
setwd("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/geneontology_processes/")
summary_GO_selected <- read.csv("summary_GO_selected_manuallyAdapted.csv") 


# prepare frame for plotting 
positiveFDR <- data.frame(summary_GO_selected[,2],summary_GO_selected[,9:10],summary_GO_selected[,12])
negativeFDR <- data.frame(summary_GO_selected[,2],summary_GO_selected[,9],summary_GO_selected[,11:12])
colnames(positiveFDR) <-c("Networks" ,  "FDR_genes", "FDR_genes_value" ,"category" )
colnames(negativeFDR) <-c("Networks" ,  "FDR_genes", "FDR_genes_value" ,"category")
negativeFDR$FDR_genes_value = negativeFDR$FDR_genes_value*(-1)
summaryPlot <- rbind(positiveFDR,negativeFDR)


#summaryPlot2<- arrange(desc(summaryPlot$FDR_genes_value)) 
summaryPlot$colour <- ifelse(summaryPlot$FDR_genes_value < 0, 'downregulated','upregulated')
summaryPlot$colourPath <- summaryPlot$category
#a <- unique(summaryPlot$colourPath)
summaryPlot$colourPath <- str_replace_all(summaryPlot$colourPath, 'neuronal differentiation', 'darkcyan')
summaryPlot$colourPath <- str_replace_all(summaryPlot$colourPath, 'neuronal projection', 'deeppink3')
summaryPlot$colourPath <- str_replace_all(summaryPlot$colourPath, 'synapse regulation', 'darkmagenta')


# 
# library(ggplot2)
# ggplot(data = summaryPlot) + 
#   geom_bar(aes(x=reorder(Networks, desc(FDR_genes_value)),y=FDR_genes_value, fill=colourPath),stat="identity",position="identity") +
#   geom_text(aes(x=Networks,y=FDR_genes_value,label=abs(FDR_genes_value)))+
#   #geom_text(aes(label=abs(FDR_genes_value)), vjust=1.5, colour="white", size=3.5) + #,vjust = ifelse(summaryPlot2$FDR_genes_value >= 0, 0, 1)) +
#   #scale_y_continuous(labels=abs)+
#   coord_flip() +
#   theme_classic() 
# 
# 
summaryPlot2 <-summaryPlot[with(summaryPlot, order(summaryPlot$FDR_genes)),] 
summaryPlot2$order <- 1:nrow(summaryPlot2)
summaryPlot2[,c("order",setdiff(names(summaryPlot2),"order"))]

co <- data.frame(summaryPlot2$Networks,summaryPlot2$colourPath)
co2 <-unique(co)
color1 <- as.vector(co2[,2]) 
colors2 <- rev(color1)

# # octain vector colors 
# summary_GO_selected2 <-summary_GO_selected[with(summary_GO_selected, order(desc(summary_GO_selected$FDR_genes))),] 
# summary_GO_selected2$order <- 1:nrow(summary_GO_selected2)
# summary_GO_selected2[,c("order",setdiff(names(summary_GO_selected2),"order"))]
# summary_GO_selected2$colourPath <- summary_GO_selected2$category
# summary_GO_selected2$colourPath <- str_replace_all(summary_GO_selected2$colourPath, 'neuronal differentiation', 'darkcyan')
# summary_GO_selected2$colourPath <- str_replace_all(summary_GO_selected2$colourPath, 'neuronal projection', 'deeppink3')
# summary_GO_selected2$colourPath <- str_replace_all(summary_GO_selected2$colourPath, 'synapse regulation', 'darkmagenta')
# 
# color <- c(summary_GO_selected2$colourPath)



ggplot(data = summaryPlot2) + 
  geom_bar(aes(x=reorder(Networks, desc(FDR_genes)),y=FDR_genes_value, fill=colour),stat="identity",position="identity") +
  scale_fill_manual(values=c('gray57','goldenrod2')) +
  geom_text(aes(x=Networks,y=FDR_genes_value,label=abs(FDR_genes_value))) +
  #geom_text(aes(label=abs(FDR_genes_value)), vjust=1.5, colour="white", size=3.5) + #,vjust = ifelse(summaryPlot2$FDR_genes_value >= 0, 0, 1)) +
  #scale_y_continuous(labels=abs)+
  coord_flip() +
   labs(x   ="significantly deregualted pathways",
     y     = "number of deregulated genes (>0< logFC)",
     fill  = "category",
     title = "Total deregulated genes per pathway" ) +
  theme_classic() +
  theme(
    axis.line = element_line(colour = 'grey', size = 0.1) ,
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=9),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size=9, color=colors2),  #summaryPlot$colourPath
    axis.ticks.y = element_line(),
    legend.position="right",
    legend.text = element_text(size = 11, face = "bold") ,
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.key.width = unit(0.6,"cm"),
    plot.title = element_text(size = 15, hjust=0.5, vjust= 1, face = "bold"),
    plot.subtitle = element_blank(),#element_text(size = 2, hjust=0.5)
    strip.text = element_text(size=12, vjust=0.5),
    strip.background = element_rect(fill="lightgray"),
   # panel.border = element_rect(fill = NA, color = "black"),
    panel.spacing.y = unit(0.8, "lines"),
    strip.switch.pad.wrap=unit(20, "lines")
  )  -> p




################################ PATHWAY ANALYSIS (networks) summay per pathway
## plotting of all pathways by sub-clasification 


summary_GO_processes_manuallyAdapted <- read.csv("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis genego/geneontology_processes/summary_GO_processes_manuallyAdapted.csv", row.names=1)
View(summary_GO_processes_manuallyAdapted)


summary_GO_processes_manuallyAdapted2<-summary_GO_processes_manuallyAdapted %>%
  arrange(desc(FDR_genes_media))

ggplot(data=summary_GO_processes_manuallyAdapted2, aes(x=Networks, y=FDR_genes_media, fill = pathway)) +
  geom_bar(stat="identity")+
  #geom_text(aes(label=Gene_symbol))+
  scale_colour_gradient2()+
  coord_flip()+
  #ylim(0, 7)+
  scale_x_discrete(limits = summary_GO_processes_manuallyAdapted2$Networks)+
  labs( x = "Neuronal associated and deregulated GO processes ", y="media logFC genes within pathway (FDR>0.05)") +
  labs(subtitle="deregulated GO processes (FDR<0.05)", 
       title= "Neuronal relevant pathwyas") +
  scale_fill_manual(values=c('aquamarine4','darkgoldenrod4','darkslategray4','deeppink','darkslategray3','aquamarine3','darkgoldenrod3','navyblue','darkgoldenrod3','deeppink4'))+
  theme_classic()

###################################################################################################
# load previous data 
load("//atlas.uni.lux/users/jennifer.modamio/DVB_group/2021/transcriptomics/post-analysis logFoldChange/R scrip gene level analysis 05252021 cleaned INPUT DATA.RData")

###################################################################################################
# re-plotting  
# Synaptic_vesicle_cycle -> Synaptic_vesicle_cycle_FDR
# pathcards_PPI_at_synapse_genes ->  PPI_at_synapse_FDR


PPI_at_synapse_FDR%>%
  arrange(desc(logFC)) ->PPI_at_synapse_FDR2


ggplot(data=PPI_at_synapse_FDR2, aes(x=`Gene symbol`, y=logFC, fill = `Gene symbol` )) +
  geom_bar(stat="identity")+
  #geom_text(aes(label=Gene_symbol))+
  scale_colour_gradient2()+
  coord_flip()+
  #ylim(0, 7)+
  scale_x_discrete(limits = PPI_at_synapse_FDR2$`Gene symbol`)+
  labs( x = "genes LogFC (FDR<0.05)", y="Synapse related GO processes") +
  labs(subtitle="logFC for genes with  FDR < 0.05", 
       title= "Synapse") +
  scale_fill_manual(values=c('red4','navyblue'))+
  theme_classic()




###  synapse pathway
synapse2<-PPI_at_synapse_FDR %>%
  synapse2<-arrange(desc(synapse2))

ggplot(data=synapse2, aes(x=`Gene symbol`, y=logFC, fill = `Gene symbol` )) +
  geom_bar(stat="identity")+
  #geom_text(aes(label=Gene_symbol))+
  scale_colour_gradient2()+
  coord_flip()+
  #ylim(0, 7)+
  scale_x_discrete(limits = synapse2$`Gene symbol`)+
  labs( x = "Synapse associated and deregulated GO processes ", y="media logFC genes within pathway (FDR>0.05)") +
  labs(subtitle="deregulated GO processes (FDR<0.05)", 
       title= "Synapse relevant pathwyas") +
  scale_fill_manual(values=c('navyblue'))+
  theme_classic()