library(tidyverse)
library(DESeq2)
library(ggplot2)
library(gprofiler2)
library(Cairo)
library(rstatix)
library(ggsignif)
library(ggpubr)




## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
ncol(raw)

counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts

## Organize patient IDs
ids <- names(raw)[-1]
ids
ids <- data.frame(sample_id = ids,
                  patient = str_extract(ids, "[0-9]{4}"),
                  visit = str_extract(ids, "_[0-9]_"),
                  type = str_extract(ids, "_[A-Z]{2}_"),
                  stringsAsFactors = F)
ids
ids$visit <- str_extract(ids$visit, "[0-9]")
ids$type <- str_extract(ids$type, "[A-Z]+")
ids$pairing <- paste0(ids$patient, ids$type)
ids
ids1<-ids # for ISS analysis

ids%>%filter(visit >5)
ids%>%filter(patient %in% c("1056", "1201", "1024"))
ids%>%filter(duplicated(patient))

## Filter for samples from newly diagnosed and every new bone marrow biopsy
ids<-ids[ids$visit %in% c("1", "2", "3", "4", "5", "6"), ]
ids
counts <- counts[, colnames(counts) %in% ids$sample_id]
counts

## Transform every sample as factor
ids$visit <- factor(ids$visit, levels = c("1", "2", "3", "4", "5", "6"))
ids$pairing <- factor(ids$pairing)
ids%>%filter(patient == 1201)

## Filter for paired samples
ids <- ids[ids$pairing %in% names(table(ids$pairing)[table(ids$pairing) >1]), ]
ids

ids%>%filter(patient == 2695) # Examples
ids%>%filter(patient == 1024)

counts2 <- counts[, colnames(counts) %in% ids$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(ids)

counts2<- as.matrix(counts2)
ncol(counts2)
nrow(ids)

counts2

## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = ids, design = ~ visit)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC_2vs1 <- lfcShrink(dds, coef = "visit_2_vs_1", type = "apeglm")   # Fold change shrinkage for visualization
resLFC_3vs1 <- lfcShrink(dds, coef = "visit_3_vs_1", type = "apeglm")
resLFC_4vs1 <- lfcShrink(dds, coef = "visit_4_vs_1", type = "apeglm")
resLFC_5vs1 <- lfcShrink(dds, coef = "visit_5_vs_1", type = "apeglm")
resLFC_6vs1 <- lfcShrink(dds, coef = "visit_6_vs_1", type = "apeglm")

#pdf(" MA plot shrink.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot Newly diagnosed vs Relapse")
plotMA(resLFC_2vs1, ylim=c(-2,2), main = "MA plot visit_2_vs_1")
plotMA(resLFC_3vs1, ylim=c(-2,2), main = "MA plot visit_3_vs_1")
plotMA(resLFC_4vs1, ylim=c(-2,2), main = "MA plot visit_4_vs_1")
plotMA(resLFC_5vs1, ylim=c(-2,2), main = "MA plot visit_5_vs_1")
plotMA(resLFC_6vs1, ylim=c(-2,2), main = "MA plot visit_6_vs_1")
#dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("visit"))

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))
temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)

# Select names by extracting number of the sample
x<-str_extract(names(temp), "_[0-9]_")[1:203] 
x<-str_extract(x, "[0-9]")[1:203]
x

temp_bar <- data.frame(y = as.numeric(temp[temp$GENE %in% "PAX5", 1:203]), 
                       x = x,                
                       stringsAsFactors = F)
temp_bar

# Determing significant p-values for showing it up in the plot only the significant ones
stat_pvalue <- temp_bar %>% 
  rstatix::wilcox_test(y ~ x) %>%
  #filter(p < 0.05) %>% 
  #rstatix::add_significance("p") %>% 
  rstatix::add_y_position() %>% 
  mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))



# PLOT variance of PAX5 expression across every sample

temp_bar%>%ggplot(aes(x = x, y = y)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "Patient Status",
                   breaks = 1:7,
                   labels = c("Baseline","2nd BM", "3rd BM", "4th BM", "5th BM", "6th BM", "7th BM")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% c("PAX5")],2),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 3))) +
  #ggpubr::stat_pvalue_manual(stat_pvalue, label = "p.signif") +
  theme_minimal()

#############################

#############################
## Filter for samples from newly diagnosed and Relapse
ids_ND_relapse<-ids[ids$visit %in% c("1", "2", "3", "4", "5", "6"), ]
ids_ND_relapse$visit<-as.numeric(ids_ND_relapse$visit)
ids_ND_relapse$visit<-replace(ids_ND_relapse$visit, ids_ND_relapse$visit>2,2)
nrow(ids_ND_relapse)

counts <- counts[, colnames(counts) %in% ids_ND_relapse$sample_id]
counts

## newly diagnosed and relapse as a factors
ids_ND_relapse$visit <- factor(ids_ND_relapse$visit, levels = c("1", "2"))
ids_ND_relapse$pairing <- factor(ids_ND_relapse$pairing)
ids_ND_relapse%>%filter(patient == 1201)

## Filter for paired samples
ids_ND_relapse<- ids_ND_relapse[ids_ND_relapse$pairing %in% names(table(ids_ND_relapse$pairing)[table(ids_ND_relapse$pairing) ==2]), ]
ids_ND_relapse
nrow(ids_ND_relapse)

ids_ND_relapse%>%filter(patient == 2695) # Examples
ids_ND_relapse%>%filter(patient == 1024)

counts2 <- counts[, colnames(counts) %in% ids_ND_relapse$sample_id]
counts
nrow(counts2)
ncol(counts2)

counts2<- as.matrix(counts2)
ncol(counts2)
nrow(ids)


## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = ids_ND_relapse, design = ~ visit)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC_2vs1 <- lfcShrink(dds, coef = "visit_2_vs_1", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink Newly diagnosed vs Relapse.pdf", width=6.72, height=4.07)
plotMA(res, ylim=c(-2,2), main = "MA plot Newly diagnosed vs Relapse")
plotMA(resLFC_2vs1, ylim=c(-2,2), main = "MA plot Newly diagnosed vs Relapse resLFC")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("visit"))

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization for Newly Diagnosed and relapsed
temp <- as.data.frame(assay(vsd))
temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
temp

temp%>%filter(GENE == "PAX5")
names(temp)
temp[temp$GENE == "PAX5", 1:162]
x<-str_extract(names(temp), "_[0-9]_")[1:162]
x<-str_extract(x, "[0-9]")[1:162]
x

temp_bar <- data.frame(y = as.numeric(temp[temp$GENE %in% c("PAX5"), 1:162]), 
                       x = x,                
                       stringsAsFactors = F)

temp_bar

temp_bar1<-temp_bar
temp_bar1
temp_bar1$x<-replace(temp_bar1$x, temp_bar1$x>2, 2)
temp_bar1

res$pvalue[res$GENE == "PAX5"]

## BOX PLOT FOR NEWLY DIAGNOSED VS RELAPSED
temp_bar1%>%ggplot(aes(x = x, y = y)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "Patient Status",
                   breaks = 1:2,
                   labels = c("Newly Diagnosed", "Relapsed")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"], 3),
                            "  P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 5))) +
  theme_minimal()




##########################
## Selecting ISS stage ###
##########################


## ISS 1 vs 2 Vs 3

ISS<-read_tsv("MMRF_CoMMpass_IA19_PER_PATIENT.tsv")
ISS<-data.frame(ISS)
ISS<-ISS%>%select(PUBLIC_ID, D_PT_iss)
ISS
nrow(ISS)
ISS$PUBLIC_ID <-str_extract(ISS$PUBLIC_ID, "[0-9]{4}")
colnames(ISS)<-c("patient", "ISS")
ISS

nrow(ids)
ids1
unique(ids$patient)
ids1<-ids%>%filter(visit==1) # base line
ids1
nrow(ids1)
nrow(ISS)

ids$patient%in%ISS$patient
ISS$patient%in%ids$patient
head(ids1)
head(ISS)

ISS1<-merge(ids1, ISS, by="patient")
ISS1%>%filter(ISS == 3)
ISS1
nrow(ISS1)

ISS_1vs2<-ISS1%>%filter(ISS %in% c("1", "2", "3")) 
nrow(ISS1)

ISS_1vs2<-na.omit(ISS_1vs2)
nrow(ISS_1vs2)
ISS_1vs2$ISS <- factor(ISS_1vs2$ISS, levels = c("1", "2", "3"))
ISS_1vs2$ISS
ids

summary(colnames(counts) %in% ISS_1vs2$sample_id)
counts3 <-counts[,colnames(counts) %in% ISS_1vs2$sample_id]
counts3
counts3<- as.matrix(counts3)
ncol(counts3)
nrow(ids1)

#translocations_df2<-translocations_df2%>%filter(sample_id %in% colnames(counts2)) ###############


## Create DDS object for ISS
dds1 <- DESeqDataSetFromMatrix(countData = counts3, colData = ISS_1vs2, design = ~ ISS)   # Effect of visit, account for paired samples

## Analyze expression differences between ISS stage
dds1 <- DESeq(dds1)   # Run this analysis on Fillmore workstation
res1 <- results(dds1)
dds1
res1

## Plot MA plots for ISS
resultsNames(dds1)
resLFC1 <- lfcShrink(dds1, coef = c("ISS_2_vs_1"), type = "apeglm")   # Fold change shrinkage for visualization
plotMA(res1, ylim=c(-2,2), main="MA plot ISS 2 vs 1")
plotMA(resLFC1, ylim=c(-2,2), main = " MA plot ISS 2 vs 1 LFC shrink")

## Plot PCA for ISS
vsd1 <- vst(dds1, blind = FALSE)
plotPCA(vsd1, intgroup = c("ISS"))
temp <- as.data.frame(assay(vsd1))

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res1), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

res1$GENE <- convert_df$target
res1$GENE_NAME <- convert_df$description
res1$GENE_NAME <- sub("\\[.*", "", res1$GENE_NAME)
convert_df

ISS_1vs2

temp <- as.data.frame(assay(vsd))
temp$GENE <- res1$GENE
temp$GENE_NAME <- res1$GENE_NAME
ncol(temp)
head(temp)

temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)
colnames(temp)%in%ISS_1vs2$sample_id
rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

nrow(ISS_1vs2)
nrow(temp_flip2)
temp_bar2<-merge(temp_flip2, ISS_1vs2, by="sample_id")
head(temp_bar2)

temp_bar2$ISS <- factor(temp_bar2$ISS, levels = c("1", "2", "3"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2

stat_pvalue <- temp_bar2 %>% 
  rstatix::wilcox_test(temp_flip ~ ISS) %>%
  filter(p < 0.05) %>% 
  rstatix::add_significance("p") %>% 
  rstatix::add_y_position() %>% 
  mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

stat_pvalue
temp_bar2 %>% rstatix::wilcox_test(temp_flip ~ ISS)

temp_bar2%>%ggplot(aes(ISS, temp_flip)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "Patient Status", 
                   breaks = 1:3,
                   labels = c("ISS 1","ISS 2", "ISS 3")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "CD70 expression based on ISS stage")+
  ggpubr::stat_pvalue_manual(stat_pvalue, label = "p.signif") +
  theme_minimal()


### MUTATIONS #####

translocations_df<-read_tsv("SeqFISH Files_MMRF_CoMMpass_IA19_genome_tumor_only_mm_igtx_pairoscope.tsv")
translocations_df<-data.frame(translocations_df)

translocations_df1<-data.frame(sample_id = translocations_df$SAMPLE,
                               patient = str_extract(translocations_df$SAMPLE,"[0-9]{4}"),
                               visit = str_extract(translocations_df$SAMPLE, "_[0-9]_"),
                               type = str_extract(translocations_df$SAMPLE, "_[A-Z]{2}_"),
                               NDS2 = translocations_df$NSD2_CALL,
                               CCND3 = translocations_df$CCND3_CALL,
                               MYC = translocations_df$MYC_CALL,
                               MAFA = translocations_df$MAFA_CALL,
                               CCND1 = translocations_df$CCND1_CALL,
                               CCN2 = translocations_df$CCND2_CALL,
                               MAF = translocations_df$MAF_CALL,
                               MAFB = translocations_df$MAFB_CALL)

translocations_df1  
translocations_df1$visit<-str_extract(translocations_df1$visit, "[0-9]")
translocations_df1$type<-str_extract(translocations_df1$type, "[A-Z]+")
translocations_df1<-translocations_df1%>%filter(type == "BM")  
nrow(translocations_df1)



loci414<-data_frame(sample_id = translocations_df$SAMPLE,
                    patient = str_extract(translocations_df$SAMPLE,"[0-9]{4}"),
                    visit = str_extract(translocations_df$SAMPLE, "_[0-9]_"),
                    type = str_extract(translocations_df$SAMPLE, "_[A-Z]{2}_"),
                    NDS2 = translocations_df$NSD2_CALL,
                    NDS2_Ig = translocations_df$NSD2_IGSOURCE)

loci414
loci414$visit<-str_extract(loci414$visit, "[0-9]")
loci414$type<-str_extract(loci414$type, "[A-Z]+")
loci414<-data.frame(loci414)
loci414<-loci414%>%filter(type == "BM")
nrow(loci414)

seqFISH<-read_tsv("SeqFISH Files_MMRF_CoMMpass_IA19_genome_gatk_cna_seqFISH.tsv")
seqFISH<-data.frame(seqFISH)

seqFISH<-data.frame(sample_id = seqFISH$SAMPLE,
                    patient = str_extract(seqFISH$SAMPLE,"[0-9]{4}"),
                    visit = str_extract(seqFISH$SAMPLE, "_[0-9]_"),
                    type = str_extract(seqFISH$SAMPLE, "_[A-Z]{2}_"),
                    "del17p13" = seqFISH$SeqWGS_Cp_17p13_20percent,
                    "TP53" = seqFISH$SeqWGS_Cp_TP53_20percent,
                    "Gain1q21" = seqFISH$SeqWGS_Cp_1q21_20percent,
                    "del1p22" = seqFISH$SeqWGS_Cp_1p22_20percent,
                    "MYC" = seqFISH$SeqWGS_Cp_MYC_20percent)

seqFISH
seqFISH$visit <- str_extract(seqFISH$visit, "[0-9]")
seqFISH$type = str_extract(seqFISH$type, "[A-Z]+")
seqFISH<-seqFISH%>%filter(type == "BM")
seqFISH
nrow(seqFISH)

#####################################
########## NDS2 or t(4:14) ##########
#####################################

## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts


## Filter for translocation
translocations_df1
counts <- counts[, colnames(counts) %in% translocations_df1$sample_id]
counts

## Transform every sample as factor for NDS2
translocations_df2<-translocations_df1
translocations_df2$NDS2 <- factor(translocations_df2$NDS2, levels = c("0", "1"))
translocations_df2$sample_id

colnames(counts) %in% translocations_df2$sample_id

counts2 <- counts[, colnames(counts) %in% translocations_df2$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(translocations_df2)

counts2<- as.matrix(counts2)
ncol(counts2)
counts2
nrow(translocations_df2)

translocations_df2<-translocations_df2%>%filter(sample_id %in% colnames(counts2))
translocations_df2



## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = translocations_df2, design = ~ NDS2)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "NDS2_1_vs_0", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink NDS2.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot NDS2 ")
plotMA(resLFC, ylim=c(-2,2), main = "MA plot NDS2 shrink")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("NDS2"))
res

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df
res

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))

temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)



# Select names by extracting number of the sample
temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)

rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

head(translocations_df2)
temp_bar2<-merge(temp_flip2, translocations_df2, by="sample_id")
head(temp_bar2)
nrow(temp_bar2)

temp_bar2$NDS2 <- factor(temp_bar2$NDS2, levels = c("0", "1"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2<-data.frame(temp_bar2)
nrow(temp_bar2)
head(temp_bar2)

class(temp_flip2)

# PLOT variance of PAX5 expression across every sample
temp_bar2%>%ggplot(aes( x = NDS2, y = temp_flip))+ 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "t(4:14) Status",
                   breaks = 0:1,
                   labels = c("Negative", "Positive")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5 expression based on t(4:14) [NDS2] status", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"], 3),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 3))) +
  theme_minimal()



##################################
#########CCND1 t(11:14)###########
##################################

## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts


## Filter for translocation
translocations_df1
counts <- counts[, colnames(counts) %in% translocations_df1$sample_id]
counts

## Transform every sample as factor for NDS2
translocations_df2<-translocations_df1
head(translocations_df2)
translocations_df2$CCND1 <- factor(translocations_df2$CCND1, levels = c("0", "1"))

counts2 <- counts[, colnames(counts) %in% translocations_df2$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(translocations_df2)

counts2<- as.matrix(counts2)
ncol(counts2)
counts2
nrow(translocations_df2)

translocations_df2<-translocations_df2%>%filter(sample_id %in% colnames(counts2))
translocations_df2
nrow(translocations_df2)


## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = translocations_df2, design = ~ CCND1)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "CCND1_1_vs_0", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink CCND1.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot CCND1 ")
plotMA(resLFC, ylim=c(-2,2), main = "MA plot CCND1 shrink")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("CCND1"))
res

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))

temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)



# Select names by extracting number of the sample
temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)

rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

head(translocations_df2)
temp_bar2<-merge(temp_flip2, translocations_df2, by="sample_id")
head(temp_bar2)
nrow(temp_bar2)

temp_bar2$CCND1 <- factor(temp_bar2$CCND1, levels = c("0", "1"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2

# PLOT variance of PAX5 expression across every sample
temp_bar2%>%ggplot(aes( x = CCND1, y = temp_flip)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "t(11:14) Status",
                   breaks = 0:1,
                   labels = c("Negative", "Positive")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5 expression based on t(11:14) [CCND1] status", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"],3),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"],3))) +
  theme_minimal()

options(scipen = 100)
res$pvalue[res$GENE == "PAX5"]

################################
#######CCND3 or t(6:14) #########
################################


## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts


## Filter for translocation
translocations_df1
counts <- counts[, colnames(counts) %in% translocations_df1$sample_id]
counts

## Transform every sample as factor for NDS2
translocations_df2<-translocations_df1
head(translocations_df2)
translocations_df2$CCND3 <- factor(translocations_df2$CCND3, levels = c("0", "1"))

counts2 <- counts[, colnames(counts) %in% translocations_df2$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(translocations_df2)

counts2<- as.matrix(counts2)
ncol(counts2)
counts2
nrow(translocations_df2)

translocations_df2<-translocations_df2%>%filter(sample_id %in% colnames(counts2))
translocations_df2
nrow(translocations_df2)


## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = translocations_df2, design = ~ CCND3)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "CCND3_1_vs_0", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink CCND3.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot CCND3 ")
plotMA(resLFC, ylim=c(-2,2), main = "MA plot CCND3 shrink")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("CCND3"))
res

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))

temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)



# Select names by extracting number of the sample
temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)

rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

head(translocations_df2)
temp_bar2<-merge(temp_flip2, translocations_df2, by="sample_id")
head(temp_bar2)
nrow(temp_bar2)

temp_bar2$CCND3 <- factor(temp_bar2$CCND3, levels = c("0", "1"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2


# PLOT variance of PAX5 expression across every sample
temp_bar2%>%ggplot(aes( x = CCND3, y = temp_flip)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "t(6:14) Status",
                   breaks = 0:1,
                   labels = c("Negative", "Positive")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5 expression based on t(6:14) [CCND3] status", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"], 3),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 3))) +
  theme_minimal()

res$pvalue[res$GENE %in% "PAX5"]



############################
########### MYC ############
############################

## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts


## Filter for translocation
translocations_df1
counts <- counts[, colnames(counts) %in% translocations_df1$sample_id]
counts

## Transform every sample as factor for NDS2
translocations_df2<-translocations_df1
head(translocations_df2)
translocations_df2$MYC <- factor(translocations_df2$MYC, levels = c("0", "1"))

counts2 <- counts[, colnames(counts) %in% translocations_df2$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(translocations_df2)

counts2<- as.matrix(counts2)
ncol(counts2)
counts2
nrow(translocations_df2)

translocations_df2<-translocations_df2%>%filter(sample_id %in% colnames(counts2))
translocations_df2
nrow(translocations_df2)


## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = translocations_df2, design = ~ MYC)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "MYC_1_vs_0", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink MYC.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot MYC ")
plotMA(resLFC, ylim=c(-2,2), main = "MA plot MYC shrink")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("MYC"))
res

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))

temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)



# Select names by extracting number of the sample
temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)

rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

head(translocations_df2)
temp_bar2<-merge(temp_flip2, translocations_df2, by="sample_id")
head(temp_bar2)
nrow(temp_bar2)

temp_bar2$MYC <- factor(temp_bar2$MYC, levels = c("0", "1"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2


# PLOT variance of PAX5 expression across every sample
temp_bar2%>%ggplot(aes( x = MYC, y = temp_flip)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "MYC mutation",
                   breaks = 0:1,
                   labels = c("Negative", "Positive")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5 expression based on MYC mutation", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"], 3),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 3))) +
  theme_minimal()

res$pvalue[res$GENE == "PAX5"]


####################################
######### c-MAF or t(14:16)#########
#####################################



## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts


## Filter for translocation
translocations_df1
counts <- counts[, colnames(counts) %in% translocations_df1$sample_id]
counts

## Transform every sample as factor for NDS2
translocations_df2<-translocations_df1
head(translocations_df2)
translocations_df2$MAF <- factor(translocations_df2$MAF, levels = c("0", "1"))

counts2 <- counts[, colnames(counts) %in% translocations_df2$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(translocations_df2)

counts2<- as.matrix(counts2)
ncol(counts2)
counts2
nrow(translocations_df2)

translocations_df2<-translocations_df2%>%filter(sample_id %in% colnames(counts2))
translocations_df2
nrow(translocations_df2)


## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = translocations_df2, design = ~ MAF)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "MAF_1_vs_0", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink MAF.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot MAF ")
plotMA(resLFC, ylim=c(-2,2), main = "MA plot MAF shrink")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("MAF"))
res

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

convert_df
res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))

temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)



# Select names by extracting number of the sample
temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)

rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

head(translocations_df2)
temp_bar2<-merge(temp_flip2, translocations_df2, by="sample_id")
head(temp_bar2)
nrow(temp_bar2)

temp_bar2$MAF <- factor(temp_bar2$MAF, levels = c("0", "1"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2


# PLOT variance of PAX5 expression across every sample
temp_bar2%>%ggplot(aes( x = MAF, y = temp_flip)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "t(14;16) Status",
                   breaks = 0:1,
                   labels = c("Negative", "Positive")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5 expression based on t(14;16) [c-MAF] status", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"], 3),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 3))) +
  theme_minimal()

res$pvalue[res$GENE == "PAX5"]



###########################################
############# MAFB or t(14:20)#############
###########################################


## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts


## Filter for translocation
translocations_df1
counts <- counts[, colnames(counts) %in% translocations_df1$sample_id]
counts

## Transform every sample as factor for NDS2
translocations_df2<-translocations_df1
head(translocations_df2)
translocations_df2$MAFB <- factor(translocations_df2$MAFB, levels = c("0", "1"))

counts2 <- counts[, colnames(counts) %in% translocations_df2$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(translocations_df2)

counts2<- as.matrix(counts2)
ncol(counts2)
counts2
nrow(translocations_df2)

translocations_df2<-translocations_df2%>%filter(sample_id %in% colnames(counts2))
translocations_df2
nrow(translocations_df2)


## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = translocations_df2, design = ~ MAFB)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "MAFB_1_vs_0", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink MAFB.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot MAFB ")
plotMA(resLFC, ylim=c(-2,2), main = "MA plot MAFB shrink")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("MAFB"))
res

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))

temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)



# Select names by extracting number of the sample
temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)

rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

head(translocations_df2)
temp_bar2<-merge(temp_flip2, translocations_df2, by="sample_id")
head(temp_bar2)
nrow(temp_bar2)

temp_bar2$MAFB <- factor(temp_bar2$MAFB, levels = c("0", "1"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2

# PLOT variance of PAX5 expression across every sample
temp_bar2%>%ggplot(aes( x = MAFB, y = temp_flip)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "t(14;20) Status",
                   breaks = 0:1,
                   labels = c("Negative", "Positive")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5 expression based on t(14;20) [MAFB] status", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"], 3),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 3))) +
  theme_minimal()

res$pvalue[res$GENE == "PAX5"]




#####################################
############# del17p13###############
#####################################

seqFISH

## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts


## Filter for translocation
head(seqFISH)
counts <- counts[, colnames(counts) %in% seqFISH$sample_id]
counts

## Transform every sample as factor for NDS2
seqFISH1<-seqFISH
seqFISH1
seqFISH1$del17p13 <- factor(seqFISH1$del17p13, levels = c("0", "1"))
seqFISH1$del17p13

colnames(counts) %in% seqFISH1$sample_id

counts2 <- counts[, colnames(counts) %in% seqFISH1$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(translocations_df2)

counts2<- as.matrix(counts2)
ncol(counts2)
counts2
nrow(seqFISH1)

seqFISH1<-seqFISH1%>%filter(sample_id %in% colnames(counts2))
seqFISH1

ncol(counts2)
nrow(seqFISH1)


## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = seqFISH1, design = ~ del17p13)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "del17p13_1_vs_0", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink del17p13.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot del17p13 ")
plotMA(resLFC, ylim=c(-2,2), main = "MA plot del17p13 shrink")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("del17p13"))
res

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))

temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)



# Select names by extracting number of the sample
temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)

rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

head(seqFISH1)
temp_bar2<-merge(temp_flip2, seqFISH1, by="sample_id")
head(temp_bar2)
nrow(temp_bar2)

temp_bar2$del17p13 <- factor(temp_bar2$del17p13, levels = c("0", "1"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2


# PLOT variance of PAX5 expression across every sample
temp_bar2%>%ggplot(aes( x = del17p13, y = temp_flip)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "del17p13 Status",
                   breaks = 0:1,
                   labels = c("Negative", "Positive")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5 expression based on del17p13 status", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"], 3),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 3))) +
  theme_minimal()

res$pvalue[res$GENE == "PAX5"]


#####################################
#############  TP53    ##############
#####################################

head(seqFISH)

## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts


## Filter for translocation
head(seqFISH)
counts <- counts[, colnames(counts) %in% seqFISH$sample_id]
counts

## Transform every sample as factor for NDS2
seqFISH1<-seqFISH
seqFISH1
seqFISH1$TP53 <- factor(seqFISH1$TP53, levels = c("0", "1"))
seqFISH1$TP53

colnames(counts) %in% seqFISH1$sample_id

counts2 <- counts[, colnames(counts) %in% seqFISH1$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(seqFISH1)

counts2<- as.matrix(counts2)
ncol(counts2)
counts2
nrow(seqFISH1)

seqFISH1<-seqFISH1%>%filter(sample_id %in% colnames(counts2))
seqFISH1

ncol(counts2)
nrow(seqFISH1)


## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = seqFISH1, design = ~ TP53)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "TP53_1_vs_0", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink TP53.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot TP53 ")
plotMA(resLFC, ylim=c(-2,2), main = "MA plot TP53 shrink")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("TP53"))
res

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

convert_df
res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))

temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)



# Select names by extracting number of the sample
temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)

rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

head(seqFISH1)
temp_bar2<-merge(temp_flip2, seqFISH1, by="sample_id")
head(temp_bar2)
nrow(temp_bar2)

temp_bar2$TP53 <- factor(temp_bar2$TP53, levels = c("0", "1"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2


# PLOT variance of PAX5 expression across every sample
temp_bar2%>%ggplot(aes( x = TP53, y = temp_flip)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "TP53 Status",
                   breaks = 0:1,
                   labels = c("Negative", "Positive")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5 expression based on TP53 status", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"], 3),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 3))) +
  theme_minimal()

res$pvalue[res$GENE == "PAX5"]




#####################################
#############  Gain1q21   ###########
#####################################

head(seqFISH)

## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts


## Filter for translocation
head(seqFISH)
counts <- counts[, colnames(counts) %in% seqFISH$sample_id]
counts

## Transform every sample as factor for NDS2
seqFISH1<-seqFISH
seqFISH1
seqFISH1$Gain1q21 <- factor(seqFISH1$Gain1q21, levels = c("0", "1"))
seqFISH1$Gain1q21

colnames(counts) %in% seqFISH1$sample_id

counts2 <- counts[, colnames(counts) %in% seqFISH1$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(seqFISH1)

counts2<- as.matrix(counts2)
ncol(counts2)
counts2
nrow(seqFISH1)

seqFISH1<-seqFISH1%>%filter(sample_id %in% colnames(counts2))
seqFISH1

ncol(counts2)
nrow(seqFISH1)


## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = seqFISH1, design = ~ Gain1q21)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "Gain1q21_1_vs_0", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink Gain1q21.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot Gain1q21 ")
plotMA(resLFC, ylim=c(-2,2), main = "MA plot Gain1q21 shrink")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("Gain1q21"))
res

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

convert_df
res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))

temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)



# Select names by extracting number of the sample
temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)

rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

head(seqFISH1)
temp_bar2<-merge(temp_flip2, seqFISH1, by="sample_id")
head(temp_bar2)
nrow(temp_bar2)

temp_bar2$Gain1q21 <- factor(temp_bar2$Gain1q21, levels = c("0", "1"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2
nrow(temp_bar2)


# PLOT variance of PAX5 expression across every sample
temp_bar2%>%ggplot(aes( x = Gain1q21, y = temp_flip)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "Gain1q21 Status",
                   breaks = 0:1,
                   labels = c("Negative", "Positive")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5 expression based on Gain1q21 status", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"], 3),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 3))) +
  theme_minimal()

res$pvalue[res$GENE == "PAX5"]


####################################
############ del1p22 ###############
###################################

head(seqFISH)

## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts


## Filter for translocation
head(seqFISH)
counts <- counts[, colnames(counts) %in% seqFISH$sample_id]
counts

## Transform every sample as factor for NDS2
seqFISH1<-seqFISH
seqFISH1
seqFISH1$del1p22 <- factor(seqFISH1$del1p22, levels = c("0", "1"))
seqFISH1$del1p22

colnames(counts) %in% seqFISH1$sample_id

counts2 <- counts[, colnames(counts) %in% seqFISH1$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(seqFISH1)

counts2<- as.matrix(counts2)
ncol(counts2)
counts2
nrow(seqFISH1)

seqFISH1<-seqFISH1%>%filter(sample_id %in% colnames(counts2))
seqFISH1

ncol(counts2)
nrow(seqFISH1)


## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = seqFISH1, design = ~ del1p22)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "del1p22_1_vs_0", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink del1p22.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot del1p22 ")
plotMA(resLFC, ylim=c(-2,2), main = "MA plot del1p22 shrink")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("del1p22"))
res

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

convert_df
res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))

temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)



# Select names by extracting number of the sample
temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)

rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

head(seqFISH1)
temp_bar2<-merge(temp_flip2, seqFISH1, by="sample_id")
head(temp_bar2)
nrow(temp_bar2)

temp_bar2$del1p22 <- factor(temp_bar2$del1p22, levels = c("0", "1"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2


# PLOT variance of PAX5 expression across every sample
temp_bar2%>%ggplot(aes( x = del1p22, y = temp_flip)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "del1p22 Status",
                   breaks = 0:1,
                   labels = c("Negative", "Positive")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5 expression based on del1p22 status", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"], 3),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 3))) +
  theme_minimal()

res$pvalue[res$GENE == "PAX5"]


###############################################
########## MYC in seqFISH #####################
###############################################

head(seqFISH)

## Read HtSeq data
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
counts <- as.matrix(raw[-1])
counts
rownames(counts) <- raw$Gene
counts


## Filter for translocation
head(seqFISH)
counts <- counts[, colnames(counts) %in% seqFISH$sample_id]
counts

## Transform every sample as factor for NDS2
seqFISH1<-seqFISH
seqFISH1
seqFISH1$MYC <- factor(seqFISH1$MYC, levels = c("0", "1"))
seqFISH1$MYC

colnames(counts) %in% seqFISH1$sample_id

counts2 <- counts[, colnames(counts) %in% seqFISH1$sample_id]
counts2
nrow(counts2)
ncol(counts2)
nrow(seqFISH1)

counts2<- as.matrix(counts2)
ncol(counts2)
counts2
nrow(seqFISH1)

seqFISH1<-seqFISH1%>%filter(sample_id %in% colnames(counts2))
seqFISH1

ncol(counts2)
nrow(seqFISH1)


## Create DDS object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = seqFISH1, design = ~ MYC)   # Effect of visit, account for paired samples

## Analyze expression differences between 1st relapse vs newly diagnosed
dds <- DESeq(dds)   # Run this analysis on Fillmore workstation
res <- results(dds)
dds
res

## Plot MA plots 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "MYC_1_vs_0", type = "apeglm")   # Fold change shrinkage for visualization

pdf(" MA plot shrink MYC in seqFISH dataset.pdf")
plotMA(res, ylim=c(-2,2), main = "MA plot MYC in seqFISH dataset ")
plotMA(resLFC, ylim=c(-2,2), main = "MA plot MYC in seqFISH dataset shrink")
dev.off()

## Plot PCA 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("MYC"))
res

## Convert ENSG to gene name
convert_df <- gconvert(rownames(res), organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = F)

convert_df
res$GENE <- convert_df$target
res$GENE_NAME <- convert_df$description
res$GENE_NAME <- sub("\\[.*", "", res$GENE_NAME)
convert_df

## Output flat file for data visualization
temp <- as.data.frame(assay(vsd))

temp$GENE <- res$GENE
temp$GENE_NAME <- res$GENE_NAME
ncol(temp)
head(temp)



# Select names by extracting number of the sample
temp_flip<-temp%>%filter(GENE == "PAX5")
head(temp_flip)

rownames(temp_flip)<-NULL
temp_flip<-t(temp_flip)

temp_flip2<-data.frame(temp_flip)
temp_flip2
temp_flip2$sample_id<-rownames(temp_flip)

rownames(temp_flip2)<-NULL
head(temp_flip2)

head(seqFISH1)
temp_bar2<-merge(temp_flip2, seqFISH1, by="sample_id")
head(temp_bar2)
nrow(temp_bar2)

temp_bar2$MYC <- factor(temp_bar2$MYC, levels = c("0", "1"))
temp_bar2$temp_flip<-as.numeric(temp_bar2$temp_flip)

temp_bar2


# PLOT variance of PAX5 expression across every sample
temp_bar2%>%ggplot(aes( x = MYC, y = temp_flip)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2),
              alpha = 0.4,
              color = "red") +
  scale_x_discrete(name = "MYC Status",
                   breaks = 0:1,
                   labels = c("Negative", "Positive")) +
  scale_y_continuous(name = expression("Normalized log "[2]*" (Counts)")) +
  ggtitle(label = "PAX5 expression based on MYC status (SeqFISH)", 
          subtitle = paste0("log2FC = ",
                            round(res$log2FoldChange[res$GENE %in% "PAX5"], 3),
                            "   P-value = ", round(res$pvalue[res$GENE %in% "PAX5"], 3))) +
  theme_minimal()

res$pvalue[res$GENE == "PAX5"]
