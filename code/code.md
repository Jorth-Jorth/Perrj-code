# code
```R
############################################################################===
############################################################################===
rm(list=ls())
options(stringsAsFactors = F)


library(tidyverse)
library(GEOquery)

#### 1 three microarray datasets and one rna-seq_count ----
###  1.1 GSE31821 ----
sm_filename <-  "batch/GSE19804/GSE19804_series_matrix.txt.gz"
GSE19804_sm <- getGEO(filename = sm_filename, getGPL = F) 
GSE19804_pd <- pData(GSE19804_sm) 
GSE19804_targets <- GSE19804_pd %>%
  dplyr::select(2, 1) %>%#取第2列和第1列
  separate(title, into = paste0("x", 1:3), sep = " ") %>%#把title按照空格分割成为三个部分
  mutate(sample_id = geo_accession,#添加一列
         file_name = paste0(sample_id, ".CEL.gz"),
         sample_name = x3,
         patient_id = str_sub(x3, 1, -2),#str_sub----Extract and replace substrings from a character vector.从第一个到倒数第二个字符
         GSE_num = "GSE19804",
         group = paste0(x2, 1)) %>%
  dplyr::select(sample_id:group)#选择这些列

table(GSE19804_targets$group)



###  1.2 GSE115574----
sm_filename <- "batch/GSE18842/GSE18842_series_matrix.txt.gz"
GSE18842_sm <- getGEO(filename = sm_filename, getGPL = F) 
GSE18842_pd <- pData(GSE18842_sm) 
GSE18842_targets <- GSE18842_pd %>%
  dplyr::select(2, 1) %>%#选中该第二列和第一列
  separate(title, into = paste0("x", 1:2), sep = " ", remove = F) %>%#按照空格分成两列
  mutate(sample_id = geo_accession,#添加几列
         file_name = paste0(sample_id, ".CEL.gz"),
         sample_name = str_c(x1, x2, sep = "_"),# Join multiple strings into a single string.
         patient_id = x2,
         GSE_num = "GSE18842",
         group = ifelse(x1 == "Tumor", "Cancer2", "Normal2")) %>%#如果X1等于tumor，那么填cancer2，否则填normal2
  dplyr::select(sample_id:group)

table(GSE18842_targets$group)

###  1.3 GSE41177----
GSE18842_sm <- getGEO(filename = sm_filename, getGPL = F) 
GSE18842_pd <- pData(GSE18842_sm) 
GSE18842_targets <- GSE18842_pd %>%
  dplyr::select(2, 1) %>%#选中该第二列和第一列
  separate(title, into = paste0("x", 1:2), sep = " ", remove = F) %>%#按照空格分成两列
  mutate(sample_id = geo_accession,#添加几列
         file_name = paste0(sample_id, ".CEL.gz"),
         sample_name = str_c(x1, x2, sep = "_"),# Join multiple strings into a single string.
         patient_id = x2,
         GSE_num = "GSE18842",
         group = ifelse(x1 == "Tumor", "Cancer2", "Normal2")) %>%#如果X1等于tumor，那么填cancer2，否则填normal2
  dplyr::select(sample_id:group)

table(GSE18842_targets$group)

###  1.4 RNA-seq_count ----
library(data.table)
count =fread('RNA-seq_count',header = T,qute = F)
sample=colnames(count)
group=c(rep('P4',6),rep('P14',7))
count_target=data.frame(sample=sample,
              group=group)


#### 2 combined and imported ----
GSE2_targets <- bind_rows(GSE19804_targets, GSE18842_targets)
rownames(GSE2_targets) <- GSE2_targets$sample_id


library(affy)
GSE2_cel <- ReadAffy(filenames = GSE2_targets$file_name,
                     celfile.path = "batch/GSE2combined/",
                     phenoData = GSE2_targets)

GSE2_targets$batch_date <- GSE2_cel@protocolData@data$ScanDate %>% 
  str_sub(1, 8) %>% #取1到8
  as.factor()#因子化，这样GSE2_targets就多了一列batch_date

#### 3 qc ----
###  3.1 arrayQualityMetrics ----
 library(arrayQualityMetrics)
 arrayQualityMetrics(CEL2, 
                     outdir = "batch/GSE2combined/GSE2combined_QC_affybatch", 
                    do.logtransform = TRUE,
                     force = TRUE,
                     intgroup = "group")


###  3.2 RNA degradation ----
library(simpleaffy)
RNAdeg <- AffyRNAdeg(GSE2_cel)
RNAd <- summaryAffyRNAdeg(RNAdeg)
cols <- GSE2_targets$batch_date

plotAffyRNAdeg(RNAdeg, cols = cols)
box()

###  3.3 sampleaffy QC ----
library(simpleaffy)
GSE2_qc <- simpleaffy::qc(GSE2_cel)
pdf("batch/GSE2_qc.pdf", height = 100, width = 10)
plot(GSE2_qc)
dev.off()

###  3.4 boxplot ----
boxplot(GSE2_cel, 
        las = 2, 
        outline = F,
        col = as.factor(GSE2_targets$group))

boxplot(count, 
        las = 2, 
        outline = F,
        col = as.factor(count_target$group))
        
#### 4 normalized----
GSE2_rma <- rma(GSE2_cel)

count_rlog=rlog(as.matrix(count))

## boxplot
boxplot(exprs(GSE2_rma), 
        las = 2, 
        outline = F,
        col = as.factor(GSE2_targets$group))
        
boxplot(exprs(count_rlog), 
        las = 2, 
        outline = F,
        col = as.factor(count_target$group))

## PCA
PCA_new(exprs(GSE2_rma), 
        ntop = nrow(exprs(GSE2_rma)),
        group = GSE2_targets$group)

## hclust
expr <- exprs(GSE2_rma)
colnames(expr) <- paste(GSE2_targets$group, 1:211, sep = "_")

dist <- dist(t(expr))
hc = hclust(dist)

install.packages('factoextra')
library(factoextra)
library(RColorBrewer)
pdf("batch/clust.pdf", height = 10, width = 30)
fviz_dend(hc, 
          k = 6, 
          cex = 0.5,
          k_colors = brewer.pal(6, "Dark2"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())
dev.off()

#### 5 Remove batch effects ----
batch <- GSE2_targets$GSE_num
library(sva)
mod <-  model.matrix(~ as.factor(group), data= GSE2_targets)
GSE2_combat <- ComBat(dat = expr,
                       batch = batch,
                       par.prior = TRUE,
                       prior.plots = FALSE)

## boxplot
boxplot(GSE2_combat, 
        las = 2, 
        outline = F,
        col = as.factor(GSE2_targets$group))

## PCA
PCA_new(GSE2_combat, 
        ntop = nrow(GSE2_combat),
        group = GSE2_targets$group)

## hclust
expr <- GSE2_combat
colnames(expr) <- paste(GSE2_targets$group, 1:211, sep = "_")
dist <- dist(t(expr))
hc = hclust(dist)
library(factoextra)
library(RColorBrewer)
pdf("batch/clust2.pdf", height = 10, width = 30)
fviz_dend(hc, k = 6, 
          cex = 0.5,
          k_colors = brewer.pal(6, "Set2"),
          color_labels_by_k = TRUE,
          ggtheme = theme_classic())
dev.off()

#### 6 annotation ----
#BiocManager::install('hgu133plus2.db')
library(hgu133plus2.db)

ids_GPL570=toTable(hgu133plus2SYMBOL)
GSE2_combat$probe_id=as.character(rownames(GSE2_combat)) 
exprSet_GSE2_combat <- merge(ids_GPL570,GSE2_combat, by = "probe_id") 
table(is.na(exprSet_GSE2_combat))
exprSet_GSE2_combat<- exprSet_GSE2_combat[order(rowMeans(exprSet_GSE2_combat[,3:ncol(exprSet_GSE2_combat)]), 
                                           decreasing = TRUE), ]
exprSet_GSE2_combat <- exprSet_GSE2_combat[!duplicated(exprSet_GSE2_combat$symbol), ]
rownames(exprSet_GSE2_combat) <- exprSet_GSE2_combat$symbol
exprSet_GSE2_combat <- exprSet_GSE2_combat[, -c(1, 2)]


#### 6 DEG ----
## limma ##
library(limma)
data= exprSet_GSE2_combat
gl=GSE2_targets
exp="AF"
ctr="SR"


design <- model.matrix(~0+factor(gl$group))
colnames(design)=levels(factor(gl$group))
rownames(design)=colnames(data)
design

contrast.matrix <- makeContrasts(contrasts=paste0(exp,'-',ctr),  #"exp/ctrl"
                                 levels = design)
contrast.matrix 
##step1
fit <- lmFit(data,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
head(nrDEG)
write.csv(nrDEG,"./AF_limma_DEG.csv",quote = F)

## DEseq2 ##

library(DESeq2)
library("BiocParallel")

## build dds DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count,
                                  colData = count_target,
                                  design = ~ group)

dds$group_list <- relevel(dds$group_list, ref = ctr)  

keep <- rowSums(counts(dds)) >= 1.5*ncol(counts) 
dds <- dds[keep,] 
dds <- DESeq(dds,quiet = F) 
res <- results(dds,contrast=c("group_list", exp, ctr)) 
resOrdered <- res[order(res$padj),]
tempDEG <- as.data.frame(resOrdered)
DEG_DEseq2 <- na.omit(tempDEG)
head(DEG_DEseq2)
write.csv(EG_DEseq2,"./senescence_DEseq2_DEG.csv",quote = F)

#### 7 GSEA ----
count_SERPIEN1=count_rlog[,1:7]
cutoff=sum(count_SERPINE1['SERPINE1',])/2
count_SERPIEN1=as.data.frame(t(count_SERPIEN1))
count_SERPIEN1$change = as.factor(ifelse(count_SERPIEN1$SERPINE1 >cutoff  ,'High','Low'))

## limma ##
library(limma)
count_SERPIEN1=as.data.frame(t(count_SERPIEN1))
data= count_SERPIEN1
gl=count_SERPIEN1$change
exp="High"
ctr="Low"


design <- model.matrix(~0+factor(gl))
colnames(design)=levels(factor(gl))
rownames(design)=colnames(data)
design

contrast.matrix <- makeContrasts(contrasts=paste0(exp,'-',ctr),  #"exp/ctrl"
                                 levels = design)
contrast.matrix 
##step1
fit <- lmFit(data,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG_SERPINE1 = na.omit(tempOutput) 
head(nrDEG_SERPINE1)


library(GSEABase)
geneList=data.frame(gene=rownames(nrDEG_SERPINE1),
                    logFC=nrDEG_SERPINE1$logFC)
gmts <- list.files(pattern = 'c2.all.v2023.1.Hs.symbols.gmt')
lapply(gmts, function(gmtfile){
  # gmtfile=gmts[2]
  geneset <- read.gmt(file.path(d,gmts)) 
  print(paste0('Now process the ',gmtfile))
  egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE)
  head(egmt)
  # gseaplot(egmt, geneSetID = rownames(egmt[1,]))  
  return(egmt)
})



#### END ----


```