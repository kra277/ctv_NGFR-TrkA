## code to prepare `DATASET`

library(tidyverse)
library(DESeq2)
library(EnrichmentBrowser)
library(RUVSeq)
library(sva)
library(corrplot)

# Function for clean deseq results

clean_deseq_res <- function(deseq_res) {
  ## To add q-values, annotate, arrange, and clean the deseq2 results

  deg_res_dat <-
    as.data.frame(deseq_res) %>% ## convert into a dataframe
    tibble::rownames_to_column("ENTREZID") %>% ## convert rownames to column
    dplyr::mutate(qvalue = qvalue::qvalue(p=deseq_res$pvalue)$qvalues) ## Get qvalues

  library(Homo.sapiens)

  ## Get annotations
  anno <- AnnotationDbi::select(Homo.sapiens, keytype='ENTREZID',
                                keys=deg_res_dat$ENTREZID,
                                columns=c('SYMBOL', 'GENENAME'),
                                multiVals="first")
  ## Clean and arrange
  clean_res_df <-
    deg_res_dat %>%
    merge(anno) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::select(ENTREZID, SYMBOL, log2FoldChange,
                  pvalue, padj, qvalue, everything())

  ## Return the clean DEG dataframe
  return(clean_res_df)
}


# Get the phenotype file
colData <-
  read.csv("~/Documents/macbook_drylab/alab_packages/sample_data/phenotype.csv") %>%
  column_to_rownames("SampleID") %>%
  mutate_at(vars("Group"), as.factor)


# Get the genecount file
countData <-
  read.csv("~/Documents/macbook_drylab/alab_packages/sample_data/genecounts.csv") %>%
  column_to_rownames("gene_id") %>%
  dplyr::select(rownames(colData)) %>%
  as.matrix()

# Check if the pheno and genecounts match
all(rownames(colData) %in% colnames(countData))

# Perfrom Deseq2
ddsMat <- DESeqDataSetFromMatrix(countData = countData,
                                 colData = colData,
                                 design = ~ Group)

# Specifying the reference group
ddsMat$Group <- relevel(ddsMat$Group, ref = "Control")

# Filter
ddsMat_Filtered <-
  rowSums(counts(ddsMat) >= 10) >= round(ncol(ddsMat)*0.9)

dds <- ddsMat[ddsMat_Filtered,]

# annotating the dds
id_anno <- idMap(dds, org = "hsa", from = "ENSEMBL", to = "ENTREZID")

# compute the DEG
sample_deseq_data <- DESeq(id_anno, parallel = T)

# Get the DEG results
sample_deseq_res <- results(sample_deseq_data)

# Get clean DEG results
sample_clean_deg <- clean_deseq_res(sample_deseq_res)

# Get normalized (VST) counts
sample_vst_dd <- vst(sample_deseq_data , blind=FALSE)


# Perform RUV batch corrections for Hidden batch effects
set <- newSeqExpressionSet(counts(sample_deseq_data))
idx  <- rowSums(counts(set) > 10) >= 1
set  <- set[idx, ]
set <- betweenLaneNormalization(set, which="upper")
not.sig <- rownames(sample_clean_deg)[which(sample_clean_deg$pvalue > .1)]
empirical <- rownames(set)[ rownames(set) %in% not.sig ]
sample_set <- RUVg(set, empirical, k=2)


# Generating the correlation matrix

bc_corr_plot(bc_method = "ruv", var_int = "Group", font_size = 1)


# Add the surrogate variables to the design
ddsruv <- sample_deseq_data
ddsruv$W1 <- set$W_1
ddsruv$W2 <- set$W_2
design(ddsruv) <- ~ W1 + W2 + Group

# Re-perform Deseq with the ruv variables
de_ruv_dat <- DESeq(ddsruv, parallel = T)


# Perform SVA batch corrections for Hidden batch effects
dat  <- counts(sample_deseq_data, normalized = TRUE)

idx  <- rowMeans(dat) > 10
dat  <- dat[idx, ]

# Model matrix with all Co-variates and Variable of Interest
mod  <- model.matrix(~ Group, data=colData(sample_deseq_data))

# Null model matrix with all adjustment variables (Co-variates)
# Note: This doesnt have the variable of Interest
mod0 <- model.matrix(~1,
                     data=colData(sample_deseq_data))

sample_svseq <- svaseq(dat, mod, mod0)

svseq <- sample_svseq

# Generating the correlation matrix
bc_corr_plot(bc_method = "sva", var_int = "Group", font_size = 1)


ddssva <- sample_deseq_data

# Add all the surrgate variables to the DESeq Dataset
ddssva@colData <- cbind(ddssva@colData, sv_df)


design(ddssva) <-
  ~ SV1+ SV2+ SV3+ SV4+ SV5+
  SV6 + Group

# Reperform Deseq with the surrogate variables
de_sva_data <- DESeq(ddssva, parallel = T)


# PCA

vst_sv <- vst(de_sva_data, blind=FALSE)

pca_data_sva <- pca_prep(vst_sv, "Group")

rna_pca_plot(pca_data_sva)


# DEG Results
de_sva_res <- results(de_sva_data)

res_deg_sva <- clean_deseq_res(de_sva_res)


# Save the results as .Rdata to be used
usethis::use_data(sample_deseq_data, overwrite = TRUE)
usethis::use_data(sample_deseq_res, overwrite = TRUE)
usethis::use_data(sample_clean_deg, overwrite = TRUE)
usethis::use_data(sample_vst_dd, overwrite = TRUE)
