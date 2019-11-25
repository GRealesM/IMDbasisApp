## Sample data set processing ##

## Date of creation: 2019-11-25

library(cupcake)
library(annotSnpStats)
library(data.table)
library(magrittr)

#### Load data sets contained at Cupcake ####
data(burren_imd_13_traits)

SPARSE_BASIS_EXTDATA <- system.file("extdata/sparse_imd_basis", package="cupcake")
SNP_MANIFEST_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/sparse_snp_manifest.RDS')
TRAIT_MANIFEST_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/sparse_trait_manifest.tab')
SHRINKAGE_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/13_trait_sparse_shrinkage.RDS')
BASIS_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/13_trait_sparse_basis.RDS')
GWAS_DATA_DIR <- file.path(SPARSE_BASIS_EXTDATA,'/gwas_data/')
basis.gwas.DT<-get_gwas_data(TRAIT_MANIFEST_FILE,SNP_MANIFEST_FILE,GWAS_DATA_DIR)
basis <- readRDS(BASIS_FILE)
shrink.DT <- readRDS(SHRINKAGE_FILE)

# Initially, we'll use one of the summary statistics file from Astle et al., 2016 (27863252) data, corresponding to a GWAS of Lymphocyte levels in Europeans.
# This data set is publicly available at GWAS catalog (ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/AstleWJ_27863252_GCST004627/harmonised/27863252-GCST004627-EFO_0004587.h.tsv.gz)

# We downloaded and renamed this file, so we can read it in R (It is a heavy file, so might not be open in a regular computer)
input <- fread(input = "LYC_Astle_27863252.tsv.gz")

# This data set was harmonized by the GWAS catalog team, but since Basis is in hg19/GRCh37, we'll use the original part of the file, and only the relevant columns
input <- input[, c(13:18, 21:23)]

# The bp is in the first column, so some steps are required
input$variant = as.character(lapply(strsplit(as.character(input$variant), split="_"), "[", 1))
input$variant = as.character(lapply(strsplit(as.character(input$variant), split=":"), "[", 2))

# Let's rename the columns to match the basis
setnames(input, c("chromosome", "variant", "other_allele", "effect_allele", "p_value", "beta", "standard_error"),
                c("CHR", "POS", "REF", "ALT", "P", "BETA", "SE"))


# Add some cols to both files
input[,alleles:=paste(REF,ALT,sep="/")][,pid:=paste(CHR,POS,sep=":")]
SNP.manifest[,alleles:=paste(ref_a1,ref_a2,sep="/")]

# Merged dataset
M <- merge(input,SNP.manifest[,.(pid,alleles)], by='pid', suffixes=c("",".manifest"))

# Check if everything is alright
all(g.class(M$alleles.manifest, M$alleles)== "nochange")

# Since everything is OK, let's make it a beautiful sample data set by removing some columns, renaming and reordering
M[, c("pid", "base_pair_location", "alleles", "alleles.manifest") :=NULL]
setnames(M, "variant_id", "SNP_ID")
setcolorder(M, c("SNP_ID","CHR","POS","REF" ,"ALT","BETA","SE","P"))

write.table(M, file = "../00-Raw/Sample_dataset_LYC_Astle_27863252.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

