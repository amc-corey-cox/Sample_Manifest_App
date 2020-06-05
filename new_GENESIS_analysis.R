library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(SeqVarTools)
library(GENESIS)

library(magrittr)
library(tidyverse)

load("input_to_genesis/adrn_AA_GENESIS_PC_pheno_covar_data.Rdata")

pheno_file <- "adrn_AA_GENESIS_pheno.txt"
pheno <- read_tsv(pheno_file, col_types = "ccddcdddddd")
phenotypes <- c("DIAGROUP", "tIgE", "EASI", "EOS", "RL")


# testcovMat <- covMat[ScanAnnot$scanID, ScanAnnot$scanID]

ScanAnnot <- as_tibble(mypcair$vectors, rownames = "scanID",
                       .name_repair = ~ paste0("pc", 1:ncol(mypcair$vector))) %>%
  mutate(scanID = pheno$ID) %>%
  inner_join(pheno, by = c("scanID" = "ID")) %>%
  select("scanID", paste0("pc", 1:min(50, ncol(mypcair$vectors))),
         "GENDER", "AGE", all_of(phenotypes)) %>% as.data.frame() %>%
  ScanAnnotationDataFrame()

covMat <- pcrelateToMatrix(mypcrel)[ScanAnnot$scanID, ScanAnnot$scanID]

# covMatList <- list("Kin" = pcrelateToMatrix(mypcrel))

save(ScanAnnot, covMat, phenotypes, file = "myScanData.Rdata")
load("myScanData.Rdata")

get_pc_pval <- function(pc, pheno, scans, matList, varType) {
  fitNullModel(x = scans, outcome = pheno, covars = pc, cov.mat = matList, family = varType) %>%
    extract2("fixef") %>% pluck("pval", 2)
}

get_covars <- function(pheno, scans, matList, vartype) {
  paste0("pc", 1:4) %>% set_names() %>% 
    map(~ get_pc_pval(., pheno, scans, matList, varType)) %>%
    keep(~ .x < 0.05)
}

covars <- phenotypes %>% set_names() %>%
  map(~ get_covars(., ScanAnnot, covMat, varType)) %>% map(names)

# 
gdsfile <- tempfile()

bfile <- "adrn_AA_GENESIS"
snpgdsBED2GDS(bed.fn = paste0(bfile, ".bed"), bim.fn = paste0(bfile, ".bim"), fam.fn = paste0(bfile, ".fam"), out.gdsfn = gdsfile)

geno <- GdsGenotypeReader(gdsfile)
genoData <- GenotypeData(geno)
genoDataIt <- GenotypeBlockIterator(genoData)

null_model <- fitNullModel(x = ScanAnnot, outcome = "RL", covars = c(covars$RL, "GENDER", "AGE"), cov.mat = covMat, family = varType)
assoc_test <- assocTestSingle(gdsobj = genoDataIt, null.model = null_model)

write_tsv(assoc_test, paste0(bfile, "_GENESIS_results.txt"))