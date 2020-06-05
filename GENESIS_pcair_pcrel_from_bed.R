##### If you are running on Rosalind (our HPC server) it is best to use our
##### GENESIS package repository. Uncomment the line below (delete the #)
.libPaths( c( "/gpfs/share/TICR/Software/R/x86_64-pc-linux-gnu-library/3.3/" , .libPaths() ) )

##### If GENESIS is not installed on your system use these commands to install it (delete the #)

# install.packages("BiocManager")
# BiocManager::install("GENESIS")

##### If you have any errors installing GENESIS we will need to solve those first 


##### These are all of the libraries you will need to run GENESIS

library("GWASTools")
library("SNPRelate")
library("GENESIS")

##### End libraries section 

##### Update this section with the names of your data and output files #####

bed <- "CAAPA2_Ethiopia_MEGA_012220_updated_callrate_passing_QC.bed"
bim <- "CAAPA2_Ethiopia_MEGA_012220_updated_callrate_passing_QC.bim"
fam <- "CAAPA2_Ethiopia_MEGA_012220_updated_callrate_passing_QC.fam"
kin0 <- "CAAPA2_Ethiopia_MEGA_012220_updated_callrate_passing_QC.kin0"

# This is what you want your output files named.
prefix <- "CAAPA2 Ethiopian_Saliva Samples Received_10282019_JWM_test"

##### End section to update with data and output files #####

snpgdsBED2GDS(bed.fn = bed, bim.fn = bim, fam.fn = fam, out.gdsfn = paste(prefix, ".gds", sep=""))
geno <- GdsGenotypeReader(filename = paste(prefix, ".gds", sep=""))
genoData <- GenotypeData(geno)
iids <- getScanID(genoData)

Kingmat <- king2mat(file.kin0=kin0,file.kin=NULL,type="kinship",iids = iids)

pcmat_num <- as.integer(1)
mypcair <- pcair(genoData = genoData, kinMat = Kingmat,divMat = Kingmat, v=200)
mypcrel <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:pcmat_num], training.set = mypcair$unrels)

save(mypcair, mypcrel, file = paste(prefix, "_pcair_pcrel.Rdata", sep=""))