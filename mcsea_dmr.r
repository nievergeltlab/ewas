
#To install mCSEA
# source("https://bioconductor.org/biocLite.R")
# biocLite("mCSEA")



#Run as an interactive job or else just make this into an Rscript and call it from a job (or run it from your desktop)
#qsub -I -V -A sheering_fluxod -l nodes=1:ppn=2,qos=flux,walltime=2:00:00,pmem=2000mb -q fluxod

#call into working directory
cd /scratch/sheering_fluxod/adammai/methylation/

#Load R
module load R/3.5.0
R

#Load mCSEA
 library(mCSEA)

#Load EWAS summary statistics. '
#I am assuming results are formatted like they are in the standard scripts circulated by andrew
 load('pris_ruttenhc3_results_aug27_noutcontrolsnocov.Rdata')

# Load data
 load("pris_v1b.methylation") # beta matrix. CpG x Participants
 dat2 <- as.data.frame(combat_fixed_logit_back_noimputations)

#Read in covariates. 

#What I am doing here is loading PTSD phenotypes, cell types, then merging them. 
#The data I am working with is in long format. I convert it to wide format.

 batch0 <- read.table("pheno2.txt", header=T,na.strings=c("NA","#N/A","N/A","999"),stringsAsFactors=F)
 cellcounts <- read.csv('pris_celltypes.csv')
 names(cellcounts) [1] <- "methylome_id"
 batcha <- merge(batch0,cellcounts,by="methylome_id")
 keepidsa <- table(batcha$id)
 keepids <- names(keepidsa[keepidsa == 2])

 batchx <- subset(batcha, !is.na(ZIL_T39_6mnd) & Group != "3" & methylome_id %in% names(dat2) & id %in% keepids) #Remove untraumatized controls!

 batch_t0 <- subset(batchx,time == 1)
 batch_t1 <- subset(batchx,time == 2)

 batch1 <- merge(batchx,batch_t0,by="id",suffixes=c("","_v0"))
 batch <- merge(batch1,batch_t1,by="id",suffixes=c("","_v1"))


#Subset data to only subjects that are in the methylation data
 covariates_ordered <- subset(batch, methylome_id %in% names(dat2))

#Only take select data with pre and post observations (if data is not already filtered)
 #keepids <- read.table('cases_controls_prepost',header=F,stringsAsFactors=F)[,1]
 #covariates_ordered_subset0 <- subset(covariates_ordered, methylome_id %in% keepids[,1])
 
 #my data is filtered so I just proceed on without modifying it
 covariates_ordered_subset <- subset(covariates_ordered)


#order by id, because GEE functions usually need this to be ordered, otherwise factors not correct
 covariates_ordered_subsetX <- covariates_ordered_subset[order(covariates_ordered_subset$id),]


#Reorder subjects in methylation matrix to the order they are in the covariate sheet. Only take mutually overlapping subjects
 beta.norm <- as.data.frame(as.matrix(dat2[, covariates_ordered_subsetX$methylome_id]))

#rm(dat) #remove extra data to save some memory

#This is a little redundant here, I'm just renaming the covariate sheet
 phenox <- subset(covariates_ordered_subsetX) 

#Set the row names of the phenotype sheet to the subject IDs
 row.names(phenox) <- phenox$methylome_id

 
 all(rownames(phenox)==colnames(beta.norm)) # this should be TRUE - checks if phenotype and genotype matrices are lined up correctly

#Make any last alterations in phenotype

 phenox$agedif <- phenox$age_v1 - phenox$age_v0
 phenox$Bcelld <- phenox$Bcell_v1 - phenox$Bcell_v0
 phenox$Monod <- phenox$Mono_v1 - phenox$Mono_v0
 phenox$CD4Td <- phenox$CD4T_v1 - phenox$CD4T_v0
 phenox$CD8Td <- phenox$CD8T_v1 - phenox$CD8T_v0
 phenox$NKd <- phenox$NK_v1 - phenox$NK_v0

 phenoZ <- subset(phenox,select=c("agedif", "CD8Td","CD4Td","NKd","Bcelld","Monod","ZIL_T39_6mnd_v1"))

 bvec <- resultsT[-1,3]
  
#Set the seed then do the enrichment test
#Note that the test can be run on genes, cpg islands, or promoter regions. I run all 3.

 set.seed(17)
 prismoGenes <- mCSEATest(rank=bvec, methData=beta.norm, pheno=phenoZ,column="ZIL_T39_6mnd_v1", regionsTypes = "genes", platform = "450k",nproc=4,nperm=100000)
 save(prismoGenes,file='prismo_genes_gsea_100k.r')
  
 set.seed(17)
 prismoCGI <- mCSEATest(rank=bvec, methData=beta.norm, pheno=phenoZ,column="ZIL_T39_6mnd_v1", regionsTypes = "CGI", platform = "450k",nproc=4,nperm=100000)
 save(prismoCGI,file='prismo_cgi_gsea_100k.r')
 
 set.seed(17)
 prismoPromoters <- mCSEATest(rank=bvec, methData=beta.norm, pheno=phenoZ,column="ZIL_T39_6mnd_v1", regionsTypes = "promoters", platform = "450k",nproc=4,nperm=100000)
 save(prismoPromoters,file='prismo_promoters_gsea_100k.r')

#Show results
head(prismoGenes$genes)
head(prismoCGI$CGI)
head(prismoGenes$genes)


