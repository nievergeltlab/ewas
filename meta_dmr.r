#Note that this is for just CpG islands. The other options are genes and promoters
#So I suggest using find and replace on "CGI" to get the appropriate file for genes or promoters, making a copy for each version of this.

#Load DMR analysis results.
load('semper/mrs_cgi_gsea_100k.r')
load('mirror/prismo_cgi_gsea_100k.r')
load('starrs_cgi_gsea_100k.r')

#List the number of subjects in the analysis
 mrs = 126
 starrs=78
 prismo=62

#For each dataset:
#Make a dataframe with the p-value and normalized effect size (direction)
#Name the dataframe rows by the CpG island names
#Convert p-values to Z scores

#dataset 1
 g1 <- data.frame(mrsCGI$CGI$pval,mrsCGI$CGI$NES)
 names(g1) <- c('p_mrs','nes_mrs')
 row.names(g1) <- row.names(mrsCGI$CGI)
 g1$chi_mrs <- qchisq(g1$p_mrs,1,lower.tail=F)
 g1$z_mrs <- sqrt(g1$chi_mrs)* sign(g1$nes_mrs) 
 g1$CpG <- row.names(g1)

#dataset 2
 g2 <- data.frame(starrsCGI$CGI$pval,starrsCGI$CGI$NES)
 names(g2) <- c('p_starrs','nes_starrs')
 row.names(g2) <- row.names(starrsCGI$CGI)
 g2$chi_starrs <- qchisq(g2$p_starrs,1,lower.tail=F)
 g2$z_starrs <- sqrt(g2$chi_starrs)* sign(g2$nes_starrs) 
 g2$CpG <- row.names(g2)

#dataset 3
 g3 <- data.frame(prismoCGI$CGI$pval,prismoCGI$CGI$NES)
 names(g3) <- c('p_prismo','nes_prismo')
 row.names(g3) <- row.names(prismoCGI$CGI)
 g3$chi_pris <- qchisq(g3$p_prismo,1,lower.tail=F)
 g3$z_pris <- sqrt(g3$chi_pris)* sign(g3$nes_prismo) 
 g3$CpG <- row.names(g3)

#Merge the 3 sets of results together according to DMR name
Bmat0 <- merge(g1,g2,by="CpG",suffixes=c("_mrs","starrs"))
Bmat3 <- merge(Bmat0,g3,by="CpG")


#Weight Z by sample size
 Bmat3$bw1 <- Bmat3$z_starrs *sqrt(78) #31 Cases and 47 controls = 78
 Bmat3$bw2 <- Bmat3$z_mrs *sqrt(126) #63 and 63 = 126 
 Bmat3$bw3 <- Bmat3$z_pris *sqrt(62) #29 and 33 = 62


#Sum of weighted Zs (numerator for sz pool)
 Bmat3$bwsum <- Bmat3$bw1 + Bmat3$bw2   + Bmat3$bw3

#Sum of weights (denominator for sz pool)
 Bmat3$wsum <- sqrt( 78 + 126 + 62) 

#Pooled sz (sum of weighted sz, averaged by weight)
 Bmat3$szpool <- Bmat3$bwsum / Bmat3$wsum

#P-value
 Bmat3$p_meta <- pchisq(Bmat3$szpool^2,1,lower.tail=F)

#Get top p-values for each analysis (per study and meta analysis
 r1 <- row.names(Bmat3[Bmat3$p_starrs < 0.00017,])
 r1 <- c(r1,row.names(Bmat3[Bmat3$p_mrs < 0.00017,]))
 r1 <- c(r1,row.names(Bmat3[Bmat3$p_pris < 0.00017,]))
 r1 <- c(r1,row.names(Bmat3[Bmat3$p < 0.00017,]))

#Filter meta analysis data down to just these markers
 sigmarkers <- Bmat3[row.names(Bmat3) %in% r1,]

#Save the significant DMRs to a file
 write.table(sigmarkers,'gsea_CGI_dmrs_meta_100k.txt')

#Also save the entire meta analysis to a file
 save(Bmat3,file="gsea_CGI_dmr_meta_100k_ruttenhc3.R")


####
#Get annotation info for DMRs, i.e. the span of the DMR .. WORK IN PROGRESS

#
islands <- c("chr1:156814881-156815792","chr6:33048416-33048814","chr21:35831697-35832365","chr7:1885033-1885402","chr7:27169572-27170638","chr6:32551851-32552331","chr6:25882327-25882560","chr7:1885033-1885402",	"chr11:70672834-70673055"	,"chr12:9217328-9217715",	"chr7:27169572-27170638",	"chr6:33048416-33048814",	"chr10:530713-531099",	"chr17:76037074-76037323"	,"chr21:35831697-35832365",	"chr17:8702342-8702824",	"chr16:1583809-1584641"	,"chr5:191792-192544")
head(starrsCGI$CGI[row.names(starrsCGI$CGI) %in% islands,])

c("HLA-DPB1","HLA-DRB1","IFT140","PRDM16","HLA-DPB1","KIF25" )
cpgs_assoc <- starrsCGI$CGI_association[islands]

save(cpgs_assoc , file="gsea_CGI_dmr_meta_100k_ruttenhc3_CpGs.R")
load("gsea_CGI_dmr_meta_100k_ruttenhc3_CpGs.R")



 library(data.table)
 annots <- fread('misc/humanmethylation450_15017482_v1-2.csv',data.table=F,skip=7,nrow=485553)
 row.names(annots) <- annots$IlmnID 

 cmat <- as.data.frame(matrix(nrow=length(cpgs_assoc),ncol=2))

 cmat$length <- NA
 cmat$nprobes <- NA
 cmat$start <-NA 
 cmat$stop <- NA
 cmat$Gene <- NA

 for (num in 1:length(names(cpgs_assoc)))
 {
   ma <- subset(annots, IlmnID %in% unlist(cpgs_assoc[num][1]))
   cmat[num,]$Gene <- ma$UCSC_RefGene_Name[1]
   cmat[num,]$start <- min(ma$MAPINFO)
   cmat[num,]$stop <- max(ma$MAPINFO)
   cmat[num,]$nprobes=dim(ma)[1]
   
 }

 cmat


#for islands, just need nprobes
 dmat <- as.data.frame(names(cpgs_assoc))
 dmat$length <- NA
 for ( i in 1:length(cpgs_assoc))
 {
  dmat[i,]$length <- length(unlist(cpgs_assoc[i]))
 }


#For a given study, plot the dmr stuff (WORK IN PROGRESS, mCSEA plot needs to be hacked for this to work!!)

#Load EWAS statistics
load('semper/MRS_ruttenhc3_results.Rdata')

#Plot DMR
pdf('test_HLACGI.pdf',7,7)
mCSEAPlot(mrsCGI, regionType = "CGI",
dmrName = "chr6:33048416-33048814",
transcriptAnnotation = "symbol", makePDF = TRUE)
dev.off()

#PLot enrichment
pdf('test_HLACGI.pdf',7,7)
mCSEAPlotGSEA(resultsT[-1,3], mrsCGI, regionType = "CGI", dmrName = "chr6:33048416-33048814")
dev.off()


