#These DMR plots are fucking trash. I need to load my results, load the changes per study, take the annotations and plot that..

#I already pulled the start-stop. I need to load those and the individual study results


#Make a file with the relevant DMRs,
#Should look like this 
# Gene	chr	start_pos	stop_pos	n_cpgs_in_dmr
# HLA-DPB1	6	33043976	33054001	56
# HLA-DRB1	6	32547019	32557404	32
# IFT140	16	1561036	1652552	121
# KIF25	6	168433191	168445614	18
# PRDM16	1	2986362	3349982	608
# KCNE1	21	35827824	35884508	23
# TRMT12	8	125461772	125464547	9
# HOXA4	7	27169740	27171528	24
# MFSD6L	17	8700574	8703341	12
# LOC144571	12	9217079	9217769	9
# SLC17A3	6	25882327	25882560	8
# NTRK1	1	156814881	156815792	9
# TNRC63	17	76037074	76037323	9
# LINC00612	12	9217328	9217715	13
# MFSD6L	17	8702342	8702824	14
# SHANK2	11	70672834	70673055	14
# KCNE1	21	35831697	35832365	18
# MAD1L1	7	1885033	1885402	18
# IFI140	16	1583809	1584641	19
# LRRC14B	5	191792	192544	19
# HLA-DRB1	6	32551851	32552331	28
# DIP2C	10	530713	531099	28
# HOXA4	7	27169572	27170638	41
# HLA-DBP1	6	33048416	33048814	46


library(data.table)
#Load position of DMRS
dmrs <- fread('top_dmrs_for_plotting.txt',data.table=F)

#Load the average changes in methylation for cases and controls
#these files should have the following columns
#mean cases 
#sd cases
#mean controls,
#sd controls
#IlmnID of CpG site

#With a row for each CpG site


load('STARRS_changes.Rdata')
resultsChanges <- as.data.frame(resultsChanges)
resultsChanges$IlmnID <- row.names(resultsChanges)
starrsChanges <- resultsChanges
load('semper/MRS_changes.Rdata')
resultsChanges <- as.data.frame(resultsChanges)
resultsChanges$IlmnID <- row.names(resultsChanges)
mrsChanges <- resultsChanges
load('mirror/pris_changes.Rdata')
resultsChanges <- as.data.frame(resultsChanges)
resultsChanges$IlmnID <- row.names(resultsChanges)
prisChanges <- resultsChanges




####
#Get annotation info (Illumina manifest) to get the span of the DM

annots <- fread('misc/humanmethylation450_15017482_v1-2.csv',data.table=F,skip=7,nrow=485553)
row.names(annots) <- annots$IlmnID 

#For each DMR, plot all 3 study results in a 3 paneled figure
library(plotrix)

for (dmr in 1:dim(dmrs)[1])
{
 dmrname=dmrs[dmr,]$Gene
 startpos=dmrs[dmr,]$start_pos
 stoppos=dmrs[dmr,]$stop_pos
 chrom=dmrs[dmr,]$chr
 annotsub <- subset(annots, CHR== chrom & MAPINFO >= startpos & MAPINFO <= stoppos)
 annotids <- annotsub$IlmnID
 
 
 
 
 pdf(paste('dmrplots2/',dmrname,'_',startpos,'.pdf',sep=''))
 par(mfrow=c(3,1), oma = c(5,2,0,2) + 0.1, mar = c(1,.5,5,.5) + 0.1)
 
 #STARRS
 changesubstarrs <- merge(starrsChanges,annotsub,by="IlmnID",suffixes=c("","ano"))
 changesubstarrs <- changesubstarrs[order(changesubstarrs$MAPINFO),]
 
 meancases=changesubstarrs$mean_cases
 meancontrols=changesubstarrs$mean_controls
 positions=changesubstarrs$MAPINFO
  
 lcica=c(changesubstarrs$mean_cases - 1.96*(changesubstarrs$sd_cases/sqrt(78)))
 ucica=c(changesubstarrs$mean_cases + 1.96*(changesubstarrs$sd_cases/sqrt(78)))
 lcico=c(changesubstarrs$mean_controls - 1.96*(changesubstarrs$sd_controls/sqrt(78)))
 ucico=c(changesubstarrs$mean_controls + 1.96*(changesubstarrs$sd_controls/sqrt(78)))
 minplot <- min(lcica,lcico,ucica,ucico)       
 maxplot <- max(lcica,lcico,ucica,ucico)  
      
 plotCI(positions,meancases,li=lcica,ui=ucica,ylim=c(minplot,maxplot),pch=19,col=rgb(1,.3,.3,1),scol=rgb(1,0,0,.5),lwd=3,sfrac=0,cex.axis=1.45,cex.main=2)
 plotCI(positions,meancontrols,li=lcico,ui=ucico,ylim=c(minplot,maxplot),add=TRUE,pch=19,col=rgb(0,.3,1,1),scol=rgb(0,.3,1,.5),lwd=3,sfrac=0,cex.axis=1.45) 
 mtext("STARRS",1,line=3,cex=1.3)
 
 #MRS
 changesubmrs <- merge(mrsChanges,annotsub,by="IlmnID",suffixes=c("","ano"))
 changesubmrs <- changesubmrs[order(changesubmrs$MAPINFO),]
 
 meancases=changesubmrs$mean_cases
 meancontrols=changesubmrs$mean_controls
 positions=changesubmrs$MAPINFO
  
 lcica=c(changesubmrs$mean_cases - 1.96*(changesubmrs$sd_cases/sqrt(78)))
 ucica=c(changesubmrs$mean_cases + 1.96*(changesubmrs$sd_cases/sqrt(78)))
 lcico=c(changesubmrs$mean_controls - 1.96*(changesubmrs$sd_controls/sqrt(78)))
 ucico=c(changesubmrs$mean_controls + 1.96*(changesubmrs$sd_controls/sqrt(78)))
 minplot <- min(lcica,lcico,ucica,ucico)       
 maxplot <- max(lcica,lcico,ucica,ucico)  
      
 plotCI(positions,meancases,li=lcica,ui=ucica,ylim=c(minplot,maxplot),pch=19,col=rgb(1,.3,.3,1),scol=rgb(1,0,0,.5),lwd=3,sfrac=0,cex.axis=1.45,cex.main=2)
 plotCI(positions,meancontrols,li=lcico,ui=ucico,ylim=c(minplot,maxplot),add=TRUE,pch=19,col=rgb(0,.3,1,1),scol=rgb(0,.3,1,.5),lwd=3,sfrac=0,cex.axis=1.45) 
 mtext("MRS",1,line=3,cex=1.3)
 
 #PRISMO
 
 changesubprismo <- merge(prisChanges,annotsub,by="IlmnID",suffixes=c("","ano"))
 changesubprismo <- changesubprismo[order(changesubprismo$MAPINFO),]


 meancases=changesubprismo$mean_cases
 meancontrols=changesubprismo$mean_controls
 positions=changesubprismo$MAPINFO
  
 lcica=c(changesubprismo$mean_cases - 1.96*(changesubprismo$sd_cases/sqrt(78)))
 ucica=c(changesubprismo$mean_cases + 1.96*(changesubprismo$sd_cases/sqrt(78)))
 lcico=c(changesubprismo$mean_controls - 1.96*(changesubprismo$sd_controls/sqrt(78)))
 ucico=c(changesubprismo$mean_controls + 1.96*(changesubprismo$sd_controls/sqrt(78)))
 minplot <- min(lcica,lcico,ucica,ucico)       
 maxplot <- max(lcica,lcico,ucica,ucico)  
      
 plotCI(positions,meancases,li=lcica,ui=ucica,ylim=c(minplot,maxplot),pch=19,col=rgb(1,.3,.3,1),scol=rgb(1,0,0,.5),lwd=3,sfrac=0,cex.axis=1.45,cex.main=2)
 plotCI(positions,meancontrols,li=lcico,ui=ucico,ylim=c(minplot,maxplot),add=TRUE,pch=19,col=rgb(0,.3,1,1),scol=rgb(0,.3,1,.5),lwd=3,sfrac=0,cex.axis=1.45) 
 mtext("PRISMO",1,line=3,cex=1.3)
 
 dev.off()
 

 }
 
 

#I may have to plot changes in methylation, not just betas...
 #MAD1L1 plot only shows 2 CpGs. 
 #There is a logical error in that i am plotting based on the intersection of 3 datasets, but the DMR analysis may have slightly different CPGs
 #There should be 5+ CpGs in the mad1l1 analysis but there arent. Unless that is by their criteria and not mine .
 
 #p-values arent really giving me much apart from some kind of gene set test. I like this one https://openi.nlm.nih.gov/detailedresult.php?img=PMC4253444_oncotarget-05-9425-g006&req=4 
 #which means I need to calculate the average changes in controls and in cases
 #Is this going to be covariate adjusted??
 