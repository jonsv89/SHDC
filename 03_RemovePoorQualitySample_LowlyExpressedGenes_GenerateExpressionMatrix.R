## CODE FOR THE DIFFERENTIAL EXPRESSION ANALYSIS OF EACH CASE SAMPLE AGAINST ALL THE CONTROL SAMPLES FROM THE SAME STUDY ##
## Developed by Jon Sanchez-Valle
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es

## load the essential packages 
library("data.table")
library("affy")
library(frma)
library("limma")
require(graphics)
library("hgu133plus2frmavecs")
library("massiR")
library("hgu133plus2.db")
library("cluster")
library("NbClust")
library(e1071)
library("GEOmetadb")
library(Biobase)
library(hgu133plus2CellScore)
## This package is not working ##
library("genefilter")
gpl570<-fread("Microarrays/Data/GPL570-55999.txt",stringsAsFactors = F,sep="\t",skip = 15,header=T)
probeandgenesymbol<-gpl570[,c(1,11)]
probeandgenesymbol2<-as.matrix(probeandgenesymbol)
rownames(probeandgenesymbol2)<-probeandgenesymbol$ID
## Identify probes pointing to several gene symbols ##
morethan1genesymbol<-probeandgenesymbol$ID[grep("///",probeandgenesymbol$`Gene Symbol`)]
## Identify probes with no gene symbol ##
nogenesymbol<-probeandgenesymbol$ID[which(probeandgenesymbol$`Gene Symbol`=="")]
## Probes to remove ##
probesremoves<-c(morethan1genesymbol,nogenesymbol)

## Create a directory to write the table indicating the samples to be removed ##
if("Targets"%in%list.files("Microarrays/Data/")==FALSE){dir.create("Microarrays/Data/Targets")}

# ## Load each study and identify low quality samples ##
# fileswithknownsex<-list.files("Microarrays/Data/Raw_data_predictedsex/")
# resumen<-c()
# for(a in 1:length(fileswithknownsex)){
#   # a<-1
#   sextab<-read.csv2(paste("Microarrays/Data/Raw_data_predictedsex/",fileswithknownsex[a],sep=""),stringsAsFactors = F,sep="\t")
#   ## Create the targets table ##
#   targets<-cbind(paste("Microarrays/Data/Raw_data/",gsub(".txt","",fileswithknownsex[a]),"/",sextab$Category,"/",sextab$Sample,sep=""),sextab,gsub(".+_","",gsub(".txt","",fileswithknownsex[a])))
#   colnames(targets)<-c("Path","Sample","Category","Sex","Study")
#   ## Read the affy files ##
#   affyBatch<-ReadAffy(filenames=targets$Path)
#   ## Normalize ##
#   expSet<-frma(affyBatch)
#   ## Calculate GNUSE median values for each sample, and remove those with a value higher than 1.25
#   ## What means that you remove those arrays whose precission is, by mean, a 25% worst that the usual array
#   gnuseres<-GNUSE(expSet,type="stats")
#   gnuseres<-gnuseres[,targets$Sample]
#   removeornot<-rep("no",length(targets$Sample))
#   if(length(which(gnuseres[1,]>1.25))>0){removeornot[which(gnuseres[1,]>1.25)]<-"yes"}
#   targets<-cbind(targets,removeornot)
#   colnames(targets)[6]<-"RemoveThem"
#   targets<-cbind(targets,gnuseres[1,])
#   colnames(targets)[7]<-"MedianGNUSE"
#   write.table(targets,paste("Microarrays/Data/Targets/",fileswithknownsex[a],sep=""),quote=F,sep="\t",row.names=F)
#   resumen<-rbind(resumen,targets)
#   print(a/length(fileswithknownsex))
# }
# pdf(file="Microarrays/Plots/Histogram_GNUSE.pdf")
#   hist(as.numeric(resumen$MedianGNUSE),breaks=100,xlim=c(0.5,2.5),xlab="median GNUSE",main="Histogram of median GNUSE\nmicroarray HG-U-133 plus 2")
# dev.off()
# 
# ## Remove the samples that should not be used ##
# ## @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ ##
# if("FinalTargets"%in%list.files("Microarrays/Data/")==FALSE){dir.create("Microarrays/Data/FinalTargets")}
# targetdatasets<-list.files("Microarrays/Data/Targets/")
# for(a in 1:length(targetdatasets)){
#   # a<-1
#   tt<-read.csv2(paste("Microarrays/Data/Targets/",targetdatasets[a],sep=""),stringsAsFactors = F,sep="\t")
#   if(length(which(tt$RemoveThem=="yes"))>0){tt<-tt[-which(tt$RemoveThem=="yes"),]}
#   if(dim(tt)[1]>0){
#     write.table(tt,paste("Microarrays/Data/FinalTargets/",targetdatasets[a],sep=""),quote=F,sep="\t",row.names=F)
#   }
#   if(dim(tt)[1]==0){print(paste("We lost the study:",targetdatasets[a]))}
# }

#### Identify the number of females and males by cases and controls: ####
## We need at least 3 samples per category and sex ##
fileswithknownsex<-list.files("Microarrays/Data/FinalTargets/")
resultstab<-c()
for(a in 1:length(fileswithknownsex)){
  # a<-1
  sextab<-read.csv2(paste("Microarrays/Data/FinalTargets/",fileswithknownsex[a],sep=""),stringsAsFactors = F,sep="\t")
  cases<-sextab$Sex[which(sextab$Category=="Case")]
  controles<-sextab$Sex[which(sextab$Category=="Control")]
  resultstab<-rbind(resultstab,c(length(which(cases=="female")),length(which(controles=="female")),length(which(cases=="male")),length(which(controles=="male"))))
}
rownames(resultstab)<-gsub(".txt","",fileswithknownsex)
colnames(resultstab)<-c("fem_case","fem_control","mal_case","mal_control")
# ## Datasets with at least 3 cases and 3 controls for one gender or both
# femdatasets3<-rownames(resultstab)[intersect(which(resultstab[,1]>=3),which(resultstab[,2]>=3))] # 137 datasets
# maldatasets3<-rownames(resultstab)[intersect(which(resultstab[,3]>=3),which(resultstab[,4]>=3))] # 157 datasets

# cogerparalgo<-sort(unique(c(intersect(which(resultstab[,1]>=3),which(resultstab[,2]>=3)),intersect(which(resultstab[,3]>=3),which(resultstab[,4]>=3)),
#                             intersect(which(resultstab[,1]>=3),which(resultstab[,3]>=3)),intersect(which(resultstab[,2]>=3),which(resultstab[,4]>=3)))))
# thestudieswekeep<-resultstab[cogerparalgo,]
# 
# casconwomen<-rep("",length(thestudieswekeep[,1])) ; casconwomen[intersect(which(thestudieswekeep[,1]>=3),which(thestudieswekeep[,2]>=3))]<-"case-control women"
# casconmen<-rep("",length(thestudieswekeep[,1])) ; casconmen[intersect(which(thestudieswekeep[,3]>=3),which(thestudieswekeep[,3]>=3))]<-"case-control men"
# 
# womenmencas<-rep("",length(thestudieswekeep[,1])) ; womenmencas[intersect(which(thestudieswekeep[,1]>=3),which(thestudieswekeep[,3]>=3))]<-"women-men case"
# womenmencon<-rep("",length(thestudieswekeep[,1])) ; womenmencon[intersect(which(thestudieswekeep[,2]>=3),which(thestudieswekeep[,4]>=3))]<-"women-men control"
# 
# thestudieswekeep<-cbind(thestudieswekeep,casconwomen,casconmen,womenmencas,womenmencon)
# colnames(thestudieswekeep)[5:8]<-c("cascon_women","cascon_men","womenmen_cas","womenmen_con")
# write.table(thestudieswekeep,"Microarrays/Data/Samples_per_study_including_potential_comparisons.txt",quote=F,sep="\t",row.names=F)

icdcode<-"ICD9"
comparison<-"Women"
removepoorqualityandlowexpression<-function(icdcode,comparison,resultstab){
  if(icdcode%in%list.files("Microarrays/")==FALSE){dir.create(paste("Microarrays/",icdcode,sep=""))}
  if(comparison%in%list.files(paste("Microarrays/",icdcode,sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/",comparison,sep=""))}
  if("Targets"%in%list.files(paste("Microarrays/",icdcode,"/",comparison,sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/",comparison,"/Targets",sep=""))}
  if("GeneExpressionMatrix"%in%list.files(paste("Microarrays/",icdcode,"/",comparison,sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/",comparison,"/GeneExpressionMatrix",sep=""))}
  if("DifferentialExpressions"%in%list.files(paste("Microarrays/",icdcode,"/",comparison,sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/",comparison,"/DifferentialExpressions",sep=""))}
  
  ## Read the table with the ICD9 information ##
  tabla<-read.csv2("Microarrays/Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
  
  if(icdcode=="ICD9"){
    if(length(which(tabla$ICD9=="ICD9_-"))>0){tabla<-tabla[-which(tabla$ICD9=="ICD9_-"),]}
    nametoicd<-tabla$ICD9 ; names(nametoicd)<-tabla$Disease
  }
  
  if(icdcode=="ICD10"){
    if(length(which(tabla$ICD10=="-"))>0){tabla<-tabla[-which(tabla$ICD10=="-"),]}
    nametoicd<-tabla$ICD10 ; names(nametoicd)<-tabla$Disease
  }
  
  if(comparison=="Women"){
    sexo<-"female"
    compadatasets3<-rownames(resultstab)[intersect(which(resultstab[,1]>=3),which(resultstab[,2]>=3))]
  }
  if(comparison=="Men"){
    compadatasets3<-rownames(resultstab)[intersect(which(resultstab[,3]>=3),which(resultstab[,4]>=3))]
    sexo<-"male"
  }
  if(comparison=="Cases"){
    compadatasets3<-rownames(resultstab)[intersect(which(resultstab[,1]>=3),which(resultstab[,3]>=3))]
    categoria<-"Case"
  }
  if(comparison=="Controls"){
    compadatasets3<-rownames(resultstab)[intersect(which(resultstab[,2]>=3),which(resultstab[,4]>=3))]
    categoria<-"Control"
  }
  if(comparison=="Adjusted"){
    compadatasets3<-unique(c(rownames(resultstab)[intersect(which(resultstab[,1]>=3),which(resultstab[,2]>=3))],
                              rownames(resultstab)[intersect(which(resultstab[,3]>=3),which(resultstab[,4]>=3))]))
  }
  
  
  ## Transform disease name datasets into ICD9 datasets
  compaicd<-as.character(nametoicd[gsub("_.+","",compadatasets3)])
  compadatasets<-paste(gsub("_","-",compaicd),gsub(".+_","",compadatasets3),sep="_")
  if(length(which(is.na(compaicd)))>0){
    compadatasets<-compadatasets[-which(is.na(compaicd))]
    compadatasets2<-compadatasets3[-which(is.na(compaicd))]
  }
  print(paste("Starting ",icdcode," - ",comparison,":",sep=""))
  compadiseases<-unique(gsub("_.+","",compadatasets))
  for(a in 1:length(compadiseases)){
    # a<-35
    datasets<-compadatasets2[grep(paste("^",compadiseases[a],"_",sep=""),compadatasets)]
    targets<-data.frame()
    for(b in 1:length(datasets)){
      # b<-1
      target<-read.csv2(paste("Microarrays/Data/FinalTargets/",datasets[b],".txt",sep=""),stringsAsFactors = F,sep="\t")
      if(comparison=="Women" || comparison=="Men"){
        targets<-rbind(targets,target[which(target$Sex==sexo),c(1:3,5,4)])
      }
      if(comparison=="Cases" || comparison=="Controls"){
        targets<-rbind(targets,target[which(target$Category==categoria),c(1:3,5,4)])
      }
      if(comparison=="Adjusted"){
        ## There are enough females
        femind<-c()
        if(length(intersect(which(target$Sex=="female"),which(target$Category=="Case")))>=3 && length(intersect(which(target$Sex=="female"),which(target$Category=="Control")))>=3){
          femind<-which(target$Sex=="female")
        }
        ## There are enough males
        malind<-c()
        if(length(intersect(which(target$Sex=="male"),which(target$Category=="Case")))>=3 && length(intersect(which(target$Sex=="male"),which(target$Category=="Control")))>=3){
          malind<-which(target$Sex=="male")
        }
        indice<-sort(c(femind,malind))
        targets<-rbind(targets,target[indice,c(1:3,5,4)])
      }
    }
    colnames(targets)<-c("Path","Sample","Category","Study","Sex")
    write.table(targets,paste("Microarrays/",icdcode,"/",comparison,"/Targets/",compadiseases[a],".txt",sep=""),quote=F,sep="\t",row.names=F)
    ## If there are duplicated samples remove them ##
    if(length(which(duplicated(targets$Sample)))>0){
      print(paste(a,"tiene rows duplicadas:",compadiseases[a]))
      targets<-targets[-which(duplicated(targets$Sample)),]
      ## Write targets file again ##
      write.table(targets,paste("Microarrays/",icdcode,"/",comparison,"/Targets/",compadiseases[a],".txt",sep=""),quote=F,sep="\t",row.names=F)
    }
    
    ## Evaluate if controlling for the experiment while having just one experiment gives different results
    affyBatch<-ReadAffy(filenames=targets$Path)
    ## Identify lowly expressed probes ##
    co <- mas5calls(affyBatch)
    co <- assayData(co)[["se.exprs"]] #extract detection p-values pvalue.detection.cutoff <- 0.05
    pvalue.detection.cutoff <- 0.05
    calls.matrix <- co < pvalue.detection.cutoff
    extranumbers<-function(vector){valor<-length(which(vector));return(valor)}
    totalexpressed<-apply(calls.matrix, 1, extranumbers)
    ## Identify the lowly expressed probes ##
    seleccionar<-which(totalexpressed>min(table(targets[,3:4])))
    expSet<-frma(affyBatch)
    comoda2<-exprs(expSet)
    ## Remove the lowly expressed probes ##
    comoda3<-comoda2[seleccionar,]
    ## Remove the probes not mapping a gene or mapping to more than one gene ##
    comoda<-comoda3[setdiff(rownames(comoda3),probesremoves),]
    coluna<-colnames(comoda)
    laseleccion<-probeandgenesymbol2[rownames(comoda),]
    genes<-unique(laseleccion[,2])
    expression<-c()
    ## convert the probes' expression matrix into a gene expression matrix, median values are calculated when several probes refer to the same gene symbol
    for (c in 1:length(genes)){
      com<-as.numeric(which(laseleccion[,2]==as.character(genes[c])))
      if (length(com)>1){
        iq<-apply(comoda[com,],2,median)
      }
      if (length(com)==1){
        iq<-comoda[com,]
      }
      expression<-rbind(expression,iq)
    }
    colnames(expression)<-coluna
    rownames(expression)<-genes
    write.table(expression,paste("Microarrays/",icdcode,"/",comparison,"/GeneExpressionMatrix/",compadiseases[a],".txt",sep=""),quote=F,sep="\t")
    print(paste(round((a/length(compadiseases))*100,3),"%",sep=""))
  }
  print("Finished")
}

## ICD9 ##
removepoorqualityandlowexpression("ICD9","Women",resultstab)
removepoorqualityandlowexpression("ICD9","Men",resultstab)
removepoorqualityandlowexpression("ICD9","Cases",resultstab)
removepoorqualityandlowexpression("ICD9","Controls",resultstab)
removepoorqualityandlowexpression("ICD9","Adjusted",resultstab)

## ICD10 ##
removepoorqualityandlowexpression("ICD10","Women",resultstab)
removepoorqualityandlowexpression("ICD10","Men",resultstab)
removepoorqualityandlowexpression("ICD10","Cases",resultstab)
removepoorqualityandlowexpression("ICD10","Controls",resultstab)
removepoorqualityandlowexpression("ICD10","Adjusted",resultstab)








