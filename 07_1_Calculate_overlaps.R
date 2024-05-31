## CODE FOR CALCULATING THE OVERLAPS BETWEEN TRANSCRIPTOMIC NETWORKS AND EPIDEMIOLOGY ##
## Developed by Jon Sanchez-Valle
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es


#### Compare the obtained networks with epidemiology by gender ####
## load the essential packages 
library("igraph")
library("corrplot")
library("rstatix")
library("scales")
library("dendextend")
library("ggplot2")
library("gplots")
library("forestplot")
library("dplyr")
require("ggpubr")

args = commandArgs(trailingOnly=TRUE)

## Order the disease codes and save them ordered with the columns of interest for comparison purposes ##
# westergaard2<-read.csv2("Epidemiology/41467_2019_8475_MOESM6_ESM.txt",stringsAsFactors = F,sep="\t")
## We made a test and we verified that the interactions described by Westergaard et al. 2019, 
## Supplementary Data 3 (the one previously loaded) are the ones mentioned in table 1
## RR2.5 and RR97.5 are the lower and upper confidence intervals (R33  N30  2.18 [2.12-2.24]  1.82 [1.74-1.92])
## For our analysis, we are going to select all those interactions with a RR2.5 > 1.01
## Adjusted ##
# length(which(as.numeric(westergaard2$adjustedRR2.5)>1.01))

# westergaard<-c() ; for(a in 1:length(westergaard2[,1])){
  # westergaard<-rbind(westergaard,c(sort(as.character(westergaard2[a,1:2]),decreasing = F),as.character(westergaard2[a,c(7,13,19,25)])))
# }
# colnames(westergaard)<-colnames(westergaard2)[c(1,2,7,13,19,25)]
# write.table(westergaard,"Epidemiology/Epidemiology_ICD10.txt",quote=F,sep="\t",row.names=F)

# hidalgo<-read.csv2("Epidemiology/PDN_3_digits.net",stringsAsFactors = F,sep="\t",header=F)
# hidalgodiseases<-unique(c(hidalgo[,1],hidalgo[,2]))
# write.table(hidalgodiseases,"Epidemiology/Hidalgo_diseases_ICD9.txt",quote=F,sep="\t",row.names=F,col.names=F)
# hidalgo<-hidalgo[which(as.numeric(as.character(hidalgo[,7]))>1),c(1,2,6)]
# colnames(hidalgo)<-c("Disease1","Disease2","RR")
# write.table(hidalgo,"Epidemiology/Epidemiology_ICD9.txt",quote=F,sep="\t",row.names=F)

#### Calculate overlaps with epidemiology ####
# parte<-"calculateoverlaps"
parte<-args[1]
## @@ @@ @@ ##
## Function ##
## @@ @@ @@ ##
# experiment<-"Microarrays"
# icdcode<-"ICD10"
# comparison<-"Adjusted"
# thenet2<-read.csv2("Microarrays/Generated_Networks/ICD10_Adjusted_binarized_network.txt",stringsAsFactors = F,sep="\t")
# ## Sort the diseases ##
# if(icdcode=="ICD10"){thenet<-c() ; for(a in 1:length(thenet2[,1])){thenet<-rbind(thenet,c(sort(as.character(thenet2[a,1:2])),as.numeric(thenet2[a,3:dim(thenet2)[2]])))}}
# if(icdcode=="ICD9"){thenet<-c() ; for(a in 1:length(thenet2[,1])){thenet<-rbind(thenet,c(sort(as.numeric(thenet2[a,1:2])),as.numeric(thenet2[a,3:dim(thenet2)[2]])))}}
# colnames(thenet)<-colnames(thenet2)
# thenet[,1]<-as.character(thenet[,1]) ; thenet[,2]<-as.character(thenet[,2])
# numberrandom<-10000
## The function ##
calculateoverlaps<-function(experiment,icdcode,comparison,thenet,numberrandom){
  icds2<-unique(c(thenet[,1],thenet[,2]))
  nodesthatshouldbethere_epidemiology<-c()
  if(icdcode=="ICD9"){
    hidalgodiseases<-read.csv2("Epidemiology/Hidalgo_diseases_ICD9.txt",stringsAsFactors = F,header=F)[,1]
    ## Diseases in both networks ##
    icds<-intersect(icds2,hidalgodiseases)
    epidemiology2<-read.csv2("Epidemiology/Epidemiology_ICD9.txt",stringsAsFactors = F,sep="\t")
    epidemiology2[,1]<-as.character(epidemiology2[,1]) ; epidemiology2[,2]<-as.character(epidemiology2[,2])
    ## Select the comorbidities between shared diseases ##
    ind1<-c() ; ind2<-c() ; for(a in 1:length(icds)){ind1<-c(ind1,which(epidemiology2[,1]==icds[a])) ; ind2<-c(ind2,which(epidemiology2[,2]==icds[a]))}
    epidemiology<-epidemiology2[intersect(ind1,ind2),]
    west<-epidemiology
    ## Which are the diseases that should be included as isolated nodes in the network? ##
    nodesthatshouldbethere_epidemiology<-setdiff(icds,unique(c(west[,1],west[,2])))
  }
  
  if(icdcode=="ICD10"){
    epidemiology2<-read.csv2("Epidemiology/Epidemiology_ICD10.txt",stringsAsFactors = F,sep="\t")
    ## Diseases in both networks ##
    icds<-intersect(icds2,unique(c(epidemiology2$A,epidemiology2$B)))
    ## Select the comorbidities between shared diseases ##
    ind1<-c() ; ind2<-c() ; indl<-list()
    for(a in 1:length(icds)){
      ind1<-c(ind1,which(epidemiology2[,1]==icds[a]))
      ind2<-c(ind2,which(epidemiology2[,2]==icds[a]))
      indl[[icds[a]]]<-unique(c(which(epidemiology2[,1]==icds[a]),which(epidemiology2[,2]==icds[a])))
    }
    ## Is any of the common diseases sex-specific and so they might not have comorbidities for some comparisons?
    isallnas<-c()
    for(a in 1:length(names(indl))){
      # a<-17
      if(length(which(is.na(epidemiology2[indl[[a]],grep(paste("^",comparison,sep=""),colnames(epidemiology2),ignore.case = T)])==FALSE))==0){
        isallnas<-c(isallnas,names(indl)[a])
      }
    }
    ## Select comorbidities between shared diseases ##
    epidemiology<-epidemiology2[intersect(ind1,ind2),]
    ## Select the co-occurrences that are significant for the gender of interest in the epidemiology: adjusted, women or men ##
    west<-epidemiology[which(as.numeric(as.character(epidemiology[,grep(paste("^",comparison,sep=""),colnames(epidemiology),ignore.case = T)]))>1.01),] ### Changed from 1 to 1.01
    if(length(which(duplicated(west[,1:2])))>0){west<-west[-which(duplicated(west[,1:2])),]}
    ## Which are the diseases that should be included as isolated nodes in the network? ##
    nodesthatshouldbethere_epidemiology<-setdiff(icds,unique(c(west[,1],west[,2],isallnas)))
    ## If there are nodes that should not be taken into consideration as they are not in the epidemiology remove them ##
    if(length(isallnas)>0){
      ind1<-c() ; ind2<-c() ; for(a in 1:length(isallnas)){ind1<-c(ind1,which(thenet[,1]==isallnas[a])) ; ind2<-c(ind2,which(thenet[,2]==isallnas[a]))}
      inds<-unique(c(ind1,ind2))
      if(length(inds)>0){thenet<-thenet[-inds,]}
    }
  }
  
  ## Calculate overlaps ##
  ## @@ @@ @@  @@ @@ @@ ##
  print(paste("Calculating significance of the overlap in ",icdcode," - ",comparison,":",sep=""))
  if(length(nodesthatshouldbethere_epidemiology)>0){
    allwest<-graph_from_data_frame(west,directed = FALSE,vertices = unique(c(west[,1],west[,2],nodesthatshouldbethere_epidemiology)))
  }
  if(length(nodesthatshouldbethere_epidemiology)==0){
    allwest<-graph_from_data_frame(west,directed = FALSE)
  }
  ## Calculate the overlap between both networks ##
  alllist<-list()
  tamanos<-c()
  overlappingepidem<-list()
  overlappingepidem[["Epidemiology"]]<-E(allwest)
  for(a in 3:length(thenet[1,])){
    # a<-3
    netunderanalysis<-thenet[which(as.numeric(thenet[,a])!=0),c(1,2,a)]
    ## Positive interactions ##
    if(length(which(netunderanalysis[,3]=="1"))>0){
      pos<-netunderanalysis[which(netunderanalysis[,3]=="1"),1:2]
      ## Are we missing any node because they didn't have enough sDEGs for some of the comparisons? ##
      missingnodes<-setdiff(unique(c(thenet[,1],thenet[,2])),unique(c(pos[,1],pos[,2])))
      if(length(missingnodes)==0){posi<-graph_from_data_frame(pos,directed = FALSE)}
      if(length(missingnodes)>0){posi<-graph_from_data_frame(pos,directed = FALSE,vertices = unique(c(thenet[,1],thenet[,2])))}
      ## Overlapping edges - Positive ##
      foverlapingpositiveedges<-E(intersection(allwest,posi))
      ## Number of overlapping edges - Positive ##
      overlapallp<-length(foverlapingpositiveedges)
      ## Save overlapping edges - Positive ##
      overlappingepidem[[colnames(thenet)[a]]]$Positive<-foverlapingpositiveedges
      ## Save all the positive transcriptomic similarities ##
      overlappingepidem[[colnames(thenet)[a]]]$AllPositive<-E(posi)
      ## Random similarities ##
      overlapallp2<-c()
      for(b in 1:numberrandom){
        posi2<-rewire(posi, with=keeping_degseq(loops=FALSE, niter=ecount(posi)*10))
        overlapallp2<-c(overlapallp2,length(E(intersection(allwest,posi2))))
      }
      alllist[[colnames(thenet)[a]]]$Positive$Original<-overlapallp
      alllist[[colnames(thenet)[a]]]$Positive$Random<-overlapallp2
    }
    ## Negative interactions ##
    if(length(which(netunderanalysis[,3]=="-1"))>0){
      neg<-netunderanalysis[which(netunderanalysis[,3]=="-1"),1:2]
      ## Are we missing any node because they didn't have enough sDEGs for some of the comparisons? ##
      missingnodes<-setdiff(unique(c(thenet[,1],thenet[,2])),unique(c(neg[,1],neg[,2])))
      if(length(missingnodes)==0){negi<-graph_from_data_frame(neg,directed = FALSE)}
      if(length(missingnodes)>0){negi<-graph_from_data_frame(neg,directed = FALSE,vertices = unique(c(thenet[,1],thenet[,2])))}
      ## Overlapping edges - Negative ##
      foverlapingnegativeedges<-E(intersection(allwest,negi))
      ## Number of overlapping edges - Negative ##
      overlapalln<-length(E(intersection(allwest,negi)))
      ## Save overlapping edges - Negative ##
      overlappingepidem[[colnames(thenet)[a]]]$Negative<-foverlapingnegativeedges
      ## Save all the negative transcriptomic similarities ##
      overlappingepidem[[colnames(thenet)[a]]]$AllNegative<-E(negi)
      ## Random similarities ##
      overlapalln2<-c()
      for(b in 1:numberrandom){
        negi2<-rewire(negi, with=keeping_degseq(loops=FALSE, niter=ecount(negi)*10))
        overlapalln2<-c(overlapalln2,length(E(intersection(allwest,negi2))))
      }
      alllist[[colnames(thenet)[a]]]$Negative$Original<-overlapalln
      alllist[[colnames(thenet)[a]]]$Negative$Random<-overlapalln2
    }
    tamanos<-rbind(tamanos,c(length(E(posi)),length(E(negi))))
    print(paste(round(((a-2)/(length(thenet[1,])-2))*100,2),"%",sep=""))
  }
  ## Obtain the values of the overlap ##
  allrest<-c()
  for(a in 1:length(names(alllist))){
    pvalpos<-length(which(alllist[[a]]$Positive$Random>=alllist[[a]]$Positive$Original))/length(alllist[[a]]$Positive$Random)
    pvalneg<-length(which(alllist[[a]]$Negative$Random>=alllist[[a]]$Negative$Original))/length(alllist[[a]]$Negative$Random)
    allrest<-rbind(allrest,c(alllist[[a]]$Positive$Original,pvalpos,alllist[[a]]$Negative$Original,pvalneg,tamanos[a,],length(E(allwest))))
  }
  ## Put it beautiful ##
  if(length(allrest[1,])==7){
    allrest2<-cbind(allrest[,c(5,1)],round((allrest[,1]/allrest[,7])*100,2),allrest[,2],round((allrest[,1]/allrest[,5])*100,2),
                    allrest[,c(6,3)],round((allrest[,3]/allrest[,7])*100,2),allrest[,4],round((allrest[,3]/allrest[,6])*100,2))
    colnames(allrest2)<-c("total_pos","overlap_pos","per_epidem_pos","pval_pos","per_molec_pos","total_neg","overlap_neg","percentage_neg","pval_neg","per_molec_neg")
    rownames(allrest2)<-colnames(thenet)[3:length(colnames(thenet))]
    
    # if(length(allrest2[,1])==19){rownames(allrest2)<-colnames(thenet)[3:length(colnames(thenet))]}
    # if(length(allrest2[,1])==13){rownames(allrest2)<-colnames(thenet)[3:length(colnames(thenet))]}
    # if(length(allrest2[,1])==10){rownames(allrest2)<-colnames(thenet)[3:length(colnames(thenet))]}
    # if(length(allrest2[,1])==7){rownames(allrest2)<-colnames(thenet)[3:length(colnames(thenet))]}
    ## Order makes sense ##
    allrest<-allrest2[c(grep("All$",rownames(allrest2)),grep("Uni$",rownames(allrest2)),grep("Int$",rownames(allrest2)),7),]
  }
  
  resultados<-list("tablaall"=allrest,"OverlappingInteractions"=overlappingepidem)
  # resultados<-list("OverlappingInteractions"=overlappingepidem)
  return(resultados)
}

#### (1) Calculate overlaps with epidemiology for each of the comparisons (women, men, and adjusting for sex) ####
if(parte=="calculateoverlaps"){
  ## ICD9 Women ##
  womennet2<-read.csv2("Microarrays/Generated_Networks/ICD9_Women_binarized_network.txt",stringsAsFactors = F,sep="\t")
  womennet<-c() ; for(a in 1:length(womennet2[,1])){womennet<-rbind(womennet,c(sort(as.numeric(womennet2[a,1:2])),as.numeric(womennet2[a,3:dim(womennet2)[2]])))}
  colnames(womennet)<-colnames(womennet2)
  womennet[,1]<-as.character(womennet[,1]) ; womennet[,2]<-as.character(womennet[,2])
  icd9women<-calculateoverlaps("Microarrays","ICD9","Women",womennet,10000)
  write.table(icd9women$tablaall,"Microarrays/Results/Overlaps_ICD9_Women.txt",quote=F,sep="\t")
  save(icd9women,file="Microarrays/Results/Overlapping_interactions_ICD9_Women.RData")

  ## ICD9 Men ##
  mennet2<-read.csv2("Microarrays/Generated_Networks/ICD9_Men_binarized_network.txt",stringsAsFactors = F,sep="\t")
  mennet<-c() ; for(a in 1:length(mennet2[,1])){mennet<-rbind(mennet,c(sort(as.numeric(mennet2[a,1:2])),as.numeric(mennet2[a,3:dim(mennet2)[2]])))}
  colnames(mennet)<-colnames(mennet2)
  mennet[,1]<-as.character(mennet[,1]) ; mennet[,2]<-as.character(mennet[,2])
  icd9men<-calculateoverlaps("Microarrays","ICD9","men",mennet,10000)
  write.table(icd9men$tablaall,"Microarrays/Results/Overlaps_ICD9_Men.txt",quote=F,sep="\t")
  save(icd9men,file="Microarrays/Results/Overlapping_interactions_ICD9_Men.RData")
  print("Men ICD9 finished!")

  ## ICD9 adjusted ##
  adjustednet2<-read.csv2("Microarrays/Generated_Networks/ICD9_Adjusted_binarized_network.txt",stringsAsFactors = F,sep="\t")
  adjustednet<-c() ; for(a in 1:length(adjustednet2[,1])){adjustednet<-rbind(adjustednet,c(sort(as.numeric(adjustednet2[a,1:2])),as.numeric(adjustednet2[a,3:dim(adjustednet2)[2]])))}
  colnames(adjustednet)<-colnames(adjustednet2)
  adjustednet[,1]<-as.character(adjustednet[,1]) ; adjustednet[,2]<-as.character(adjustednet[,2])
  icd9adjusted<-calculateoverlaps("Microarrays","ICD9","adjusted",adjustednet,10000)
  write.table(icd9adjusted$tablaall,"Microarrays/Results/Overlaps_ICD9_Adjusted.txt",quote=F,sep="\t")
  save(icd9adjusted,file="Microarrays/Results/Overlapping_interactions_ICD9_Adjusted.RData")
  print("Adjusted ICD9 finished!")
  
  ## ICD10 Women ##
  womennet2<-read.csv2("Microarrays/Generated_Networks/ICD10_Women_binarized_network.txt",stringsAsFactors = F,sep="\t")
  womennet<-c() ; for(a in 1:length(womennet2[,1])){womennet<-rbind(womennet,c(sort(as.character(womennet2[a,1:2])),as.character(womennet2[a,3:dim(womennet2)[2]])))}
  colnames(womennet)<-colnames(womennet2)
  icd10women<-calculateoverlaps("Microarrays","ICD10","Women",womennet,10000)
  write.table(icd10women$tablaall,"Microarrays/Results/Overlaps_ICD10_Women.txt",quote=F,sep="\t")
  save(icd10women,file="Microarrays/Results/Overlapping_interactions_ICD10_Women.RData")
  print("Women ICD10 finished!")
  
  ## ICD10 Men ##
  mennet2<-read.csv2("Microarrays/Generated_Networks/ICD10_Men_binarized_network.txt",stringsAsFactors = F,sep="\t")
  mennet<-c() ; for(a in 1:length(mennet2[,1])){mennet<-rbind(mennet,c(sort(as.character(mennet2[a,1:2])),as.character(mennet2[a,3:dim(mennet2)[2]])))}
  colnames(mennet)<-colnames(mennet2)
  icd10men<-calculateoverlaps("Microarrays","ICD10","men",mennet,10000)
  write.table(icd10men$tablaall,"Microarrays/Results/Overlaps_ICD10_Men.txt",quote=F,sep="\t")
  save(icd10men,file="Microarrays/Results/Overlapping_interactions_ICD10_Men.RData")
  print("Men ICD10 finished!")
  
  ## ICD10 adjusted ##
  adjustednet2<-read.csv2("Microarrays/Generated_Networks/ICD10_Adjusted_binarized_network.txt",stringsAsFactors = F,sep="\t")
  adjustednet<-c() ; for(a in 1:length(adjustednet2[,1])){adjustednet<-rbind(adjustednet,c(sort(as.character(adjustednet2[a,1:2])),as.character(adjustednet2[a,3:dim(adjustednet2)[2]])))}
  colnames(adjustednet)<-colnames(adjustednet2)
  icd10adjusted<-calculateoverlaps("Microarrays","ICD10","adjusted",adjustednet,10000)
  write.table(icd10adjusted$tablaall,"Microarrays/Results/Overlaps_ICD10_Adjusted.txt",quote=F,sep="\t")
  save(icd10adjusted,file="Microarrays/Results/Overlapping_interactions_ICD10_Adjusted.RData")
  print("Adjusted ICD10 finished!")
}






