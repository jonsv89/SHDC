## CODE FOR ANALYZING THE SIGNIFICANTLY ENRICHED PATHWAYS IDENTIFIED IN THE DIFFERENT COMPARISONS ##
## Developed by Jon Sanchez-Valle
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es

## Load functions ##
library("forestplot")
library("dplyr")
library(tidyverse)
library(ggforestplot)

## Load the disease information table ##
info<-read.csv2("Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
sindupcol<-info[-(which(duplicated(info[,c(2,5)]))),c(2,5)]
colorcat9<-sindupcol[,2] ; names(colorcat9)<-gsub("ICD9_","",sindupcol[,1])
sindupcol<-info[-(which(duplicated(info[,c(6,5)]))),c(6,5)]
colorcat10<-sindupcol[,2] ; names(colorcat10)<-sindupcol[,1]
## The ICD9 names ##
info2<-info[,c(2,7)]
if(length(which(duplicated(info2)))>0){info2<-info2[-which(duplicated(info2)),]}
icd9nombres<-info2$ICD9_our_names ; names(icd9nombres)<-as.character(gsub("ICD9_","",info2$ICD9))
## The ICD10 names ##
info2<-info[,c(6,9)]
if(length(which(duplicated(info2)))>0){info2<-info2[-which(duplicated(info2)),]}
icd10nombres<-info2$ICD10_our_names ; names(icd10nombres)<-as.character(info2$ICD10)
## Colors for the pathways ##
patcatcol<-read.csv2("Data/final_colors_reactome_category.txt",stringsAsFactors = F,sep="\t")
reaccolors<-patcatcol[,2] ; names(reaccolors)<-patcatcol[,1]
## Identify the reactome parents ##
reactomeparents<-read.csv2("Data/Reactome_parents.txt",stringsAsFactors = F,sep="\t")
reactomeparents<-reactomeparents[order(reactomeparents[,2]),]
factreactomeparents<-reactomeparents[,2] ; names(factreactomeparents)<-reactomeparents[,1]
allcats<-sort(unique(reactomeparents$Parent))
totalpathspercat<-as.numeric(as.character(table(reactomeparents$Parent)[allcats]))


#### Create a table for each condition ####
comparison<-"Cases"
icdcode<-"ICD9"
experiment<-"Microarrays"
filtro<-"REACTOME"
pathwaytables<-function(comparison,icdcode,experiment,filtro,colorcat,icdnombres){
  ## Create directory ##
  if("ForestPlotsByDisease"%in%list.files(experiment)==FALSE){dir.create(paste(experiment,"/ForestPlotsByDisease",sep=""))}
  if("ForestPlotsByReactomeCategory"%in%list.files(experiment)==FALSE){dir.create(paste(experiment,"/ForestPlotsByReactomeCategory",sep=""))}
  
  files<-list.files(paste(experiment,"/",icdcode,"/",comparison,"/PathwayTables",sep=""))
  allpathnames<-c() ; allpathlist<-list() ; sexes<-c() ; genesaltered<-c()
  numeros<-read.csv2(paste(experiment,"/",icdcode,"/",comparison,"/Number_DEGs.txt",sep=""),stringsAsFactors = F,sep="\t")
  ## Save the tables into a list ##
  for(a in 1:length(files)){
    # a<-9
    tt<-read.csv2(paste(experiment,"/",icdcode,"/",comparison,"/PathwayTables/",files[a],sep=""),stringsAsFactors = F,sep="\t")
    targets<-read.csv2(paste(experiment,"/",icdcode,"/",comparison,"/Targets/",gsub("__","-",files[a]),sep=""),stringsAsFactors = F,sep="\t")
    ## If we are plotting the men-women comparison ##
    if(comparison=="Cases" || comparison=="Controls"){
      sexes<-c(sexes,paste(as.numeric(table(targets$Sex)[c("female","male")])[1],"W/",as.numeric(table(targets$Sex)[c("female","male")])[2],"M",sep=""))
    }
    numis<-numeros[which(rownames(numeros)==gsub(".txt","",gsub("__","-",files[a]))),]
    genesaltered<-c(genesaltered,paste(as.character(numis[2]),"D/",as.character(numis[3]),"U",sep=""))
    if(filtro!=""){
      if(filtro=="REACTOME"){
        tt<-tt[grep(paste("^",filtro,"_",sep=""),tt$NAME),]
        if(length(grep("REACTOME_GLYCOGEN_BREAKDOWN_GLYCOGENOLYSIS",tt$NAME))>0){
          tt<-tt[-grep("REACTOME_GLYCOGEN_BREAKDOWN_GLYCOGENOLYSIS",tt$NAME),]
        }
      }
      
    }
    allpathnames<-c(allpathnames,unique(tt$NAME))
    if(length(which(as.numeric(as.character(tt$FDR.q.val))<0.05))>0){
      ttsign<-tt[which(as.numeric(as.character(tt$FDR.q.val))<0.05),]
      signpaths<-as.numeric(as.character(ttsign$NES)) ; names(signpaths)<-as.character(ttsign$NAME)
      allpathlist[[gsub(".txt","",gsub("ICD9__","",files[a]))]]<-signpaths
    }
    if(length(which(as.numeric(as.character(tt$FDR.q.val))<0.05))==0){
      allpathlist[[gsub(".txt","",gsub("ICD9__","",files[a]))]]<-NA
    }
  }
  print("All the tables have been saved")
  allpaths<-unique(allpathnames)
  ## Plot each disease separately ##
  lfortheforest<-list()
  for(a in 1:length(names(allpathlist))){
    # a<-9
    pathways<-allpathlist[[a]]
    if(length(which(is.na(pathways)))==1){
      allcatsmeans<-cbind(allcats,0,NA,NA)
    }
    if(length(which(is.na(pathways)))==0){
      cats<-as.character(factreactomeparents[names(pathways)])
      catsnum<-as.numeric(pathways)
      allcatsmeans<-c()
      for(b in 1:length(allcats)){
        # b<-1
        cual<-which(cats==allcats[b])
        if(length(cual)==0){
          allcatsmeans<-rbind(allcatsmeans,c(allcats[b],0,NA,NA))
        }
        if(length(cual)>0){
          allcatsmeans<-rbind(allcatsmeans,c(allcats[b],length(cual),mean(catsnum[cual]),sd(catsnum[cual])))
        }
      }
    }
    ## Generate the files needed for the generation of the table ##
    nombres<-paste(allcatsmeans[,1]," (",allcatsmeans[,2],"/",totalpathspercat,")",sep="")
    fortheforest<-as.data.frame(cbind(nombres,"",allcatsmeans[,3],allcatsmeans[,4]))
    colnames(fortheforest)<-c("name","gender","mean","se")
    fortheforest$mean<-as.numeric(fortheforest$mean)
    fortheforest$se<-as.numeric(fortheforest$se)
    if(length(intersect(which(fortheforest$mean!=0),which(is.na(fortheforest$se))))>0){
      fortheforest$se[intersect(which(fortheforest$mean!=0),which(is.na(fortheforest$se)))]<-0
    }
    if(comparison=="Cases" || comparison=="Controls"){
      fortheforest$gender[which(as.numeric(fortheforest$mean)>0)]<-"women"
      fortheforest$gender[which(as.numeric(fortheforest$mean)<0)]<-"men"
      ## Add the color to women and men lines ##
      collabels<-names(table(fortheforest$gender))
      # collabels<-collabels[-which(collabels=="")]
      colvalues<-rep("#FFFFFF",length(collabels))
      if(length(which(collabels=="women"))>0){colvalues[which(collabels=="women")]<-"#C81E17"}
      if(length(which(collabels=="men"))>0){colvalues[which(collabels=="men")]<-"#3D4D8B"}
      lfortheforest[[names(allpathlist)[a]]]<-fortheforest
      thexlab<-"Women vs. Men comparison (NES)\n(positive values => up-regulated in women)\n(negative values => up-regulated in men)"
    }
    if(comparison!="Cases" && comparison!="Controls"){
      fortheforest$gender[which(as.numeric(fortheforest$mean)>0)]<-"up-regulated"
      fortheforest$gender[which(as.numeric(fortheforest$mean)<0)]<-"down-regulated"
      ## Add the color to women and men lines ##
      collabels<-names(table(fortheforest$gender))
      # collabels<-collabels[-which(collabels=="")]
      colvalues<-rep("#FFFFFF",length(collabels))
      if(length(which(collabels=="up-regulated"))>0){colvalues[which(collabels=="up-regulated")]<-"#C81E17"}
      if(length(which(collabels=="down-regulated"))>0){colvalues[which(collabels=="down-regulated")]<-"#3D4D8B"}
      lfortheforest[[names(allpathlist)[a]]]<-fortheforest
      thexlab<-"Case vs. Control comparison (NES)\n(positive values => up-regulated in cases)\n(negative values => down-regulated in cases)"
    }
    ## Plot ##
    p<-forestplot(
      df = fortheforest,
      estimate = mean,
      logodds = FALSE,
      colour=gender,
      title = paste("ICD (",names(allpathlist)[a],")\n(",genesaltered[a]," -- ",sexes[a],")\n",as.character(icdnombres[names(allpathlist)[a]]),sep=""),
      xlab = thexlab
    )+
      ggplot2::scale_colour_manual(values=colvalues,labels=collabels)+
      ggplot2::theme(plot.title=element_text(hjust=0.5))+
      ggplot2::theme(plot.title = element_text(size=10))+
      ggplot2::theme(axis.title.x = element_text(size = 10))
    ## Save the plot ##
    pdf(file=paste(experiment,"/ForestPlotsByDisease/",icdcode,"_",comparison,"_",names(allpathlist)[a],".pdf",sep=""))
      print(p)
    dev.off()
   }
  ## Generate a forest plot by pathway category ##
  for(a in 1:length(allcats)){
    fortheforestcat<-c()
    for(b in 1:length(names(lfortheforest))){
      fortheforestcat<-rbind(fortheforestcat,lfortheforest[[b]][a,])
    }
    fortheforestcat$name<-paste(names(lfortheforest),gsub(allcats[a],"",fortheforestcat$name,fixed = T))
    fortheforestcat<-fortheforestcat[order(as.numeric(gsub(" .+","",fortheforestcat$name))),]
    ## Add the color to women and men lines ##
    collabels<-names(table(fortheforestcat$gender))
    colvalues<-rep("#FFFFFF",length(collabels))
    
    
    ## Plot ##
    if(comparison=="Cases" || comparison=="Controls"){
      if(length(which(collabels=="women"))>0){colvalues[which(collabels=="women")]<-"#C81E17"}
      if(length(which(collabels=="men"))>0){colvalues[which(collabels=="men")]<-"#3D4D8B"}
      thexlab<-"Women vs. Men comparison (NES)\n(positive values => up-regulated in women)\n(negative values => up-regulated in men)"
    }
    if(comparison!="Cases" && comparison!="Controls"){
      if(length(which(collabels=="up-regulated"))>0){colvalues[which(collabels=="up-regulated")]<-"#C81E17"}
      if(length(which(collabels=="down-regulated"))>0){colvalues[which(collabels=="down-regulated")]<-"#3D4D8B"}
      thexlab<-"Case vs. Control comparison (NES)\n(positive values => up-regulated in cases)\n(negative values => down-regulated in cases)"
    }
    p<-forestplot(
      df = fortheforestcat,
      estimate = mean,
      logodds = FALSE,
      colour=gender,
      title = allcats[a],
      xlab = thexlab
    )+
      ggplot2::scale_colour_manual(values=colvalues,labels=collabels)+
      ggplot2::theme(plot.title=element_text(hjust=0.5))+
      ggplot2::theme(axis.text.y=element_text(color=as.character(colorcat[gsub(" .+","",fortheforestcat$name)])[length(fortheforestcat$name):1],size=5))+
      ggplot2::theme(plot.title = element_text(size=10))+
      ggplot2::theme(axis.title.x = element_text(size = 10))
    ## Save the plot ##
    pdf(file=paste(experiment,"/ForestPlotsByReactomeCategory/",icdcode,"_",comparison,"_",allcats[a],".pdf",sep=""))
      print(p)
    dev.off()
  }
  resultados<-list("allpaths"=allpaths,"signif_paths"=allpathlist)
  return(resultados)
}

## Women vs. Men ##
cases9<-pathwaytables("Cases","ICD9","Microarrays","REACTOME",colorcat9)
controls9<-pathwaytables("Controls","ICD9","Microarrays","REACTOME",colorcat9)
cases10<-pathwaytables("Cases","ICD10","Microarrays","REACTOME",colorcat10,icd10nombres)
controls10<-pathwaytables("Controls","ICD10","Microarrays","REACTOME",colorcat10,icd10nombres)
women10<-pathwaytables("Women","ICD10","Microarrays","REACTOME",colorcat10,icd10nombres)
men10<-pathwaytables("Men","ICD10","Microarrays","REACTOME",colorcat10,icd10nombres)


























