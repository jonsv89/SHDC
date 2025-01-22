## CODE FOR CALCULATING THE OVERLAPS BETWEEN TRANSCRIPTOMIC NETWORKS AND EPIDEMIOLOGY ##
## Developed by Jon Sanchez-Valle
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es


#### Compare the obtained networks with epidemiology by gender ####
## load the essential packages 
library("igraph")
library("UpSetR")
library("dendextend")
library("ggplot2")
library("gplots")
library("forestplot")
library("dplyr")
if("ResultsPathways"%in%list.files("Microarrays/")==FALSE){dir.create("Microarrays/ResultsPathways/")}

args = commandArgs(trailingOnly=TRUE)
categorycolors<-read.csv2("Microarrays/Data/final_colors_reactome_category.txt",stringsAsFactors = F,sep="\t")
categorycolors$Parent<-gsub("\\)","",gsub("\\(","",categorycolors$Parent))
categorycolors<-rbind(categorycolors,c("Entire set of genes","#000000"))
catcol<-categorycolors$Color ; names(catcol)<-categorycolors$Parent

## Function ##
calculateoverlaps<-function(experiment,icdcode,comparison,thenet,numberrandom,thepathway){
  ## Object specifying that this is the first time the function is run ##
  firstime<-"yes"
  if(icdcode=="ICD9"){
    epidemiology2<-read.csv2("Epidemiology/Epidemiology_ICD9.txt",stringsAsFactors = F,sep="\t")
    epidemiology2[,1]<-as.character(epidemiology2[,1]) ; epidemiology2[,2]<-as.character(epidemiology2[,2])
    icds<-unique(c(thenet[,1],thenet[,2]))
    ## Select the significant interactions between shared diseases ##
    ind1<-c() ; ind2<-c() ; for(a in 1:length(icds)){ind1<-c(ind1,which(epidemiology2[,1]==icds[a])) ; ind2<-c(ind2,which(epidemiology2[,2]==icds[a]))}
    epidemiology<-epidemiology2[intersect(ind1,ind2),]
    west<-epidemiology
  }
  
  if(icdcode=="ICD10"){
    epidemiology2<-read.csv2("Epidemiology/Epidemiology_ICD10.txt",stringsAsFactors = F,sep="\t")
    icds<-unique(c(thenet[,1],thenet[,2]))
    ## Select the significant interactions between shared diseases ##
    ind1<-c() ; ind2<-c() ; for(a in 1:length(icds)){ind1<-c(ind1,which(epidemiology2[,1]==icds[a])) ; ind2<-c(ind2,which(epidemiology2[,2]==icds[a]))}
    epidemiology<-epidemiology2[intersect(ind1,ind2),]
    west<-epidemiology[which(as.numeric(as.character(epidemiology[,grep(paste("^",comparison,sep=""),colnames(epidemiology),ignore.case = T)]))>1.01),] ### Changed from 1 to 1.01
    if(length(which(duplicated(west[,1:2])))>0){west<-west[-which(duplicated(west[,1:2])),]}
  }
  
  ## Calculate overlaps ##
  ## @@ @@ @@  @@ @@ @@ ##
  print(paste("Calculating significance of the overlap in ",icdcode," - ",comparison,":",sep=""))
  allwest<-graph_from_data_frame(west,directed = FALSE)
  ## Calculate the overlap between both networks ##
  alllist<-list()
  tamanos<-c()
  overlappingepidem<-list()
  overlappingepidem[["Epidemiology"]]<-E(allwest)
  if(length(grep(paste(icdcode,"_",comparison,"_",gsub("_ICD.+","",thepathway),".txt",sep=""),list.files("Microarrays/ResultsPathways/"),ignore.case = T))>0){
    firstime<-"no"
    fileexists<-list.files("Microarrays/ResultsPathways/")[grep(paste(icdcode,"_",comparison,"_",gsub("_ICD.+","",thepathway),".txt",sep=""),list.files("Microarrays/ResultsPathways/"),ignore.case = T)]
    previousmetrics<-rownames(read.csv2(paste("Microarrays/ResultsPathways/",fileexists,sep=""),stringsAsFactors = F,sep="\t"))
    ## Get the indexes of the comparisons that have been already done ## 
    indexes<-c() ; for(zz in 1:length(previousmetrics)){indexes<-c(indexes,which(colnames(thenet)==previousmetrics[zz]))}
    thenet<-thenet[,-indexes]
  }
  if(length(thenet[1,])>2){
    for(a in 3:length(thenet[1,])){
      # a<-3
      netunderanalysis<-thenet[which(as.numeric(thenet[,a])!=0),c(1,2,a)]
      if(class(netunderanalysis)[1]=="character"){
        netunderanalysis2<-cbind(netunderanalysis[1],netunderanalysis[2],netunderanalysis[3])
        colnames(netunderanalysis2)<-names(netunderanalysis)
        netunderanalysis<-netunderanalysis2
      }
      ## We don't find a single connection ##
      ## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
      if(length(which(netunderanalysis[,3]=="1"))==0 && length(which(netunderanalysis[,3]=="-1"))==0){
        ## Positive
        overlappingepidem[[colnames(thenet)[a]]]$Positive<-NA
        overlappingepidem[[colnames(thenet)[a]]]$AllPositive<-NA
        alllist[[colnames(thenet)[a]]]$Positive$Original<-0
        alllist[[colnames(thenet)[a]]]$Positive$Random<-NA
        ## Negative
        overlappingepidem[[colnames(thenet)[a]]]$Negative<-NA
        overlappingepidem[[colnames(thenet)[a]]]$AllNegative<-NA
        alllist[[colnames(thenet)[a]]]$Negative$Original<-0
        alllist[[colnames(thenet)[a]]]$Negative$Random<-NA
        ## Tamanos
        tamanos<-rbind(tamanos,c(0,0))
      }
      ## We only find positive connections ##
      ## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
      if(length(which(netunderanalysis[,3]=="1"))>0 && length(which(netunderanalysis[,3]=="-1"))==0){
        ## Positive ##
        pos<-netunderanalysis[which(netunderanalysis[,3]=="1"),1:2]
        if(length(dim(pos))==0){
          pos<-cbind(pos[1],pos[2])
          colnames(pos)<-c("Disease1","Disease2")
        }
        posi<-graph_from_data_frame(pos,directed = FALSE)
        overlapallp<-length(E(intersection(allwest,posi)))
        foverlapingpositiveedges<-E(intersection(allwest,posi))
        overlappingepidem[[colnames(thenet)[a]]]$Positive<-foverlapingpositiveedges
        overlappingepidem[[colnames(thenet)[a]]]$AllPositive<-E(posi)
        alllist[[colnames(thenet)[a]]]$Positive$Original<-overlapallp
        ## Si solo dos interacciones no se puede hacer rewire
        if(length(which(netunderanalysis[,3]=="1"))>3){
          # if(length(which(netunderanalysis[,3]=="1"))>0){
          overlapallp2<-c()
          for(b in 1:numberrandom){
            posi2<-rewire(posi, with=keeping_degseq(loops=FALSE, niter=ecount(posi)*10))
            overlapallp2<-c(overlapallp2,length(E(intersection(allwest,posi2))))
            alllist[[colnames(thenet)[a]]]$Positive$Random<-overlapallp2
          }
        }
        if(length(which(netunderanalysis[,3]=="1"))<=3){
          alllist[[colnames(thenet)[a]]]$Positive$Random<-NA
        }
        ## Negative
        overlappingepidem[[colnames(thenet)[a]]]$Negative<-NA
        overlappingepidem[[colnames(thenet)[a]]]$AllNegative<-NA
        alllist[[colnames(thenet)[a]]]$Negative$Original<-0
        alllist[[colnames(thenet)[a]]]$Negative$Random<-NA
        ## Tamanos
        tamanos<-rbind(tamanos,c(length(E(posi)),0))
      }
      ## We only find negative connections ##
      ## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
      if(length(which(netunderanalysis[,3]=="1"))==0 && length(which(netunderanalysis[,3]=="-1"))>0){
        ## Positive
        overlappingepidem[[colnames(thenet)[a]]]$Positive<-NA
        overlappingepidem[[colnames(thenet)[a]]]$AllPositive<-NA
        alllist[[colnames(thenet)[a]]]$Positive$Original<-0
        alllist[[colnames(thenet)[a]]]$Positive$Random<-NA
        ## Negative
        neg<-netunderanalysis[which(netunderanalysis[,3]=="-1"),1:2]
        if(length(dim(neg))==0){
          neg<-cbind(neg[1],neg[2])
          colnames(neg)<-c("Disease1","Disease2")
        }
        negi<-graph_from_data_frame(neg,directed = FALSE)
        overlapalln<-length(E(intersection(allwest,negi)))
        foverlapingnegativeedges<-E(intersection(allwest,negi))
        overlappingepidem[[colnames(thenet)[a]]]$Negative<-foverlapingnegativeedges
        overlappingepidem[[colnames(thenet)[a]]]$AllNegative<-E(negi)
        alllist[[colnames(thenet)[a]]]$Negative$Original<-overlapalln
        ## Si solo dos interacciones no se puede hacer rewire
        if(length(which(netunderanalysis[,3]=="-1"))>3){
          # if(length(which(netunderanalysis[,3]=="-1"))>0){
          overlapalln2<-c()
          for(b in 1:numberrandom){
            negi2<-rewire(negi, with=keeping_degseq(loops=FALSE, niter=ecount(negi)*10))
            overlapalln2<-c(overlapalln2,length(E(intersection(allwest,negi2))))
            alllist[[colnames(thenet)[a]]]$Negative$Random<-overlapalln2
          }
        }
        if(length(which(netunderanalysis[,3]=="-1"))<=3){
          alllist[[colnames(thenet)[a]]]$Negative$Random<-NA
        }
        ## Tamanos
        tamanos<-rbind(tamanos,c(0,length(E(negi))))
      }
      ## We find positive and negative connections ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
      if(length(which(netunderanalysis[,3]=="1"))>0 && length(which(netunderanalysis[,3]=="-1"))>0){
        ## Positive
        pos<-netunderanalysis[which(netunderanalysis[,3]=="1"),1:2]
        if(length(dim(pos))==0){
          pos<-cbind(pos[1],pos[2])
          colnames(pos)<-c("Disease1","Disease2")
        }
        posi<-graph_from_data_frame(pos,directed = FALSE)
        overlapallp<-length(E(intersection(allwest,posi)))
        foverlapingpositiveedges<-E(intersection(allwest,posi))
        overlappingepidem[[colnames(thenet)[a]]]$Positive<-foverlapingpositiveedges
        overlappingepidem[[colnames(thenet)[a]]]$AllPositive<-E(posi)
        alllist[[colnames(thenet)[a]]]$Positive$Original<-overlapallp
        ## Si solo dos interacciones no se puede hacer rewire
        if(length(which(netunderanalysis[,3]=="1"))>3){
          # if(length(which(netunderanalysis[,3]=="1"))>0){
          overlapallp2<-c()
          for(b in 1:numberrandom){
            posi2<-rewire(posi, with=keeping_degseq(loops=FALSE, niter=ecount(posi)*10))
            overlapallp2<-c(overlapallp2,length(E(intersection(allwest,posi2))))
            alllist[[colnames(thenet)[a]]]$Positive$Random<-overlapallp2
          }
        }
        if(length(which(netunderanalysis[,3]=="1"))<=3){
          alllist[[colnames(thenet)[a]]]$Positive$Random<-NA
        }
        ## Negative
        neg<-netunderanalysis[which(netunderanalysis[,3]=="-1"),1:2]
        if(length(dim(neg))==0){
          neg<-cbind(neg[1],neg[2])
          colnames(neg)<-c("Disease1","Disease2")
        }
        negi<-graph_from_data_frame(neg,directed = FALSE)
        overlapalln<-length(E(intersection(allwest,negi)))
        foverlapingnegativeedges<-E(intersection(allwest,negi))
        overlappingepidem[[colnames(thenet)[a]]]$Negative<-foverlapingnegativeedges
        overlappingepidem[[colnames(thenet)[a]]]$AllNegative<-E(negi)
        alllist[[colnames(thenet)[a]]]$Negative$Original<-overlapalln
        ## Si solo dos interacciones no se puede hacer rewire
        if(length(which(netunderanalysis[,3]=="-1"))>3){
          # if(length(which(netunderanalysis[,3]=="-1"))>0){
          overlapalln2<-c()
          for(b in 1:numberrandom){
            negi2<-rewire(negi, with=keeping_degseq(loops=FALSE, niter=ecount(negi)*10))
            overlapalln2<-c(overlapalln2,length(E(intersection(allwest,negi2))))
            alllist[[colnames(thenet)[a]]]$Negative$Random<-overlapalln2
          }
        }
        if(length(which(netunderanalysis[,3]=="-1"))<=3){
          alllist[[colnames(thenet)[a]]]$Negative$Random<-NA
        }
        ## Tamanos
        tamanos<-rbind(tamanos,c(length(E(posi)),length(E(negi))))
      }
      print(paste(round(((a-2)/(length(thenet[1,])-2))*100,2),"%",sep=""))
    }
    ## Obtain the values of the overlap ##
    allrest<-c()
    for(a in 1:length(names(alllist))){
      if(length(alllist[[a]]$Positive$Random)==1){pvalpos<-NA}
      if(length(alllist[[a]]$Positive$Random)>1){pvalpos<-length(which(alllist[[a]]$Positive$Random>=alllist[[a]]$Positive$Original))/length(alllist[[a]]$Positive$Random)}
      if(length(alllist[[a]]$Negative$Random)==1){pvalneg<-NA}
      if(length(alllist[[a]]$Negative$Random)>1){pvalneg<-length(which(alllist[[a]]$Negative$Random>=alllist[[a]]$Negative$Original))/length(alllist[[a]]$Negative$Random)}
      
      allrest<-rbind(allrest,c(alllist[[a]]$Positive$Original,pvalpos,alllist[[a]]$Negative$Original,pvalneg,tamanos[a,],length(E(allwest))))
    }
    ## Put it beautiful ##
    if(length(allrest[1,])==7){
      allrest2<-cbind(allrest[,c(5,1)],round((allrest[,1]/allrest[,7])*100,2),allrest[,2],round((allrest[,1]/allrest[,5])*100,2),
                      allrest[,c(6,3)],round((allrest[,3]/allrest[,7])*100,2),allrest[,4],round((allrest[,3]/allrest[,6])*100,2))
      colnames(allrest2)<-c("total_pos","overlap_pos","per_epidem_pos","pval_pos","per_molec_pos","total_neg","overlap_neg","percentage_neg","pval_neg","per_molec_neg")
      rownames(allrest2)<-colnames(thenet)[3:length(colnames(thenet))]
      ## Order makes sense ##
      if(firstime=="yes"){
        allrest<-allrest2[c(grep("All$",rownames(allrest2)),grep("Uni$",rownames(allrest2)),grep("Int$",rownames(allrest2)),grep("Fisher",rownames(allrest2))),]
      }
      if(firstime=="no"){
        ## Integrate both tables
        previousfile<-read.csv2(paste("Microarrays/ResultsPathways/",fileexists,sep=""),stringsAsFactors = F,sep="\t")
        allrest2<-rbind(previousfile,allrest2)
        allrest<-allrest2[c(grep("All$",rownames(allrest2)),grep("Uni$",rownames(allrest2)),grep("Int$",rownames(allrest2)),grep("Fisher",rownames(allrest2))),]
        ## Integrate both objects 
        ## Adjusted ICD9
        if(icdcode=="ICD9" && comparison=="adjusted"){
          load(paste("Microarrays/ResultsPathways/",gsub("Overlaps","Overlapping_interactions",gsub(".txt",".RData",fileexists)),sep=""))
          for(zz in 2:length(names(overlappingepidem))){
            icd9adjusted$OverlappingInteractions[[names(overlappingepidem)[zz]]]<-overlappingepidem[[zz]]
          }
          overlappingepidem<-icd9adjusted$OverlappingInteractions
        }
        ## Women ICD9
        if(icdcode=="ICD9" && comparison=="women"){
          load(paste("Microarrays/ResultsPathways/",gsub("Overlaps","Overlapping_interactions",gsub(".txt",".RData",fileexists)),sep=""))
          for(zz in 2:length(names(overlappingepidem))){
            icd9women$OverlappingInteractions[[names(overlappingepidem)[zz]]]<-overlappingepidem[[zz]]
          }
          overlappingepidem<-icd9women$OverlappingInteractions
        }
        ## Men ICD9
        if(icdcode=="ICD9" && comparison=="men"){
          load(paste("Microarrays/ResultsPathways/",gsub("Overlaps","Overlapping_interactions",gsub(".txt",".RData",fileexists)),sep=""))
          for(zz in 2:length(names(overlappingepidem))){
            icd9men$OverlappingInteractions[[names(overlappingepidem)[zz]]]<-overlappingepidem[[zz]]
          }
          overlappingepidem<-icd9men$OverlappingInteractions
        }
        ## Adjusted ICD10
        if(icdcode=="ICD10" && comparison=="adjusted"){
          load(paste("Microarrays/ResultsPathways/",gsub("Overlaps","Overlapping_interactions",gsub(".txt",".RData",fileexists)),sep=""))
          for(zz in 2:length(names(overlappingepidem))){
            icd10adjusted$OverlappingInteractions[[names(overlappingepidem)[zz]]]<-overlappingepidem[[zz]]
          }
          overlappingepidem<-icd10adjusted$OverlappingInteractions
        }
        ## Women ICD10
        if(icdcode=="ICD10" && comparison=="women"){
          load(paste("Microarrays/ResultsPathways/",gsub("Overlaps","Overlapping_interactions",gsub(".txt",".RData",fileexists)),sep=""))
          for(zz in 2:length(names(overlappingepidem))){
            icd10women$OverlappingInteractions[[names(overlappingepidem)[zz]]]<-overlappingepidem[[zz]]
          }
          overlappingepidem<-icd10women$OverlappingInteractions
        }
        ## Men ICD10
        if(icdcode=="ICD10" && comparison=="men"){
          load(paste("Microarrays/ResultsPathways/",gsub("Overlaps","Overlapping_interactions",gsub(".txt",".RData",fileexists)),sep=""))
          for(zz in 2:length(names(overlappingepidem))){
            icd10men$OverlappingInteractions[[names(overlappingepidem)[zz]]]<-overlappingepidem[[zz]]
          }
          overlappingepidem<-icd10men$OverlappingInteractions
        }
      }
    }
    resultados<-list("tablaall"=allrest,"OverlappingInteractions"=overlappingepidem)
    # resultados<-list("OverlappingInteractions"=overlappingepidem)
    return(resultados)
  }
  if(length(thenet[1,])==2){print("Pathway already analyzed!")}
}

## ICD9 adjusted ##
## @@ @@ @ @@ @@ ##
if(args[1]=="adjusted9"){
  adj9<-list.files("Microarrays/Generated_Networks_Pathways/")[grep("ICD9_Adjusted_binarized",list.files("Microarrays/Generated_Networks_Pathways/"))]
  for(tt in 1:length(adj9)){
    print(paste(tt,adj9[tt],sep=" - "))
    adjustednet2<-read.csv2(paste("Microarrays/Generated_Networks_Pathways/",adj9[tt],sep=""),stringsAsFactors = F,sep="\t")
    adjustednet<-c() ; for(a in 1:length(adjustednet2[,1])){adjustednet<-rbind(adjustednet,c(sort(as.numeric(adjustednet2[a,1:2])),as.numeric(adjustednet2[a,3:dim(adjustednet2)[2]])))}
    colnames(adjustednet)<-colnames(adjustednet2)
    adjustednet[,1]<-as.character(adjustednet[,1]) ; adjustednet[,2]<-as.character(adjustednet[,2])
    icd9adjusted<-calculateoverlaps("Microarrays","ICD9","adjusted",adjustednet,10000,adj9[tt])
    write.table(icd9adjusted$tablaall,paste("Microarrays/ResultsPathways/Overlaps_ICD9_Adjusted_",gsub("_ICD9_Adjusted_binarized_network.txt","",adj9[tt]),".txt",sep=""),quote=F,sep="\t")
    save(icd9adjusted,file=paste("Microarrays/ResultsPathways/Overlapping_interactions_ICD9_Adjusted_",gsub("_ICD9_Adjusted_binarized_network.txt","",adj9[tt]),".RData",sep=""))
  }
}

## ICD9 Women ##
## @@ @@@@ @@ ##
if(args[1]=="women9"){
  wom9<-list.files("Microarrays/Generated_Networks_Pathways/")[grep("ICD9_Women_binarized",list.files("Microarrays/Generated_Networks_Pathways/"))]
  for(tt in 1:length(wom9)){
    print(paste(tt,wom9[tt],sep=" - "))
    womennet2<-read.csv2(paste("Microarrays/Generated_Networks_Pathways/",wom9[tt],sep=""),stringsAsFactors = F,sep="\t")
    womennet<-c() ; for(a in 1:length(womennet2[,1])){womennet<-rbind(womennet,c(sort(as.numeric(womennet2[a,1:2])),as.numeric(womennet2[a,3:dim(womennet2)[2]])))}
    colnames(womennet)<-colnames(womennet2)
    womennet[,1]<-as.character(womennet[,1]) ; womennet[,2]<-as.character(womennet[,2])
    icd9women<-calculateoverlaps("Microarrays","ICD9","women",womennet,10000,wom9[tt])
    write.table(icd9women$tablaall,paste("Microarrays/ResultsPathways/Overlaps_ICD9_Women_",gsub("_ICD9_Women_binarized_network.txt","",wom9[tt]),".txt",sep=""),quote=F,sep="\t")
    save(icd9women,file=paste("Microarrays/ResultsPathways/Overlapping_interactions_ICD9_Women_",gsub("_ICD9_Women_binarized_network.txt","",wom9[tt]),".RData",sep=""))
  }
}

## ICD9 Men ##
## @@ @@ @@ ##
if(args[1]=="men9"){
  men9<-list.files("Microarrays/Generated_Networks_Pathways/")[grep("ICD9_Men_binarized",list.files("Microarrays/Generated_Networks_Pathways/"))]
  for(tt in 1:length(men9)){
    print(paste(tt,men9[tt],sep=" - "))
    mennet2<-read.csv2(paste("Microarrays/Generated_Networks_Pathways/",men9[tt],sep=""),stringsAsFactors = F,sep="\t")
    mennet<-c() ; for(a in 1:length(mennet2[,1])){mennet<-rbind(mennet,c(sort(as.numeric(mennet2[a,1:2])),as.numeric(mennet2[a,3:dim(mennet2)[2]])))}
    colnames(mennet)<-colnames(mennet2)
    mennet[,1]<-as.character(mennet[,1]) ; mennet[,2]<-as.character(mennet[,2])
    icd9men<-calculateoverlaps("Microarrays","ICD9","men",mennet,10000,men9[tt])
    write.table(icd9men$tablaall,paste("Microarrays/ResultsPathways/Overlaps_ICD9_Men_",gsub("_ICD9_Men_binarized_network.txt","",men9[tt]),".txt",sep=""),quote=F,sep="\t")
    save(icd9men,file=paste("Microarrays/ResultsPathways/Overlapping_interactions_ICD9_Men_",gsub("_ICD9_Men_binarized_network.txt","",men9[tt]),".RData",sep=""))
  }
}

## ICD10 adjusted ##
## @@ @@ @ @@ @@ ##
if(args[1]=="adjusted10"){
  adj10<-list.files("Microarrays/Generated_Networks_Pathways/")[grep("ICD10_Adjusted_binarized",list.files("Microarrays/Generated_Networks_Pathways/"))]
  for(tt in 1:length(adj10)){
    print(paste(tt,adj10[tt],sep=" - "))
    adjustednet2<-read.csv2(paste("Microarrays/Generated_Networks_Pathways/",adj10[tt],sep=""),stringsAsFactors = F,sep="\t")
    adjustednet<-c() ; for(a in 1:length(adjustednet2[,1])){adjustednet<-rbind(adjustednet,c(sort(as.character(adjustednet2[a,1:2])),as.numeric(adjustednet2[a,3:dim(adjustednet2)[2]])))}
    colnames(adjustednet)<-colnames(adjustednet2)
    icd10adjusted<-calculateoverlaps("Microarrays","ICD10","adjusted",adjustednet,10000,adj10[tt])
    write.table(icd10adjusted$tablaall,paste("Microarrays/ResultsPathways/Overlaps_ICD10_Adjusted_",gsub("_ICD10_Adjusted_binarized_network.txt","",adj10[tt]),".txt",sep=""),quote=F,sep="\t")
    save(icd10adjusted,file=paste("Microarrays/ResultsPathways/Overlapping_interactions_ICD10_Adjusted_",gsub("_ICD10_Adjusted_binarized_network.txt","",adj10[tt]),".RData",sep=""))
  }
}

## ICD10 Women ##
## @@ @@@@ @@ ##
if(args[1]=="women10"){
  wom10<-list.files("Microarrays/Generated_Networks_Pathways/")[grep("ICD10_Women_binarized",list.files("Microarrays/Generated_Networks_Pathways/"))]
  for(tt in 1:length(wom10)){
    print(paste(tt,wom10[tt],sep=" - "))
    womennet2<-read.csv2(paste("Microarrays/Generated_Networks_Pathways/",wom10[tt],sep=""),stringsAsFactors = F,sep="\t")
    womennet<-c() ; for(a in 1:length(womennet2[,1])){womennet<-rbind(womennet,c(sort(as.character(womennet2[a,1:2])),as.numeric(womennet2[a,3:dim(womennet2)[2]])))}
    colnames(womennet)<-colnames(womennet2)
    icd10women<-calculateoverlaps("Microarrays","ICD10","women",womennet,10000,wom10[tt])
    write.table(icd10women$tablaall,paste("Microarrays/ResultsPathways/Overlaps_ICD10_Women_",gsub("_ICD10_Women_binarized_network.txt","",wom10[tt]),".txt",sep=""),quote=F,sep="\t")
    save(icd10women,file=paste("Microarrays/ResultsPathways/Overlapping_interactions_ICD10_Women_",gsub("_ICD10_Women_binarized_network.txt","",wom10[tt]),".RData",sep=""))
  }
}

## ICD10 Men ##
## @@ @@ @@ ##
if(args[1]=="men10"){
  men10<-list.files("Microarrays/Generated_Networks_Pathways/")[grep("ICD10_Men_binarized",list.files("Microarrays/Generated_Networks_Pathways/"))]
  for(tt in 1:length(men10)){
    print(paste(tt,men10[tt],sep=" - "))
    mennet2<-read.csv2(paste("Microarrays/Generated_Networks_Pathways/",men10[tt],sep=""),stringsAsFactors = F,sep="\t")
    mennet<-c() ; for(a in 1:length(mennet2[,1])){mennet<-rbind(mennet,c(sort(as.character(mennet2[a,1:2])),as.numeric(mennet2[a,3:dim(mennet2)[2]])))}
    colnames(mennet)<-colnames(mennet2)
    icd10men<-calculateoverlaps("Microarrays","ICD10","men",mennet,10000,men10[tt])
    write.table(icd10men$tablaall,paste("Microarrays/ResultsPathways/Overlaps_ICD10_Men_",gsub("_ICD10_Men_binarized_network.txt","",men10[tt]),".txt",sep=""),quote=F,sep="\t")
    save(icd10men,file=paste("Microarrays/ResultsPathways/Overlapping_interactions_ICD10_Men_",gsub("_ICD10_Men_binarized_network.txt","",men10[tt]),".RData",sep=""))
  }
}

## Lets put the names of the distances in the same way in the original networks and the pathway networks ##
if(args[1]=="makenetscomparablebyname"){
  rdatas<-list.files("Microarrays/ResultsPathways/")[grep("RData",list.files("Microarrays/ResultsPathways/"))]
  for(a in 1:length(rdatas)){
    # a<-1
    load(paste("Microarrays/ResultsPathways/",rdatas[a],sep=""))
    ## Adjusted ##
    if(length(grep("_Adjusted_",rdatas[a]))>0 && length(grep("ICD10",rdatas[a]))>0){
      if(length(grep("EuclideanAll",names(icd10adjusted$OverlappingInteractions)))>0){
        names(icd10adjusted$OverlappingInteractions)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",names(icd10adjusted$OverlappingInteractions))))
        rownames(icd10adjusted$tablaall)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",rownames(icd10adjusted$tablaall))))
        save(icd10adjusted,file=paste("Microarrays/ResultsPathways/",rdatas[a],sep=""))
      }
    }
    if(length(grep("_Adjusted_",rdatas[a]))>0 && length(grep("ICD9",rdatas[a]))>0){
      if(length(grep("EuclideanAll",names(icd9adjusted$OverlappingInteractions)))>0){
        names(icd9adjusted$OverlappingInteractions)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",names(icd9adjusted$OverlappingInteractions))))
        rownames(icd9adjusted$tablaall)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",rownames(icd9adjusted$tablaall))))
        save(icd9adjusted,file=paste("Microarrays/ResultsPathways/",rdatas[a],sep=""))
      }
    }
    ## Women ##
    if(length(grep("_Women_",rdatas[a]))>0 && length(grep("ICD10",rdatas[a]))>0){
      if(length(grep("EuclideanAll",names(icd10women$OverlappingInteractions)))>0){
        names(icd10women$OverlappingInteractions)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",names(icd10women$OverlappingInteractions))))
        rownames(icd10women$tablaall)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",rownames(icd10women$tablaall))))
        save(icd10women,file=paste("Microarrays/ResultsPathways/",rdatas[a],sep=""))
      }
    }
    if(length(grep("_Women_",rdatas[a]))>0 && length(grep("ICD9",rdatas[a]))>0){
      if(length(grep("EuclideanAll",names(icd9women$OverlappingInteractions)))>0){
        names(icd9women$OverlappingInteractions)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",names(icd9women$OverlappingInteractions))))
        rownames(icd9women$tablaall)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",rownames(icd9women$tablaall))))
        save(icd9women,file=paste("Microarrays/ResultsPathways/",rdatas[a],sep=""))
      }
    }
    ## Men ##
    if(length(grep("_Men_",rdatas[a]))>0 && length(grep("ICD10",rdatas[a]))>0){
      if(length(grep("EuclideanAll",names(icd10men$OverlappingInteractions)))>0){
        names(icd10men$OverlappingInteractions)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",names(icd10men$OverlappingInteractions))))
        rownames(icd10men$tablaall)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",rownames(icd10men$tablaall))))
        save(icd10men,file=paste("Microarrays/ResultsPathways/",rdatas[a],sep=""))
      }
    }
    if(length(grep("_Men_",rdatas[a]))>0 && length(grep("ICD9",rdatas[a]))>0){
      if(length(grep("EuclideanAll",names(icd9men$OverlappingInteractions)))>0){
        names(icd9men$OverlappingInteractions)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",names(icd9men$OverlappingInteractions))))
        rownames(icd9men$tablaall)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",rownames(icd9men$tablaall))))
        save(icd9men,file=paste("Microarrays/ResultsPathways/",rdatas[a],sep=""))
      }
    }
    print(paste(round((a/length(rdatas))*100,2),"%",sep=""))
  }
  print(".RDatas corrected!!")
  txts<-list.files("Microarrays/ResultsPathways/")[grep("txt",list.files("Microarrays/ResultsPathways/"))]
  for(a in 1:length(txts)){
    # a<-1
    tabla<-read.csv2(paste("Microarrays/ResultsPathways/",txts[a],sep=""),stringsAsFactors = F,sep="\t")
    if(length(grep("EuclideanAll",rownames(tabla)))>0){
      rownames(tabla)<-gsub("Canberra","CanberraDistance",gsub("Manhattan","ManhattanDistance",gsub("Euclidean","EuclideanDistance",rownames(tabla))))
      write.table(tabla,paste("Microarrays/ResultsPathways/",txts[a],sep=""),quote=F,sep="\t")
    }
    print(paste(round((a/length(txts))*100,2),"%",sep=""))
  }
  print(".txts corrected!!")
}




















