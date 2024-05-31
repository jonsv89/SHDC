## CODE FOR THE GENERATION OF THE DISEASE NETWORKS FOR MEN, WOMEN, AND ADJUSTING BY SEX BASED ON GENES ##
## Developed by Jon Sanchez-Valle & Beatriz Urda-Garcia
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es @ beatriz.urda@bsc.es

args = commandArgs(trailingOnly=TRUE)
library(lsa)
library("geometry")

#### Add the estimate of cosine similarities ####
## @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ ##
if("OtherMetrics"%in%list.files("Microarrays")==FALSE){dir.create("Microarrays/OtherMetrics")}
if("CosineSimilarities"%in%list.files("Microarrays/OtherMetrics")==FALSE){dir.create("Microarrays/OtherMetrics/CosineSimilarities")}
if("EuclideanDistance"%in%list.files("Microarrays/OtherMetrics")==FALSE){dir.create("Microarrays/OtherMetrics/EuclideanDistance")}
if("ManhattanDistance"%in%list.files("Microarrays/OtherMetrics")==FALSE){dir.create("Microarrays/OtherMetrics/ManhattanDistance")}
if("CanberraDistance"%in%list.files("Microarrays/OtherMetrics")==FALSE){dir.create("Microarrays/OtherMetrics/CanberraDistance")}

# nameout<-"Men"
# classification<-"ICD9"
# experiment<-"Microarrays"
# set.seed(1)
build_cosinesimilarity_network<-function(nameout,classification,experiment){
  ## Start the analysis!! ##
  ## @@ @@ @@ @@ @@ @@ @@ ##
  files<-list.files(paste(experiment,"/",classification,"/Comparisons/",nameout,"/DifferentialExpressions/",sep=""))
  ## Run the loop ##
  sdegs<-list() ; vsdegs<-list() ; vnum<-c() ; vname<-c()
  allgenes<-list()
  for(a in 1:length(files)){
    # a<-1
    tabla<-read.csv2(paste(experiment,"/",classification,"/Comparisons/",nameout,"/DifferentialExpressions/",files[a],sep=""),stringsAsFactors = F,sep="\t")
    ups<-rownames(tabla)[intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))>0))]
    upslogs<-as.numeric(as.character(tabla$logFC[intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))>0))]))
    downs<-rownames(tabla)[intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))<0))]
    downslogs<-as.numeric(as.character(tabla$logFC[intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))<0))]))
    names(downslogs)<-downs
    names(upslogs)<-ups
    sdegs[[gsub(".txt","",files[a])]]$Up<-ups
    sdegs[[gsub(".txt","",files[a])]]$Down<-downs
    vsdegs[[gsub(".txt","",files[a])]]$Up<-upslogs
    vsdegs[[gsub(".txt","",files[a])]]$Down<-downslogs
    vnum<-rbind(vnum,c(length(ups),length(downs)))
    vname<-c(vname,gsub(".+_","",gsub(".txt","",files[a])))
    todosgenes<-tabla$logFC ; names(todosgenes)<-rownames(tabla)
    allgenes[[gsub(".txt","",files[a])]]<-todosgenes
    print(paste(round((a/length(files))*100,2),"%",sep=""))
  }
  print("Vectors generated")
  
  ## Similarities between pairs of diseases ##
  netall<-c()
  for(a in 1:(length(names(vsdegs))-1)){
    for(b in (a+1):length(names(vsdegs))){
      ## Over all the genes ##
      ## @@ @@ @@  @@ @@ @@ ##
      commongenes<-intersect(names(allgenes[[a]]),names(allgenes[[b]]))
      avec1<-as.numeric(as.character(allgenes[[a]][commongenes]))
      avec2<-as.numeric(as.character(allgenes[[b]][commongenes]))
      acossim<-cosine(avec1,avec2)[1,1]
      rancossim<-numeric(length=10000)
      for(ran in 1:10000){
        rancossim[ran]<-cosine(sample(avec1),avec2)[1,1]
      }
      if(acossim<0){
        apval<-length(which(rancossim<acossim))/10000
      }
      if(acossim>0){
        apval<-length(which(rancossim>acossim))/10000
      }
      
      ## Over the intersection of the sDEGs ##
      ## @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ ##
      unof<-vsdegs[[a]] ; dosf<-vsdegs[[b]]
      allunof<-c(unof$Up,unof$Down) ; alldosf<-c(dosf$Up,dosf$Down)
      fintersection<-intersect(names(allunof),names(alldosf))
      if(length(fintersection)<=3){
        icossim<-NA
        ipval<-NA
      }
      if(length(fintersection)>3){
        ivec1<-as.numeric(as.character(allunof[fintersection]))
        ivec2<-as.numeric(as.character(alldosf[fintersection]))
        icossim<-cosine(ivec1,ivec2)[1,1]
        irancossim<-numeric(length=10000)
        for(ran in 1:10000){
          irancossim[ran]<-cosine(sample(ivec1),ivec2)[1,1]
        }
        if(icossim<0){
          ipval<-length(which(irancossim<icossim))/10000
        }
        if(icossim>0){
          ipval<-length(which(irancossim>icossim))/10000
        }
      }
      ## Over the union of sDEGs ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ ##
      funion<-unique(c(names(allunof),names(alldosf)))
      unioncommon<-intersect(names(which(is.na(allgenes[[a]][funion])==FALSE)),names(which(is.na(allgenes[[b]][funion])==FALSE)))
      numberifunion<-length(unioncommon)
      if(numberifunion<=3){
        ucossim<-NA
        upval<-NA
      }
      if(numberifunion>3){
        uvec1<-as.numeric(as.character(allgenes[[a]][unioncommon]))
        uvec2<-as.numeric(as.character(allgenes[[b]][unioncommon]))
        ucossim<-cosine(uvec1,uvec2)[1,1]
        urancossim<-numeric(length=10000)
        for(ran in 1:10000){
          urancossim[ran]<-cosine(sample(uvec1),uvec2)[1,1]
        }
        if(ucossim<0){
          upval<-length(which(urancossim<ucossim))/10000
        }
        if(ucossim>0){
          upval<-length(which(urancossim>ucossim))/10000
        }
      }
      ## Put all the correlations together ##
      ## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
      netall<-rbind(netall,c(names(allgenes)[c(a,b)],acossim,apval,icossim,ipval,ucossim,upval))
    }
    print(paste(round((a/length(names(allgenes)))*100,2),"%",sep=""))
  }
  ## Add column names ##
  colnames(netall)<-c("Disease1","Disease2","AllCosineSimilarity","AllPvalCosineSimilarity",
                      "IntCosineSimilarity","IntPvalCosineSimilarity",
                      "UniCosineSimilarity","UniPvalCosineSimilarity")
  ## Write the tables ##
  write.table(netall,paste(experiment,"/OtherMetrics/CosineSimilarities/",classification,"_",nameout,"_pvals_network.txt",sep=""),quote=F,sep="\t",row.names=F)
  return(netall)
}

## Run in Mare ##
if(as.numeric(args[1])==2){
  nameout<-args[2]
  classification<-args[3]
  experiment<-args[4]
  print(nameout)
  print(classification)
  print(experiment)
  salida<-build_cosinesimilarity_network(nameout,classification,experiment)
}

## Create a function to correct for multiple testing ##
# nameout<-"Women" ; classification<-"ICD9" ; experiment<-"Microarrays"
correctformultipletesting<-function(nameout,classification,experiment){
  netall1<-read.csv2(paste(experiment,"/Generated_Networks/",classification,"_",nameout,"_pvals_network.txt",sep=""),stringsAsFactors = F,sep="\t")
  if(dim(netall1)[2]==15){
    netall2<-read.csv2(paste(experiment,"/OtherMetrics/CosineSimilarities/",classification,"_",nameout,"_pvals_network.txt",sep=""),stringsAsFactors = F,sep="\t")
    ## If both tables have the same number of interactions and the interactions in the same order:
    if(dim(netall1)[1]==dim(netall2)[1] && length(which(netall1[,1]==netall2[,1]))==dim(netall1)[1] && length(which(netall1[,2]==netall2[,2]))==dim(netall2)[1]){
      netall<-cbind(netall1,netall2[,3:dim(netall2)[2]])
      write.table(netall,paste(experiment,"/Generated_Networks/",classification,"_",nameout,"_pvals_network.txt",sep=""),sep="\t",row.names = F,quote=F)
      ## Correct for multipletesting ##
      if(length(grep("Pval",colnames(netall)))>0){
        netall[,4]<-as.numeric(netall[,4])*length(netall[,1])
        netall[,6]<-as.numeric(netall[,6])*length(netall[,1])
        netall[,8]<-as.numeric(netall[,8])*length(netall[,1])
        netall[,10]<-as.numeric(netall[,10])*length(netall[,1])
        netall[,12]<-as.numeric(netall[,12])*length(netall[,1])
        netall[,14]<-as.numeric(netall[,14])*length(netall[,1])
        netall[,17]<-as.numeric(netall[,17])*length(netall[,1])
        netall[,19]<-as.numeric(netall[,19])*length(netall[,1])
        netall[,21]<-as.numeric(netall[,21])*length(netall[,1])
        colnames(netall)<-gsub("Pval","FDR",colnames(netall))
        netall[which(as.numeric(netall[,4])>1),4]<-1
        netall[which(as.numeric(netall[,6])>1),6]<-1
        netall[which(as.numeric(netall[,8])>1),8]<-1
        netall[which(as.numeric(netall[,10])>1),10]<-1
        netall[which(as.numeric(netall[,12])>1),12]<-1
        netall[which(as.numeric(netall[,14])>1),14]<-1
        netall[which(as.numeric(netall[,17])>1),17]<-1
        netall[which(as.numeric(netall[,19])>1),19]<-1
        netall[which(as.numeric(netall[,21])>1),21]<-1
      }
      ## Write the tables ##
      write.table(netall,paste(experiment,"/Generated_Networks/",classification,"_",nameout,"_fdr_network.txt",sep=""),quote=F,sep="\t",row.names=F)
      ## Convert the correlations and FDRs into 1s, 0s, and -1s ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
      palldir<-rep(0,length(netall[,1]))
      palldir[intersect(which(as.numeric(netall[,4])<=0.05),which(as.numeric(netall[,3])<0))]<-(-1)
      palldir[intersect(which(as.numeric(netall[,4])<=0.05),which(as.numeric(netall[,3])>0))]<-1
      salldir<-rep(0,length(netall[,1]))
      salldir[intersect(which(as.numeric(netall[,6])<=0.05),which(as.numeric(netall[,5])<0))]<-(-1)
      salldir[intersect(which(as.numeric(netall[,6])<=0.05),which(as.numeric(netall[,5])>0))]<-1
      pintdir<-rep(0,length(netall[,1]))
      pintdir[intersect(which(as.numeric(netall[,8])<=0.05),which(as.numeric(netall[,7])<0))]<-(-1)
      pintdir[intersect(which(as.numeric(netall[,8])<=0.05),which(as.numeric(netall[,7])>0))]<-1
      sintdir<-rep(0,length(netall[,1]))
      sintdir[intersect(which(as.numeric(netall[,10])<=0.05),which(as.numeric(netall[,9])<0))]<-(-1)
      sintdir[intersect(which(as.numeric(netall[,10])<=0.05),which(as.numeric(netall[,9])>0))]<-1
      punidir<-rep(0,length(netall[,1]))
      punidir[intersect(which(as.numeric(netall[,12])<=0.05),which(as.numeric(netall[,11])<0))]<-(-1)
      punidir[intersect(which(as.numeric(netall[,12])<=0.05),which(as.numeric(netall[,11])>0))]<-1
      sunidir<-rep(0,length(netall[,1]))
      sunidir[intersect(which(as.numeric(netall[,14])<=0.05),which(as.numeric(netall[,13])<0))]<-(-1)
      sunidir[intersect(which(as.numeric(netall[,14])<=0.05),which(as.numeric(netall[,13])>0))]<-1
      csalldir<-rep(0,length(netall[,1]))
      csalldir[intersect(which(as.numeric(netall[,17])<=0.05),which(as.numeric(netall[,16])<0))]<-(-1)
      csalldir[intersect(which(as.numeric(netall[,17])<=0.05),which(as.numeric(netall[,16])>0))]<-1
      csintdir<-rep(0,length(netall[,1]))
      csintdir[intersect(which(as.numeric(netall[,19])<=0.05),which(as.numeric(netall[,18])<0))]<-(-1)
      csintdir[intersect(which(as.numeric(netall[,19])<=0.05),which(as.numeric(netall[,18])>0))]<-1
      csunidir<-rep(0,length(netall[,1]))
      csunidir[intersect(which(as.numeric(netall[,21])<=0.05),which(as.numeric(netall[,20])<0))]<-(-1)
      csunidir[intersect(which(as.numeric(netall[,21])<=0.05),which(as.numeric(netall[,20])>0))]<-1
      
      
      binary<-cbind(netall[,1:2],palldir,salldir,pintdir,sintdir,punidir,sunidir,netall[,15],csalldir,csintdir,csunidir)
      ## Write both tables ##
      colnames(binary)<-c("Disease1","Disease2","PearsonAll","SpearmanAll","PearsonInt","SpearmanInt","PearsonUni","SpearmanUni","Fisher",
                          "CosineSimilarityAll","CosineSimilarityInt","CosineSimilarityUni")
      write.table(binary,paste(experiment,"/Generated_Networks/",classification,"_",nameout,"_binarized_network.txt",sep=""),quote=F,sep="\t",row.names=F)
      ## Which interactions are different?
      resultstab<-rbind(cbind(c(length(which(binary[,3]!=0))),
                              c(length(which(binary[,4]!=0))),
                              c(length(which(binary[,5]!=0))),
                              c(length(which(binary[,6]!=0))),
                              c(length(which(binary[,7]!=0))),
                              c(length(which(binary[,8]!=0))),
                              c(length(which(binary[,9]!=0))),
                              c(length(which(binary[,10]!=0))),
                              c(length(which(binary[,11]!=0))),
                              c(length(which(binary[,12]!=0)))))
      resultstab<-cbind(resultstab,length(binary[,1]),length(unique(c(binary$Disease1,binary$Disease2))))
      colnames(resultstab)<-c("PearsonAll","SpearmanAll","PearsonInt","SpearmanInt","PearsonUni","SpearmanUni","Fisher",
                              "CosineSimilarityAll","CosineSimilarityInt","CosineSimilarityUni","PotentialInteractions","Diseases")
      resultstab<-rbind(resultstab,round(100*(resultstab[1,]/resultstab[1,11]),2))
      rownames(resultstab)<-c("All","Percentage_All")
      write.table(resultstab,paste(experiment,"/Generated_Networks/",classification,"_",nameout,"_number_interactions.txt",sep=""),quote=F,sep="\t")
      resultados<-list("resultstab"=resultstab,"binarynet"=binary,"fdrnet"=netall)
    }
    ## If both tables have the same number of interactions and the interactions in the same order:
    if((dim(netall1)[1]==dim(netall2)[1] && length(which(netall1[,1]==netall2[,1]))==dim(netall1)[1] && length(which(netall1[,2]==netall2[,2]))==dim(netall2)[1])==FALSE){
      print("Error! Tables have different number of interactions or are ordered in a different way")
    }
  }
  print("Finished!")
  return(resultados)
}

if(as.numeric(args[1])==3){
  microicd9women<-correctformultipletesting("Women","ICD9","Microarrays")
  microicd9men<-correctformultipletesting("Men","ICD9","Microarrays")
  microicd9cases<-correctformultipletesting("Cases","ICD9","Microarrays")
  microicd9controls<-correctformultipletesting("Controls","ICD9","Microarrays")
  microicd9adjusted<-correctformultipletesting("Adjusted","ICD9","Microarrays")
  microicd10women<-correctformultipletesting("Women","ICD10","Microarrays")
  microicd10men<-correctformultipletesting("Men","ICD10","Microarrays")
  microicd10cases<-correctformultipletesting("Cases","ICD10","Microarrays")
  microicd10controls<-correctformultipletesting("Controls","ICD10","Microarrays")
  microicd10adjusted<-correctformultipletesting("Adjusted","ICD10","Microarrays")
}

#### Add the estimate of euclidean distances ####
## @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @
# nameout<-"Adjusted"
# classification<-"ICD9"
# experiment<-"Microarrays"
# themetric<-"euclidean"
# theoutput<-"EuclideanDistance"
set.seed(1)
build_euclideandistance_network<-function(nameout,classification,experiment,themetric,theoutput){
  ## Start the analysis!! ##
  ## @@ @@ @@ @@ @@ @@ @@ ##
  files<-list.files(paste(experiment,"/",classification,"/Comparisons/",nameout,"/DifferentialExpressions/",sep=""))
  ## Run the loop ##
  sdegs<-list() ; vsdegs<-list() ; vnum<-c() ; vname<-c()
  allgenes<-list()
  for(a in 1:length(files)){
    # a<-1
    tabla<-read.csv2(paste(experiment,"/",classification,"/Comparisons/",nameout,"/DifferentialExpressions/",files[a],sep=""),stringsAsFactors = F,sep="\t")
    ups<-rownames(tabla)[intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))>0))]
    upslogs<-as.numeric(as.character(tabla$logFC[intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))>0))]))
    downs<-rownames(tabla)[intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))<0))]
    downslogs<-as.numeric(as.character(tabla$logFC[intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))<0))]))
    names(downslogs)<-downs
    names(upslogs)<-ups
    sdegs[[gsub(".txt","",files[a])]]$Up<-ups
    sdegs[[gsub(".txt","",files[a])]]$Down<-downs
    vsdegs[[gsub(".txt","",files[a])]]$Up<-upslogs
    vsdegs[[gsub(".txt","",files[a])]]$Down<-downslogs
    vnum<-rbind(vnum,c(length(ups),length(downs)))
    vname<-c(vname,gsub(".+_","",gsub(".txt","",files[a])))
    todosgenes<-tabla$logFC ; names(todosgenes)<-rownames(tabla)
    allgenes[[gsub(".txt","",files[a])]]<-todosgenes
    print(paste(round((a/length(files))*100,2),"%",sep=""))
  }
  print("Vectors generated")
  
  ## Similarities between pairs of diseases ##
  netall<-c()
  for(a in 1:(length(names(vsdegs))-1)){
    for(b in (a+1):length(names(vsdegs))){
      ## Over all the genes ##
      ## @@ @@ @@  @@ @@ @@ ##
      commongenes<-intersect(names(allgenes[[a]]),names(allgenes[[b]]))
      avec1<-as.numeric(as.character(allgenes[[a]][commongenes]))
      avec2<-as.numeric(as.character(allgenes[[b]][commongenes]))
      acossim<-dist(rbind(avec1,avec2),method = themetric)[[1]]
      rancossim<-numeric(length=10000)
      for(ran in 1:10000){
        rancossim[ran]<-dist(rbind(sample(avec1),avec2),method = themetric)[[1]]
      }
      adifference<-mean(rancossim)-acossim
      if(adifference<0){
        apval<-length(which(rancossim>acossim))/10000
      }
      if(adifference>0){
        apval<-length(which(rancossim<acossim))/10000
      }
      
      ## Over the intersection of the sDEGs ##
      ## @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ ##
      unof<-vsdegs[[a]] ; dosf<-vsdegs[[b]]
      allunof<-c(unof$Up,unof$Down) ; alldosf<-c(dosf$Up,dosf$Down)
      fintersection<-intersect(names(allunof),names(alldosf))
      if(length(fintersection)<=3){
        icossim<-NA
        ipval<-NA
        idifference<-NA
        irancossim<-NA
      }
      if(length(fintersection)>3){
        ivec1<-as.numeric(as.character(allunof[fintersection]))
        ivec2<-as.numeric(as.character(alldosf[fintersection]))
        icossim<-dist(rbind(ivec1,ivec2),method = themetric)[[1]]
        irancossim<-numeric(length=10000)
        for(ran in 1:10000){
          irancossim[ran]<-dist(rbind(sample(ivec1),ivec2),method=themetric)[[1]]
        }
        idifference<-mean(irancossim)-icossim
        if(idifference<0){
          ipval<-length(which(irancossim>icossim))/10000
        }
        if(idifference>0){
          ipval<-length(which(irancossim<icossim))/10000
        }
      }
      ## Over the union of sDEGs ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ ##
      funion<-unique(c(names(allunof),names(alldosf)))
      unioncommon<-intersect(names(which(is.na(allgenes[[a]][funion])==FALSE)),names(which(is.na(allgenes[[b]][funion])==FALSE)))
      numberifunion<-length(unioncommon)
      if(numberifunion<=3){
        ucossim<-NA
        upval<-NA
        udifference<-NA
        urancossim<-NA
      }
      if(numberifunion>3){
        uvec1<-as.numeric(as.character(allgenes[[a]][unioncommon]))
        uvec2<-as.numeric(as.character(allgenes[[b]][unioncommon]))
        ucossim<-dist(rbind(uvec1,uvec2),method = themetric)[[1]]
        urancossim<-numeric(length=10000)
        for(ran in 1:10000){
          urancossim[ran]<-dist(rbind(sample(uvec1),uvec2),method=themetric)[[1]]
        }
        udifference<-mean(urancossim)-ucossim
        if(udifference<0){
          upval<-length(which(urancossim>ucossim))/10000
        }
        if(udifference>0){
          upval<-length(which(urancossim<ucossim))/10000
        }
      }
      ## Put all the correlations together ##
      ## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
      netall<-rbind(netall,c(names(allgenes)[c(a,b)],acossim,mean(rancossim),adifference,apval,icossim,mean(irancossim),idifference,ipval,ucossim,mean(urancossim),udifference,upval))
    }
    print(paste(round((a/length(names(allgenes)))*100,2),"%",sep=""))
  }
  ## Add column names ##
  colnames(netall)<-c("Disease1","Disease2",paste("All",theoutput,sep=""),paste("AllRandom",theoutput,sep=""),
                      paste("AllDifference",theoutput,sep=""),paste("AllPval",theoutput,sep=""),
                      paste("Int",theoutput,sep=""),paste("IntRandom",theoutput,sep=""),
                      paste("IntDifference",theoutput,sep=""),paste("IntPval",theoutput,sep=""),
                      paste("Uni",theoutput,sep=""),paste("UniRandom",theoutput,sep=""),
                      paste("UniDifference",theoutput,sep=""),paste("UniPval",theoutput,sep=""))
  ## Write the tables ##
  write.table(netall,paste(experiment,"/OtherMetrics/",theoutput,"/",classification,"_",nameout,"_pvals_network.txt",sep=""),quote=F,sep="\t",row.names=F)
  return(netall)
}

## Generate the myrun file ##
# myrun<-c(paste("Rscript 05_3_CalculateOtherSimilarities_Genes.R 2",
#              c(paste("Adjusted","ICD9","Microarrays"),paste("Adjusted","ICD10","Microarrays"),
#                paste("Women","ICD9","Microarrays"),paste("Men","ICD9","Microarrays"),paste("Women","ICD10","Microarrays"),paste("Men","ICD10","Microarrays"),
#                paste("Cases","ICD9","Microarrays"),paste("Controls","ICD9","Microarrays"),paste("Cases","ICD10","Microarrays"),paste("Controls","ICD10","Microarrays"))),
#          paste("Rscript 05_3_CalculateOtherSimilarities_Genes.R 4",
#                c(paste("Adjusted","ICD9","Microarrays"),paste("Adjusted","ICD10","Microarrays"),
#                  paste("Women","ICD9","Microarrays"),paste("Men","ICD9","Microarrays"),paste("Women","ICD10","Microarrays"),paste("Men","ICD10","Microarrays"),
#                  paste("Cases","ICD9","Microarrays"),paste("Controls","ICD9","Microarrays"),paste("Cases","ICD10","Microarrays"),paste("Controls","ICD10","Microarrays"))),
#          paste("Rscript 05_3_CalculateOtherSimilarities_Genes.R 5",
#                c(paste("Adjusted","ICD9","Microarrays"),paste("Adjusted","ICD10","Microarrays"),
#                  paste("Women","ICD9","Microarrays"),paste("Men","ICD9","Microarrays"),paste("Women","ICD10","Microarrays"),paste("Men","ICD10","Microarrays"),
#                  paste("Cases","ICD9","Microarrays"),paste("Controls","ICD9","Microarrays"),paste("Cases","ICD10","Microarrays"),paste("Controls","ICD10","Microarrays"))),
#          paste("Rscript 05_3_CalculateOtherSimilarities_Genes.R 6",
#                c(paste("Adjusted","ICD9","Microarrays"),paste("Adjusted","ICD10","Microarrays"),
#                  paste("Women","ICD9","Microarrays"),paste("Men","ICD9","Microarrays"),paste("Women","ICD10","Microarrays"),paste("Men","ICD10","Microarrays"),
#                  paste("Cases","ICD9","Microarrays"),paste("Controls","ICD9","Microarrays"),paste("Cases","ICD10","Microarrays"),paste("Controls","ICD10","Microarrays"))))
# write.table(myrun,"my_run_05_3_CalculateOtherSimilarities_Genes.txt",quote=F,sep="\t",row.names=F,col.names=F)



## Run in Mare ##
## Euclidean ##
if(as.numeric(args[1])==4){
  nameout<-args[2]
  classification<-args[3]
  experiment<-args[4]
  print(nameout)
  print(classification)
  print(experiment)
  salida<-build_euclideandistance_network(nameout,classification,experiment,"euclidean","EuclideanDistance")
}

## Manhattan ##
if(as.numeric(args[1])==5){
  nameout<-args[2]
  classification<-args[3]
  experiment<-args[4]
  print(nameout)
  print(classification)
  print(experiment)
  salida<-build_euclideandistance_network(nameout,classification,experiment,"manhattan","ManhattanDistance")
}

## Canberra ##
if(as.numeric(args[1])==6){
  nameout<-args[2]
  classification<-args[3]
  experiment<-args[4]
  print(nameout)
  print(classification)
  print(experiment)
  salida<-build_euclideandistance_network(nameout,classification,experiment,"canberra","CanberraDistance")
}

#### Create a function to correct for multiple testing ####
# nameout<-"Women"
# classification<-"ICD9"
# experiment<-"Microarrays"
correctformultipletesting2<-function(nameout,classification,experiment){
  netall1<-read.csv2(paste(experiment,"/Generated_Networks/",classification,"_",nameout,"_pvals_network.txt",sep=""),stringsAsFactors = F,sep="\t")
  if(dim(netall1)[2]<=21){
    netall2<-read.csv2(paste(experiment,"/OtherMetrics/EuclideanDistance/",classification,"_",nameout,"_pvals_network.txt",sep=""),stringsAsFactors = F,sep="\t")
    netall3<-read.csv2(paste(experiment,"/OtherMetrics/ManhattanDistance/",classification,"_",nameout,"_pvals_network.txt",sep=""),stringsAsFactors = F,sep="\t")
    netall4<-read.csv2(paste(experiment,"/OtherMetrics/CanberraDistance/",classification,"_",nameout,"_pvals_network.txt",sep=""),stringsAsFactors = F,sep="\t")
    ## If both tables have the same number of interactions and the interactions in the same order:
    if(dim(netall1)[1]==dim(netall2)[1] && dim(netall1)[1]==dim(netall3)[1] && dim(netall1)[1]==dim(netall4)[1] && 
       length(which(netall1[,1]==netall2[,1]))==dim(netall1)[1] && length(which(netall1[,2]==netall2[,2]))==dim(netall2)[1] &&
       length(which(netall1[,1]==netall3[,1]))==dim(netall1)[1] && length(which(netall1[,2]==netall3[,2]))==dim(netall3)[1] &&
       length(which(netall1[,1]==netall4[,1]))==dim(netall1)[1] && length(which(netall1[,2]==netall4[,2]))==dim(netall4)[1]){
      netall<-cbind(netall1,netall2[,c(5,6,9,10,13,14)],netall3[,c(5,6,9,10,13,14)],netall4[,c(5,6,9,10,13,14)])
      write.table(netall,paste(experiment,"/Generated_Networks/",classification,"_",nameout,"_pvals_network.txt",sep=""),sep="\t",row.names = F,quote=F)
      ## Correct for multipletesting ##
      if(length(grep("Pval",colnames(netall)))>0){
        secuencia<-grep("Pval",colnames(netall))
        for(a in secuencia){
          netall[,a]<-as.numeric(netall[,a])*length(netall[,1])
          if(length(which(as.numeric(netall[,a])>1))>0){
            netall[which(as.numeric(netall[,a])>1),a]<-1
          }
        }
        colnames(netall)<-gsub("Pval","FDR",colnames(netall))
      }
      ## Write the tables ##
      write.table(netall,paste(experiment,"/Generated_Networks/",classification,"_",nameout,"_fdr_network.txt",sep=""),quote=F,sep="\t",row.names=F)
      ## Convert the correlations and FDRs into 1s, 0s, and -1s ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
      palldir<-rep(0,length(netall[,1]))
      palldir[intersect(which(as.numeric(netall[,4])<=0.05),which(as.numeric(netall[,3])<0))]<-(-1)
      palldir[intersect(which(as.numeric(netall[,4])<=0.05),which(as.numeric(netall[,3])>0))]<-1
      salldir<-rep(0,length(netall[,1]))
      salldir[intersect(which(as.numeric(netall[,6])<=0.05),which(as.numeric(netall[,5])<0))]<-(-1)
      salldir[intersect(which(as.numeric(netall[,6])<=0.05),which(as.numeric(netall[,5])>0))]<-1
      pintdir<-rep(0,length(netall[,1]))
      pintdir[intersect(which(as.numeric(netall[,8])<=0.05),which(as.numeric(netall[,7])<0))]<-(-1)
      pintdir[intersect(which(as.numeric(netall[,8])<=0.05),which(as.numeric(netall[,7])>0))]<-1
      sintdir<-rep(0,length(netall[,1]))
      sintdir[intersect(which(as.numeric(netall[,10])<=0.05),which(as.numeric(netall[,9])<0))]<-(-1)
      sintdir[intersect(which(as.numeric(netall[,10])<=0.05),which(as.numeric(netall[,9])>0))]<-1
      punidir<-rep(0,length(netall[,1]))
      punidir[intersect(which(as.numeric(netall[,12])<=0.05),which(as.numeric(netall[,11])<0))]<-(-1)
      punidir[intersect(which(as.numeric(netall[,12])<=0.05),which(as.numeric(netall[,11])>0))]<-1
      sunidir<-rep(0,length(netall[,1]))
      sunidir[intersect(which(as.numeric(netall[,14])<=0.05),which(as.numeric(netall[,13])<0))]<-(-1)
      sunidir[intersect(which(as.numeric(netall[,14])<=0.05),which(as.numeric(netall[,13])>0))]<-1
      csalldir<-rep(0,length(netall[,1]))
      csalldir[intersect(which(as.numeric(netall[,17])<=0.05),which(as.numeric(netall[,16])<0))]<-(-1)
      csalldir[intersect(which(as.numeric(netall[,17])<=0.05),which(as.numeric(netall[,16])>0))]<-1
      csintdir<-rep(0,length(netall[,1]))
      csintdir[intersect(which(as.numeric(netall[,19])<=0.05),which(as.numeric(netall[,18])<0))]<-(-1)
      csintdir[intersect(which(as.numeric(netall[,19])<=0.05),which(as.numeric(netall[,18])>0))]<-1
      csunidir<-rep(0,length(netall[,1]))
      csunidir[intersect(which(as.numeric(netall[,21])<=0.05),which(as.numeric(netall[,20])<0))]<-(-1)
      csunidir[intersect(which(as.numeric(netall[,21])<=0.05),which(as.numeric(netall[,20])>0))]<-1
      ## Euclidean ##
      eualldir<-rep(0,length(netall[,1]))
      eualldir[intersect(which(as.numeric(netall[,23])<=0.05),which(as.numeric(netall[,22])<0))]<-(-1)
      eualldir[intersect(which(as.numeric(netall[,23])<=0.05),which(as.numeric(netall[,22])>0))]<-1
      euintdir<-rep(0,length(netall[,1]))
      euintdir[intersect(which(as.numeric(netall[,25])<=0.05),which(as.numeric(netall[,24])<0))]<-(-1)
      euintdir[intersect(which(as.numeric(netall[,25])<=0.05),which(as.numeric(netall[,24])>0))]<-1
      euunidir<-rep(0,length(netall[,1]))
      euunidir[intersect(which(as.numeric(netall[,27])<=0.05),which(as.numeric(netall[,26])<0))]<-(-1)
      euunidir[intersect(which(as.numeric(netall[,27])<=0.05),which(as.numeric(netall[,26])>0))]<-1
      ## Manhattan ##
      maalldir<-rep(0,length(netall[,1]))
      maalldir[intersect(which(as.numeric(netall[,29])<=0.05),which(as.numeric(netall[,28])<0))]<-(-1)
      maalldir[intersect(which(as.numeric(netall[,29])<=0.05),which(as.numeric(netall[,28])>0))]<-1
      maintdir<-rep(0,length(netall[,1]))
      maintdir[intersect(which(as.numeric(netall[,31])<=0.05),which(as.numeric(netall[,30])<0))]<-(-1)
      maintdir[intersect(which(as.numeric(netall[,31])<=0.05),which(as.numeric(netall[,30])>0))]<-1
      maunidir<-rep(0,length(netall[,1]))
      maunidir[intersect(which(as.numeric(netall[,33])<=0.05),which(as.numeric(netall[,32])<0))]<-(-1)
      maunidir[intersect(which(as.numeric(netall[,33])<=0.05),which(as.numeric(netall[,32])>0))]<-1
      ## Canberra ##
      caalldir<-rep(0,length(netall[,1]))
      caalldir[intersect(which(as.numeric(netall[,35])<=0.05),which(as.numeric(netall[,34])<0))]<-(-1)
      caalldir[intersect(which(as.numeric(netall[,35])<=0.05),which(as.numeric(netall[,34])>0))]<-1
      caintdir<-rep(0,length(netall[,1]))
      caintdir[intersect(which(as.numeric(netall[,37])<=0.05),which(as.numeric(netall[,36])<0))]<-(-1)
      caintdir[intersect(which(as.numeric(netall[,37])<=0.05),which(as.numeric(netall[,36])>0))]<-1
      caunidir<-rep(0,length(netall[,1]))
      caunidir[intersect(which(as.numeric(netall[,39])<=0.05),which(as.numeric(netall[,38])<0))]<-(-1)
      caunidir[intersect(which(as.numeric(netall[,39])<=0.05),which(as.numeric(netall[,38])>0))]<-1
      
      
      
      binary<-cbind(netall[,1:2],palldir,salldir,pintdir,sintdir,punidir,sunidir,netall[,15],csalldir,csintdir,csunidir,
                    eualldir,euintdir,euunidir,maalldir,maintdir,maunidir,caalldir,caintdir,caunidir)
      ## Write both tables ##
      colnames(binary)<-c("Disease1","Disease2","PearsonAll","SpearmanAll","PearsonInt","SpearmanInt","PearsonUni","SpearmanUni","Fisher",
                          "CosineSimilarityAll","CosineSimilarityInt","CosineSimilarityUni","EuclideanDistanceAll","EuclideanDistanceInt","EuclideanDistanceUni",
                          "ManhattanDistanceAll","ManhattanDistanceInt","ManhattanDistanceUni","CanberraDistanceAll","CanberraDistanceInt","CanberraDistanceUni")
      write.table(binary,paste(experiment,"/Generated_Networks/",classification,"_",nameout,"_binarized_network.txt",sep=""),quote=F,sep="\t",row.names=F)
      ## Which interactions are different?
      resultstab<-rbind(cbind(c(length(which(binary[,3]!=0))),
                              c(length(which(binary[,4]!=0))),
                              c(length(which(binary[,5]!=0))),
                              c(length(which(binary[,6]!=0))),
                              c(length(which(binary[,7]!=0))),
                              c(length(which(binary[,8]!=0))),
                              c(length(which(binary[,9]!=0))),
                              c(length(which(binary[,10]!=0))),
                              c(length(which(binary[,11]!=0))),
                              c(length(which(binary[,12]!=0))),
                              c(length(which(binary[,13]!=0))),
                              c(length(which(binary[,14]!=0))),
                              c(length(which(binary[,15]!=0))),
                              c(length(which(binary[,16]!=0))),
                              c(length(which(binary[,17]!=0))),
                              c(length(which(binary[,18]!=0))),
                              c(length(which(binary[,19]!=0))),
                              c(length(which(binary[,20]!=0))),
                              c(length(which(binary[,21]!=0)))))
      resultstab<-cbind(resultstab,length(binary[,1]),length(unique(c(binary$Disease1,binary$Disease2))))
      colnames(resultstab)<-c("PearsonAll","SpearmanAll","PearsonInt","SpearmanInt","PearsonUni","SpearmanUni","Fisher",
                              "CosineSimilarityAll","CosineSimilarityInt","CosineSimilarityUni",
                              "EuclideanDistanceAll","EuclideanDistanceInt","EuclideanDistanceUni",
                              "ManhattanDistanceAll","ManhattanDistanceInt","ManhattanDistanceUni",
                              "CanberraDistanceAll","CanberraDistanceInt","CanberraDistanceUni",
                              "PotentialInteractions","Diseases")
      resultstab<-rbind(resultstab,round(100*(resultstab[1,]/resultstab[1,14]),2))
      rownames(resultstab)<-c("All","Percentage_All")
      write.table(resultstab,paste(experiment,"/Generated_Networks/",classification,"_",nameout,"_number_interactions.txt",sep=""),quote=F,sep="\t")
      resultados<-list("resultstab"=resultstab,"binarynet"=binary,"fdrnet"=netall)
    }
    ## If both tables have the same number of interactions and the interactions in the same order:
    if((dim(netall1)[1]==dim(netall2)[1] && length(which(netall1[,1]==netall2[,1]))==dim(netall1)[1] && length(which(netall1[,2]==netall2[,2]))==dim(netall2)[1])==FALSE){
      print("Error! Tables have different number of interactions or are ordered in a different way")
    }
  }
  print("Finished!")
  return(resultados)
}

if(as.numeric(args[1])==7){
  microicd9women<-correctformultipletesting2("Women","ICD9","Microarrays")
  microicd9men<-correctformultipletesting2("Men","ICD9","Microarrays")
  microicd9cases<-correctformultipletesting2("Cases","ICD9","Microarrays")
  microicd9controls<-correctformultipletesting2("Controls","ICD9","Microarrays")
  microicd9adjusted<-correctformultipletesting2("Adjusted","ICD9","Microarrays")
  microicd10women<-correctformultipletesting2("Women","ICD10","Microarrays")
  microicd10men<-correctformultipletesting2("Men","ICD10","Microarrays")
  microicd10cases<-correctformultipletesting2("Cases","ICD10","Microarrays")
  microicd10controls<-correctformultipletesting2("Controls","ICD10","Microarrays")
  microicd10adjusted<-correctformultipletesting2("Adjusted","ICD10","Microarrays")
}













