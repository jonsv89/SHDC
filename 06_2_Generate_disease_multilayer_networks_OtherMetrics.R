## CODE FOR THE GENERATION OF THE DISEASE NETWORKS FOR MEN, WOMEN, AND ADJUSTING BY SEX BASED ON GENES SEPARATELY FOR EACH REACTOME CATEGORY ##
## Developed by Jon Sanchez-Valle & Beatriz Urda-Garcia
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es @ beatriz.urda@bsc.es

args = commandArgs(trailingOnly=TRUE)
library(lsa)

#### Add the estimate of cosine similarities ####
## @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ ##
## Identify the reactome parents and select their genes ##
reactomeparents<-read.csv2("Microarrays/Data/GSEAReactome_ParentNames.txt",stringsAsFactors = F,sep="\t")
reactomeparents<-reactomeparents[order(reactomeparents[,3]),]
factreactomeparents<-reactomeparents[,3] ; names(factreactomeparents)<-reactomeparents[,1]
## Get the genes associated to each pathway ##
gseareac<-read.csv2("GSEAfiles/gmt/c2.cp.reactome.v2023.1.Hs.symbols.gmt",stringsAsFactors = F)
pathwaygenelist<-list()
for(a in 1:length(gseareac[,1])){
  # a<-1
  pathgenevector<-strsplit(gseareac[a,1],"\t")[[1]]
  pathwaygenelist[[pathgenevector[1]]]<-pathgenevector[3:length(pathgenevector)]
}
## Get the genes associated to each reactome parent ##
theparents<-unique(reactomeparents$ParentName)
if(length(grep("Digestion and absorption",theparents))>0){theparents<-theparents[-grep("Digestion and absorption",theparents)]}
reactomeparentgenelist<-list()
genesperparent<-c()
for(a in 1:length(theparents)){
  # a<-1
  thepathways<-reactomeparents$GSEA[which(reactomeparents$ParentName==theparents[a])]
  thegenes<-c();for(t in 1:length(thepathways)){thegenes<-c(thegenes,pathwaygenelist[[thepathways[t]]])}
  reactomeparentgenelist[[theparents[a]]]<-unique(thegenes)
  genesperparent<-rbind(genesperparent,c(theparents[a],length(unique(thegenes))))
}
genesperparent<-genesperparent[order(as.numeric(genesperparent[,2]),decreasing = T),]

## Add all the genes associated with the mitochondria ##
if(length(grep("Mitochondria All",names(reactomeparentgenelist)))==0){
  mitogenes<-read.csv2("Microarrays/Data/MitoCarta.csv",stringsAsFactors = F)
  reactomeparentgenelist[["Mitochondria All"]]<-toupper(unique(mitogenes$Symbol))
}

if("OtherMetrics_Pathways"%in%list.files("Microarrays/")==FALSE){dir.create("Microarrays/OtherMetrics_Pathways")}

## Remove the "()" to avoid problems while running everything ##
names(reactomeparentgenelist)<-gsub("\\)","",gsub("\\(","",names(reactomeparentgenelist)))

# nameout<-"Women"
# classification<-"ICD9"
# experiment<-"Microarrays"
# thepathway<-"Mitochondria_All"
# thepathway<-gsub(" ","_",names(reactomeparentgenelist))[16]
set.seed(1)
build_cosinesimilarity_network<-function(nameout,classification,experiment,thepathway){
  ## Start the analysis!! ##
  ## @@ @@ @@ @@ @@ @@ @@ ##
  thegene<-reactomeparentgenelist[[which(names(reactomeparentgenelist)==gsub("_"," ",thepathway))]]
  files<-list.files(paste(experiment,"/",classification,"/Comparisons/",nameout,"/DifferentialExpressions/",sep=""))
  ## Run the loop ##
  sdegs<-list() ; vsdegs<-list() ; vnum<-c() ; vname<-c()
  allgenes<-list()
  for(a in 1:length(files)){
    # a<-1
    tabla<-read.csv2(paste(experiment,"/",classification,"/Comparisons/",nameout,"/DifferentialExpressions/",files[a],sep=""),stringsAsFactors = F,sep="\t")
    ## Now we select only the genes that belong to the pathway ##
    tabla<-tabla[intersect(thegene,rownames(tabla)),]
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
      ## If there is at least one gene in common keep going ##
      if(length(commongenes)>1){
        avec1<-as.numeric(as.character(allgenes[[a]][commongenes]))
        avec2<-as.numeric(as.character(allgenes[[b]][commongenes]))
        ## Distances ##
        acossim<-cosine(avec1,avec2)[1,1]
        aeucli<-dist(rbind(avec1,avec2),method = "euclidean")[[1]]
        amanha<-dist(rbind(avec1,avec2),method = "manhattan")[[1]]
        acanbe<-dist(rbind(avec1,avec2),method = "canberra")[[1]]
        ## Vectors for random values ##
        arancossim<-numeric(length=10000)
        araneucli<-numeric(length=10000)
        aranmanha<-numeric(length=10000)
        arancanbe<-numeric(length=10000)
        ## Calculate random distances ##
        for(ran in 1:10000){
          arancossim[ran]<-cosine(sample(avec1),avec2)[1,1]
          araneucli[ran]<-dist(rbind(sample(avec1),avec2),method = "euclidean")[[1]]
          aranmanha[ran]<-dist(rbind(sample(avec1),avec2),method = "manhattan")[[1]]
          arancanbe[ran]<-dist(rbind(sample(avec1),avec2),method = "canberra")[[1]]
        }
        ## pval cosine similarity ##
        if(acossim<0){apval<-length(which(arancossim<acossim))/10000}
        if(acossim>0){apval<-length(which(arancossim>acossim))/10000}
        ## euclidean - difference and pval ##
        aeuclidifference<-mean(araneucli)-aeucli
        if(aeuclidifference<0){aeuclipval<-length(which(araneucli>aeucli))/10000}
        if(aeuclidifference>0){aeuclipval<-length(which(araneucli<aeucli))/10000}
        ## manhattan - difference and pval ##
        amanhadifference<-mean(aranmanha)-amanha
        if(amanhadifference<0){amanhapval<-length(which(aranmanha>amanha))/10000}
        if(amanhadifference>0){amanhapval<-length(which(aranmanha<amanha))/10000}
        ## canberra - difference and pval ##
        acanbedifference<-mean(arancanbe)-acanbe
        if(acanbedifference<0){acanbepval<-length(which(arancanbe>acanbe))/10000}
        if(acanbedifference>0){acanbepval<-length(which(arancanbe<acanbe))/10000}
        
        ## Over the intersection of the sDEGs ##
        ## @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ ##
        unof<-vsdegs[[a]] ; dosf<-vsdegs[[b]]
        allunof<-c(unof$Up,unof$Down) ; alldosf<-c(dosf$Up,dosf$Down)
        fintersection<-intersect(names(allunof),names(alldosf))
        if(length(fintersection)<=3){
          icossim<-NA ; ipval<-NA
          ieuclidifference<-NA ; ieuclipval<-NA ; iraneucli<-NA ; ieucli<-NA
          imanhadifference<-NA ; imanhapval<-NA ; iranmanha<-NA ; imanha<-NA
          icanbedifference<-NA ; icanbepval<-NA ; irancanbe<-NA ; icanbe<-NA
        }
        if(length(fintersection)>3){
          ivec1<-as.numeric(as.character(allunof[fintersection]))
          ivec2<-as.numeric(as.character(alldosf[fintersection]))
          ## Distances ##
          icossim<-cosine(ivec1,ivec2)[1,1]
          ieucli<-dist(rbind(ivec1,ivec2),method = "euclidean")[[1]]
          imanha<-dist(rbind(ivec1,ivec2),method = "manhattan")[[1]]
          icanbe<-dist(rbind(ivec1,ivec2),method = "canberra")[[1]]
          ## Vectors for random values ##
          irancossim<-numeric(length=10000)
          iraneucli<-numeric(length=10000)
          iranmanha<-numeric(length=10000)
          irancanbe<-numeric(length=10000)
          ## Calculate random distances ##
          for(ran in 1:10000){
            irancossim[ran]<-cosine(sample(ivec1),ivec2)[1,1]
            iraneucli[ran]<-dist(rbind(sample(ivec1),ivec2),method = "euclidean")[[1]]
            iranmanha[ran]<-dist(rbind(sample(ivec1),ivec2),method = "manhattan")[[1]]
            irancanbe[ran]<-dist(rbind(sample(ivec1),ivec2),method = "canberra")[[1]]
          }
          ## pval cosine similarity ##
          if(icossim<0){ipval<-length(which(irancossim<icossim))/10000}
          if(icossim>0){ipval<-length(which(irancossim>icossim))/10000}
          ## euclidean - difference and pval ##
          ieuclidifference<-mean(iraneucli)-ieucli
          if(ieuclidifference<0){ieuclipval<-length(which(iraneucli>ieucli))/10000}
          if(ieuclidifference>0){ieuclipval<-length(which(iraneucli<ieucli))/10000}
          ## manhattan - difference and pval ##
          imanhadifference<-mean(iranmanha)-imanha
          if(imanhadifference<0){imanhapval<-length(which(iranmanha>imanha))/10000}
          if(imanhadifference>0){imanhapval<-length(which(iranmanha<imanha))/10000}
          ## canberra - difference and pval ##
          icanbedifference<-mean(irancanbe)-icanbe
          if(icanbedifference<0){icanbepval<-length(which(irancanbe>icanbe))/10000}
          if(icanbedifference>0){icanbepval<-length(which(irancanbe<icanbe))/10000}
        }
        ## Over the union of sDEGs ##
        ## @@ @@ @@ @@ @@ @@ @@ @@ ##
        funion<-unique(c(names(allunof),names(alldosf)))
        unioncommon<-intersect(names(which(is.na(allgenes[[a]][funion])==FALSE)),names(which(is.na(allgenes[[b]][funion])==FALSE)))
        numberifunion<-length(unioncommon)
        if(numberifunion<=3){
          ucossim<-NA ; upval<-NA
          ueuclidifference<-NA ; ueuclipval<-NA ; uraneucli<-NA ; ueucli<-NA
          umanhadifference<-NA ; umanhapval<-NA ; uranmanha<-NA ; umanha<-NA
          ucanbedifference<-NA ; ucanbepval<-NA ; urancanbe<-NA ; ucanbe<-NA
        }
        if(numberifunion>3){
          uvec1<-as.numeric(as.character(allgenes[[a]][unioncommon]))
          uvec2<-as.numeric(as.character(allgenes[[b]][unioncommon]))
          ## Distances ##
          ucossim<-cosine(uvec1,uvec2)[1,1]
          ueucli<-dist(rbind(uvec1,uvec2),method = "euclidean")[[1]]
          umanha<-dist(rbind(uvec1,uvec2),method = "manhattan")[[1]]
          ucanbe<-dist(rbind(uvec1,uvec2),method = "canberra")[[1]]
          ## Vectors for random values ##
          urancossim<-numeric(length=10000)
          uraneucli<-numeric(length=10000)
          uranmanha<-numeric(length=10000)
          urancanbe<-numeric(length=10000)
          ## Calculate random distances ##
          for(ran in 1:10000){
            urancossim[ran]<-cosine(sample(uvec1),uvec2)[1,1]
            uraneucli[ran]<-dist(rbind(sample(uvec1),uvec2),method = "euclidean")[[1]]
            uranmanha[ran]<-dist(rbind(sample(uvec1),uvec2),method = "manhattan")[[1]]
            urancanbe[ran]<-dist(rbind(sample(uvec1),uvec2),method = "canberra")[[1]]
          }
          ## pval cosine similarity ##
          if(ucossim<0){upval<-length(which(urancossim<ucossim))/10000}
          if(ucossim>0){upval<-length(which(urancossim>ucossim))/10000}
          ## euclidean - difference and pval ##
          ueuclidifference<-mean(uraneucli)-ueucli
          if(ueuclidifference<0){ueuclipval<-length(which(uraneucli>ueucli))/10000}
          if(ueuclidifference>0){ueuclipval<-length(which(uraneucli<ueucli))/10000}
          ## manhattan - difference and pval ##
          umanhadifference<-mean(uranmanha)-umanha
          if(umanhadifference<0){umanhapval<-length(which(uranmanha>umanha))/10000}
          if(umanhadifference>0){umanhapval<-length(which(uranmanha<umanha))/10000}
          ## canberra - difference and pval ##
          ucanbedifference<-mean(urancanbe)-ucanbe
          if(ucanbedifference<0){ucanbepval<-length(which(urancanbe>ucanbe))/10000}
          if(ucanbedifference>0){ucanbepval<-length(which(urancanbe<ucanbe))/10000}
        }
        ## Put all the correlations together ##
        ## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
        netall<-rbind(netall,c(names(allgenes)[c(a,b)],acossim,apval,icossim,ipval,ucossim,upval,
                               aeucli,mean(araneucli),aeuclidifference,aeuclipval,ieucli,mean(iraneucli),ieuclidifference,ieuclipval,ueucli,mean(uraneucli),ueuclidifference,ueuclipval,
                               amanha,mean(aranmanha),amanhadifference,amanhapval,imanha,mean(iranmanha),imanhadifference,imanhapval,umanha,mean(uranmanha),umanhadifference,umanhapval,
                               acanbe,mean(arancanbe),acanbedifference,acanbepval,icanbe,mean(irancanbe),icanbedifference,icanbepval,ucanbe,mean(urancanbe),ucanbedifference,ucanbepval))
      }
    }
    print(paste(round((a/length(names(allgenes)))*100,2),"%",sep=""))
  }
  ## Add column names ##
  colnames(netall)<-c("Disease1","Disease2","AllCosineSimilarity","AllPvalCosineSimilarity",
                      "IntCosineSimilarity","IntPvalCosineSimilarity",
                      "UniCosineSimilarity","UniPvalCosineSimilarity",
                      "AllEuclidean","AllRandomEuclidean",
                      "AllDifferenceEuclidean","AllPvalEuclidean",
                      "IntEuclidean","IntRandomEuclidean",
                      "IntDifferenceEuclidean","IntPvalEuclidean",
                      "UniEuclidean","UniRandomEuclidean",
                      "UniDifferenceEuclidean","UniPvalEuclidean",
                      "AllManhattan","AllRandomManhattan",
                      "AllDifferenceManhattan","AllPvalManhattan",
                      "IntManhattan","IntRandomManhattan",
                      "IntDifferenceManhattan","IntPvalManhattan",
                      "UniManhattan","UniRandomManhattan",
                      "UniDifferenceManhattan","UniPvalManhattan",
                      "AllCanberra","AllRandomCanberra",
                      "AllDifferenceCanberra","AllPvalCanberra",
                      "IntCanberra","IntRandomCanberra",
                      "IntDifferenceCanberra","IntPvalCanberra",
                      "UniCanberra","UniRandomCanberra",
                      "UniDifferenceCanberra","UniPvalCanberra")
  ## Write the tables ##
  write.table(netall,paste(experiment,"/OtherMetrics_Pathways/",thepathway,"_",classification,"_",nameout,"_pvals_network.txt",sep=""),quote=F,sep="\t",row.names=F)
  return(netall)
}

## Run in Mare ##
if(as.numeric(args[1])==2){
  nameout<-args[2]
  classification<-args[3]
  experiment<-args[4]
  thepathway<-args[5]
  print(nameout)
  print(classification)
  print(experiment)
  print(thepathway)
  salida<-build_cosinesimilarity_network(nameout,classification,experiment,thepathway)
}

# nameout<-"Women"
# classification<-"ICD9"
# experiment<-"Microarrays"
# thepathway<-gsub(" ","_",names(reactomeparentgenelist)[8])
## Create a function to correct for multiple testing ##
correctformultipletesting<-function(nameout,classification,experiment,thepathway){
  netall1<-read.csv2(paste(experiment,"/Generated_Networks_Pathways/",thepathway,"_",classification,"_",nameout,"_pvals_network.txt",sep=""),stringsAsFactors = F,sep="\t")
  if(dim(netall1)[2]==15){
    netall2<-read.csv2(paste(experiment,"/OtherMetrics_Pathways/",thepathway,"_",classification,"_",nameout,"_pvals_network.txt",sep=""),stringsAsFactors = F,sep="\t")
    ## If both tables have the same number of interactions and the interactions in the same order:
    if(dim(netall1)[1]==dim(netall2)[1] && length(which(netall1[,1]==netall2[,1]))==dim(netall1)[1] && length(which(netall1[,2]==netall2[,2]))==dim(netall2)[1]){
      netall<-cbind(netall1,netall2[,3:dim(netall2)[2]])
      ## Remove the random values for the distance measures ##
      if(length(grep("Random",colnames(netall)))>0){netall<-netall[,-grep("Random",colnames(netall))]}
      ## Remove the real values for the distance measures ##
      originales<-c(grep("AllEuclidean",colnames(netall)),grep("IntEuclidean",colnames(netall)),grep("UniEuclidean",colnames(netall)),
                    grep("AllManhattan",colnames(netall)),grep("IntManhattan",colnames(netall)),grep("UniManhattan",colnames(netall)),
                    grep("AllCanberra",colnames(netall)),grep("IntCanberra",colnames(netall)),grep("UniCanberra",colnames(netall)))
      if(length(originales)>0){netall<-netall[,-originales]}
      ## Remove "Difference" from the colnames
      colnames(netall)<-gsub("Difference","",colnames(netall))
      write.table(netall,paste(experiment,"/Generated_Networks_Pathways/",thepathway,"_",classification,"_",nameout,"_pvals_network.txt",sep=""),sep="\t",row.names = F,quote=F)
      ## Correct for multipletesting ##
      if(length(grep("Pval",colnames(netall)))>0){
        for(z in grep("Pval",colnames(netall))){
          netall[,z]<-as.numeric(netall[,z])*length(netall[,z])
          netall[which(as.numeric(netall[,z])>1),z]<-1
        }
        colnames(netall)<-gsub("Pval","FDR",colnames(netall))
      }
      ## Write the tables ##
      write.table(netall,paste(experiment,"/Generated_Networks_Pathways/",thepathway,"_",classification,"_",nameout,"_fdr_network.txt",sep=""),quote=F,sep="\t",row.names=F)
      ## Convert the correlations and FDRs into 1s, 0s, and -1s ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
      for(z in grep("FDR",colnames(netall))){
        binarized<-rep(0,length(netall[,1]))
        binarized[intersect(which(as.numeric(netall[,z])<=0.05),which(as.numeric(netall[,z-1])<0))]<-(-1)
        binarized[intersect(which(as.numeric(netall[,z])<=0.05),which(as.numeric(netall[,z-1])>0))]<-1
        netall[,z-1]<-binarized
      }
      binary<-netall[,-grep("FDR",colnames(netall))]
      ## Write both tables ##
      colnames(binary)<-c("Disease1","Disease2","PearsonAll","SpearmanAll","PearsonInt","SpearmanInt","PearsonUni","SpearmanUni","Fisher",
                          "CosineSimilarityAll","CosineSimilarityInt","CosineSimilarityUni",
                          "EuclideanAll","EuclideanInt","EuclideanUni",
                          "ManhattanAll","ManhattanInt","ManhattanUni",
                          "CanberraAll","CanberraInt","CanberraUni")
      write.table(binary,paste(experiment,"/Generated_Networks_Pathways/",thepathway,"_",classification,"_",nameout,"_binarized_network.txt",sep=""),quote=F,sep="\t",row.names=F)
      ## Which interactions are different?
      resultstab<-c() ; for(z in 3:length(binary)){resultstab<-cbind(resultstab,length(which(binary[,z]!=0)))}
      resultstab<-cbind(resultstab,length(binary[,1]),length(unique(c(binary$Disease1,binary$Disease2))))
      colnames(resultstab)<-c(colnames(binary[3:length(colnames(binary))]),"PotentialInteractions","Diseases")
      resultstab<-rbind(resultstab,round(100*(resultstab[1,]/resultstab[1,11]),2))
      rownames(resultstab)<-c("All","Percentage_All")
      write.table(resultstab,paste(experiment,"/Generated_Networks_Pathways/",classification,"_",nameout,"_number_interactions.txt",sep=""),quote=F,sep="\t")
      resultados<-list("resultstab"=resultstab,"binarynet"=binary,"fdrnet"=netall)
    }
    ## If both tables don't have the same number of interactions and the interactions in the same order:
    if((dim(netall1)[1]==dim(netall2)[1] && length(which(netall1[,1]==netall2[,1]))==dim(netall1)[1] && length(which(netall1[,2]==netall2[,2]))==dim(netall2)[1])==FALSE){
      print("Error! Tables have different number of interactions or are ordered in a different way")
    }
  }
  if(dim(netall1)[2]!=15){
    resultados<-c()
  }
  print("Finished!")
  return(resultados)
}
if(as.numeric(args[1])==3){
  paraque<-setdiff(names(reactomeparentgenelist),"Digestion and absorption")
  for(a in 1:length(paraque)){
    ## ICD10 ##
    microicd10adjusted<-correctformultipletesting("Adjusted","ICD10","Microarrays",gsub(" ","_",paraque[a]))
    microicd10women<-correctformultipletesting("Women","ICD10","Microarrays",gsub(" ","_",paraque[a]))
    microicd10men<-correctformultipletesting("Men","ICD10","Microarrays",gsub(" ","_",paraque[a]))
    # microicd10cases<-correctformultipletesting("Cases","ICD10","Microarrays",gsub(" ","_",paraque[a]))
    # microicd10controls<-correctformultipletesting("Controls","ICD10","Microarrays",gsub(" ","_",paraque[a]))
    
    ## ICD9 ##
    microicd9adjusted<-correctformultipletesting("Adjusted","ICD9","Microarrays",gsub(" ","_",paraque[a]))
    microicd9women<-correctformultipletesting("Women","ICD9","Microarrays",gsub(" ","_",paraque[a]))
    microicd9men<-correctformultipletesting("Men","ICD9","Microarrays",gsub(" ","_",paraque[a]))
    # microicd9cases<-correctformultipletesting("Cases","ICD9","Microarrays",gsub(" ","_",paraque[a]))
    # microicd9controls<-correctformultipletesting("Controls","ICD9","Microarrays",gsub(" ","_",paraque[a]))
    print(a)
  }
}

#### Plot the correlation between number of significant interactions and the size of the pathway ####
if(args[1]=="plots"){
  adj10<-list.files("Microarrays/Generated_Networks_Pathways/")[grep("ICD10_Adjusted_binarized",list.files("Microarrays/Generated_Networks_Pathways/"))]
  adj10names<-gsub("_"," ",gsub("_ICD10_Adjusted_binarized_network.txt","",adj10))
  numbergenespath<-c() ; for(a in 1:length(adj10names)){numbergenespath<-c(numbergenespath,length(reactomeparentgenelist[[adj10names[a]]]))}
  names(numbergenespath)<-adj10names
  pdf(file="Microarrays/Plots/Correlation_genesperpathway_significantinteractions.pdf")
  correls<-c()
  for(b in 3:21){ ## Change this to 21 once we get all the similarity measures 
    positivos<-c() ; negativos<-c()
    for(a in 1:length(adj10)){
      tt<-read.csv2(paste("Microarrays/Generated_Networks_Pathways/",adj10[a],sep=""),stringsAsFactors = F,sep="\t")
      positivos<-c(positivos,length(which(tt[,b]==1)))
      negativos<-c(negativos,length(which(tt[,b]==(-1))))
    }
    plot(positivos,numbergenespath,xlab="Positive Similarities",ylab="Number of genes per pathway",main=paste("Number of significant similarities\nby Nº of genes per parent pathway","\n",colnames(tt)[b],sep=""))
    text(300,2500,paste("cor = ",round(as.numeric(cor.test(positivos,numbergenespath)$estimate),3),"\npval = ",as.numeric(cor.test(positivos,numbergenespath)$p.value),sep=""),cex=0.6)
    correls<-rbind(correls,c(as.numeric(cor.test(positivos,numbergenespath)$estimate),as.numeric(cor.test(positivos,numbergenespath)$p.value)))
  }
  dev.off()
}

# ## Generate the myrun file ##
# myrun<-paste("Rscript 06_2_Generate_disease_multilayer_networks_OtherMetrics.R 2",
#              c(paste("Adjusted","ICD9","Microarrays"),paste("Adjusted","ICD10","Microarrays"),
#                paste("Women","ICD9","Microarrays"),paste("Men","ICD9","Microarrays"),paste("Women","ICD10","Microarrays"),paste("Men","ICD10","Microarrays"),
#                paste("Cases","ICD9","Microarrays"),paste("Controls","ICD9","Microarrays"),paste("Cases","ICD10","Microarrays"),paste("Controls","ICD10","Microarrays")))
# myrun2<-c()
# for(a in 1:length(myrun)){
#   # a<-1
#   myrun2<-c(myrun2,paste(myrun[a],gsub(" ","_",names(reactomeparentgenelist)),sep=" "))
# }
# write.table(myrun2,"my_run_05_4_CalculateOtherSimilarities_byPathCategory.txt",quote=F,sep="\t",row.names=F,col.names=F)



