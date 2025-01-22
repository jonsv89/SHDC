## CODE FOR THE GENERATION OF THE DISEASE NETWORKS FOR MEN, WOMEN, AND ADJUSTING BY SEX BASED ON GENES SEPARATELY FOR EACH REACTOME CATEGORY ##
## Developed by Jon Sanchez-Valle & Beatriz Urda-Garcia
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es & beatriz.urda@bsc.es

args = commandArgs(trailingOnly=TRUE)

#### Function for calculating similarities between diseases, adjusting by sex, separately for men and women, and between diseases and healthy differences women vs. men ####
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
colnames(genesperparent)<-c("Category","NumberGenes")
write.table(genesperparent,"Microarrays/Data/Number_of_genes_per_reactome_parents.txt",quote=F,sep="\t",row.names = F)

## Add all the genes associated with the mitochondria ##
if(length(grep("Mitochondria All",names(reactomeparentgenelist)))==0){
  mitogenes<-read.csv2("Microarrays/Data/MitoCarta.csv",stringsAsFactors = F)
  reactomeparentgenelist[["Mitochondria All"]]<-toupper(unique(mitogenes$Symbol))
}
names(reactomeparentgenelist)<-gsub("\\)","",gsub("\\(","",names(reactomeparentgenelist)))

# nameout<-"Adjusted"
# classification<-"ICD10"
# experiment<-"Microarrays"
# thepathway<-names(reactomeparentgenelist)[1]
# thegene<-reactomeparentgenelist[[1]]
build_network<-function(nameout,classification,experiment,thepathway,thegene){
  if("Generated_Networks_Pathways"%in%list.files(experiment)==FALSE){dir.create(paste(experiment,"/Generated_Networks_Pathways",sep=""))}
  ## Start the analysis!! ##
  ## @@ @@ @@ @@ @@ @@ @@ ##
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
  }
  ## Similarities between pairs of diseases ##
  netall<-c() ; fishtab<-c()
  for(a in 1:(length(names(vsdegs))-1)){
    for(b in (a+1):length(names(vsdegs))){
      ## Over all the genes ##
      ## @@ @@ @@  @@ @@ @@ ##
      commongenes<-intersect(names(allgenes[[a]]),names(allgenes[[b]]))
      if(length(commongenes)>3){
        femcorrelp<-cor.test(as.numeric(as.character(allgenes[[a]][commongenes])),as.numeric(as.character(allgenes[[b]][commongenes])),method = "pearson")
        femcorrels<-cor.test(as.numeric(as.character(allgenes[[a]][commongenes])),as.numeric(as.character(allgenes[[b]][commongenes])),method = "spearman")
      }
      if(length(commongenes)<=3){
        femcorrelp<-list("estimate"=NA,"p.value"=NA)
        femcorrels<-list("estimate"=NA,"p.value"=NA)
      }
      ## Over the intersection of the sDEGs ##
      ## @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ ##
      unof<-vsdegs[[a]] ; dosf<-vsdegs[[b]]
      allunof<-c(unof$Up,unof$Down) ; alldosf<-c(dosf$Up,dosf$Down)
      fintersection<-intersect(names(allunof),names(alldosf))
      if(length(fintersection)>3){
        ifemcorrelp<-cor.test(as.numeric(as.character(allunof[fintersection])),as.numeric(as.character(alldosf[fintersection])),method = "pearson")
        ifemcorrels<-cor.test(as.numeric(as.character(allunof[fintersection])),as.numeric(as.character(alldosf[fintersection])),method = "spearman")
      }
      if(length(fintersection)<=3){
        ifemcorrelp<-list("estimate"=NA,"p.value"=NA)
        ifemcorrels<-list("estimate"=NA,"p.value"=NA)
      }
      ## Over the union of sDEGs ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ ##
      funion<-unique(c(names(allunof),names(alldosf)))
      numberifunion<-length(intersect(names(which(is.na(allgenes[[a]][funion])==FALSE)),names(which(is.na(allgenes[[b]][funion])==FALSE))))
      if(numberifunion<=3){
        ufemcorrelp<-list("estimate"=NA,"p.value"=NA)
        ufemcorrels<-list("estimate"=NA,"p.value"=NA)
      }
      if(numberifunion>3){
        ufemcorrelp<-cor.test(as.numeric(as.character(allgenes[[a]][funion])),as.numeric(as.character(allgenes[[b]][funion])),method = "pearson")
        ufemcorrels<-cor.test(as.numeric(as.character(allgenes[[a]][funion])),as.numeric(as.character(allgenes[[b]][funion])),method = "spearman")
      }
      ## Calculating Fisher's overlaps ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
      background<-length(commongenes)
      positi1<-intersect(names(vsdegs[[a]]$Up),commongenes) ; negati1<-intersect(names(vsdegs[[a]]$Down),commongenes) ; pega<-c() ; pegados<-c()
      positi2<-intersect(names(vsdegs[[b]]$Up),commongenes) ; negati2<-intersect(names(vsdegs[[b]]$Down),commongenes)
      interspd <- length(intersect(positi1,positi2))
      kkpd <- matrix(c(interspd,length(positi1)-interspd,length(positi2)-interspd,background+interspd-length(positi1)-length(positi2)),nrow=2,ncol=2)
      fispd<-fisher.test(kkpd,alternative="greater") $p.value
      intersnd <- length(intersect(negati1,negati2))
      kknd <- matrix(c(intersnd,length(negati1)-intersnd,length(negati2)-intersnd,background+intersnd-length(negati1)-length(negati2)),nrow=2,ncol=2)
      fisnd<-fisher.test(kknd,alternative="greater") $p.value
      interspi <- length(intersect(positi1,negati2))
      kkpi <- matrix(c(interspi,length(positi1)-interspi,length(negati2)-interspi,background+interspi-length(positi1)-length(negati2)),nrow=2,ncol=2)
      fispi<-fisher.test(kkpi,alternative="greater") $p.value
      intersni <- length(intersect(negati1,positi2))
      kkni <- matrix(c(intersni,length(negati1)-intersni,length(positi2)-intersni,background+intersni-length(negati1)-length(positi2)),nrow=2,ncol=2)
      fisni<-fisher.test(kkni,alternative="greater") $p.value
      fjunt<-c(fispi,fispd,fisni,fisnd)
      ## Put all the correlations together ##
      ## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
      netall<-rbind(netall,c(names(allgenes)[c(a,b)],as.numeric(femcorrelp$estimate),as.numeric(femcorrelp$p.value),
                                         as.numeric(femcorrels$estimate),as.numeric(femcorrels$p.value),
                                         as.numeric(ifemcorrelp$estimate),as.numeric(ifemcorrelp$p.value),
                                         as.numeric(ifemcorrels$estimate),as.numeric(ifemcorrels$p.value),
                                         as.numeric(ufemcorrelp$estimate),as.numeric(ufemcorrelp$p.value),
                                         as.numeric(ufemcorrels$estimate),as.numeric(ufemcorrels$p.value)))
      fishtab<-rbind(fishtab,fjunt)
    }
    # print(paste(round((a/length(names(allgenes)))*100,2),"%",sep=""))
  }
  fishtab2<-fishtab
  ## Transform the fishtab into a 1, 0, and -1 column
  for(a in 1:4){fishtab[,a]<-fishtab[,a]*length(fishtab[,1])*4}
  for(a in 1:4){fishtab[which(fishtab[,a]>1),a]<-1}
  threshold<-0.05
  negs<-list("Pos_Neg"=which(fishtab[,1]<=threshold),"Neg_Pos"=which(fishtab[,3]<=threshold),"Pos_Pos"=which(fishtab[,2]>threshold),"Neg_Neg"=which(fishtab[,4]>threshold))
  poss<-list("Pos_Neg"=which(fishtab[,1]>threshold),"Neg_Pos"=which(fishtab[,3]>threshold),"Pos_Pos"=which(fishtab[,2]<=threshold),"Neg_Neg"=which(fishtab[,4]<=threshold))
  ffishint<-rep(0,length(fishtab[,1])) ; ffishint[Reduce(intersect,negs)]<--1 ; ffishint[Reduce(intersect,poss)]<-1
  ## Add the fisher's test columns ##
  if(length(netall[1,])==14){netall<-cbind(netall,ffishint)}
  ## Add column names ##
  colnames(netall)<-c("Disease1","Disease2","AllPearson","AllPvalPearson","AllSpearman","AllPvalSpearman",
                            "IntPearson","IntPvalPearson","IntSpearman","IntPvalSpearman",
                            "UniPearson","UniPvalPearson","UniSpearman","UniPvalSpearman","FisherTests")
  ## Write the tables ##
  write.table(netall,paste(experiment,"/Generated_Networks_Pathways/",gsub(" ","_",thepathway),"_",classification,"_",nameout,"_pvals_network.txt",sep=""),quote=F,sep="\t",row.names=F)
  ## Correct for multipletesting ##
  if(length(grep("Pval",colnames(netall)))>0){
    netall[,4]<-as.numeric(netall[,4])*length(netall[,1])
    netall[,6]<-as.numeric(netall[,6])*length(netall[,1])
    netall[,8]<-as.numeric(netall[,8])*length(netall[,1])
    netall[,10]<-as.numeric(netall[,10])*length(netall[,1])
    netall[,12]<-as.numeric(netall[,12])*length(netall[,1])
    netall[,14]<-as.numeric(netall[,14])*length(netall[,1])
    colnames(netall)<-c("Disease1","Disease2","AllPearson","AllFDRPearson","AllSpearman","AllFDRSpearman",
                              "IntPearson","IntFDRPearson","IntSpearman","IntFDRSpearman",
                              "UniPearson","UniFDRPearson","UniSpearman","UniFDRSpearman","FisherTests")
    netall[which(as.numeric(netall[,4])>1),4]<-1
    netall[which(as.numeric(netall[,6])>1),6]<-1
    netall[which(as.numeric(netall[,8])>1),8]<-1
    netall[which(as.numeric(netall[,10])>1),10]<-1
    netall[which(as.numeric(netall[,12])>1),12]<-1
    netall[which(as.numeric(netall[,14])>1),14]<-1
  }
  ## Write the tables ##
  write.table(netall,paste(experiment,"/Generated_Networks_Pathways/",gsub(" ","_",thepathway),"_",classification,"_",nameout,"_fdr_network.txt",sep=""),quote=F,sep="\t",row.names=F)
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
  binary<-cbind(netall[,1:2],palldir,salldir,pintdir,sintdir,punidir,sunidir,netall[,15])
  ## Write both tables ##
  colnames(binary)<-c("Disease1","Disease2","PearsonAll","SpearmanAll","PearsonInt","SpearmanInt","PearsonUni","SpearmanUni","Fisher")
  write.table(binary,paste(experiment,"/Generated_Networks_Pathways/",gsub(" ","_",thepathway),"_",classification,"_",nameout,"_binarized_network.txt",sep=""),quote=F,sep="\t",row.names=F)
  ## Which interactions are different?
  resultstab<-rbind(cbind(c(length(which(binary[,3]!=0))),
                          c(length(which(binary[,4]!=0))),
                          c(length(which(binary[,5]!=0))),
                          c(length(which(binary[,6]!=0))),
                          c(length(which(binary[,7]!=0))),
                          c(length(which(binary[,8]!=0))),
                          c(length(which(binary[,9]!=0)))))
  resultstab<-cbind(resultstab,length(binary[,1]),length(names(allgenes)))
  colnames(resultstab)<-c("PearsonAll","SpearmanAll","PearsonInt","SpearmanInt","PearsonUni","SpearmanUni","Fisher","PotentialInteractions","Diseases")
  resultstab<-rbind(resultstab,round(100*(resultstab[1,]/resultstab[1,8]),2))
  rownames(resultstab)<-c("All","Percentage_All")
  write.table(resultstab,paste(experiment,"/Generated_Networks_Pathways/",gsub(" ","_",thepathway),"_",classification,"_",nameout,"_number_interactions.txt",sep=""),quote=F,sep="\t")
  resultados<-list("resultstab"=resultstab,"binarynet"=binary,"fdrnet"=netall)
  print(paste(thepathway,": ",classification," ",nameout," Finished!",sep=""))
  return(resultados)
}

#### Microarrays ####
## @@ @@ @ @ @@ @@ ##
for(tt in 1:length(names(reactomeparentgenelist))){
  if(length(reactomeparentgenelist[[tt]])>=10){
    thepathway<-names(reactomeparentgenelist)[tt]
    thegene<-reactomeparentgenelist[[tt]]
    ## We add an if to ensure that if the network has been already generated is not generated again ##
    if(length(grep(paste(gsub(" ","_",thepathway),"_ICD9_Adjusted",sep=""),list.files("Microarrays/Generated_Networks_Pathways/")))==0){
      adj9net<-build_network("Adjusted","ICD9","Microarrays",thepathway,thegene)
    }
    if(length(grep(paste(gsub(" ","_",thepathway),"_ICD10_Adjusted",sep=""),list.files("Microarrays/Generated_Networks_Pathways/")))==0){
      adj10net<-build_network("Adjusted","ICD10","Microarrays",thepathway,thegene)
    }
    if(length(grep(paste(gsub(" ","_",thepathway),"_ICD9_Women",sep=""),list.files("Microarrays/Generated_Networks_Pathways/")))==0){
      women9net<-build_network("Women","ICD9","Microarrays",thepathway,thegene)
    }
    if(length(grep(paste(gsub(" ","_",thepathway),"_ICD10_Women",sep=""),list.files("Microarrays/Generated_Networks_Pathways/")))==0){
      women10net<-build_network("Women","ICD10","Microarrays",thepathway,thegene)
    }
    if(length(grep(paste(gsub(" ","_",thepathway),"_ICD9_Men",sep=""),list.files("Microarrays/Generated_Networks_Pathways/")))==0){
      men9net<-build_network("Men","ICD9","Microarrays",thepathway,thegene)
    }
    if(length(grep(paste(gsub(" ","_",thepathway),"_ICD10_Men",sep=""),list.files("Microarrays/Generated_Networks_Pathways/")))==0){
      men10net<-build_network("Men","ICD10","Microarrays",thepathway,thegene)
    }
    if(length(grep(paste(gsub(" ","_",thepathway),"_ICD9_Cases",sep=""),list.files("Microarrays/Generated_Networks_Pathways/")))==0){
      cases9net<-build_network("Cases","ICD9","Microarrays",thepathway,thegene)
    }
    if(length(grep(paste(gsub(" ","_",thepathway),"_ICD10_Cases",sep=""),list.files("Microarrays/Generated_Networks_Pathways/")))==0){
      cases10net<-build_network("Cases","ICD10","Microarrays",thepathway,thegene)
    }
    if(length(grep(paste(gsub(" ","_",thepathway),"_ICD9_Controls",sep=""),list.files("Microarrays/Generated_Networks_Pathways/")))==0){
      controls9net<-build_network("Controls","ICD9","Microarrays",thepathway,thegene)
    }
    if(length(grep(paste(gsub(" ","_",thepathway),"_ICD10_Controls",sep=""),list.files("Microarrays/Generated_Networks_Pathways/")))==0){
      controls10net<-build_network("Controls","ICD10","Microarrays",thepathway,thegene)
    }
  }
  print("")
  print(tt)
  print("")
}
print("Networks generated")






















