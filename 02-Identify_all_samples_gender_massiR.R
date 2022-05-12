## CODE FOR THE DIFFERENTIAL EXPRESSION ANALYSIS OF EACH CASE SAMPLE AGAINST ALL THE CONTROL SAMPLES FROM THE SAME STUDY ##
## Developed by Jon Sanchez-Valle
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es

## load the essential packages 
library("affy")
library("frma")
library("limma")
require("graphics")
library("hgu133plus2frmavecs")
library("massiR")
library("hgu133plus2.db")
library("cluster")
library("NbClust")
library("e1071")
## This package is not working ##
# library("genefilter")

#### Test if MassiR classifies better the samples separating cases and controls, both together, and combining all the studies ####
fileswithknownsex<-list.files("Microarrays/Data/Raw_data/")
## Extract the Y chromose probes from the massiR packages and identify the ones for HG U 133 plus2 ##
hgu133plus2yprobes <- data.frame(y.probes["affy_hg_u133_plus_2"])
theprobes<-rownames(hgu133plus2yprobes)
bystudyresults<-c()
for(a in 1:length(fileswithknownsex)){
  ## Load expression information ##
  control<-list.files(paste("./Microarrays/Data/Raw_data/",fileswithknownsex[a],"/Control",sep=""))
  path_con<-paste("./Microarrays/Data/Raw_data/",fileswithknownsex[a],"/Control/",control,sep="")
  case<-list.files(paste("./Data/Raw_data/",fileswithknownsex[a],"/Case",sep=""))
  path_ca<-paste("./Microarrays/Data/Raw_data/",fileswithknownsex[a],"/Case/",case,sep="")
  path_read<-c(path_con, path_ca)
  controls<-cbind(control,rep("control",length(control)))
  cases<-cbind(case,rep("case",length(case)))
  targets<-rbind(controls, cases)
  affyBatch<-ReadAffy(filenames=path_read)
  ## Normalize expression
  expSet<-frma(affyBatch)
  ## A data.frame
  comoda<-t(as.data.frame(expSet))
  guayes<-comoda[paste("X",theprobes,sep=""),]
  tguayes<-t(guayes) ; tguayes<-cbind(rownames(tguayes),tguayes)
  write.table(tguayes,"Microarrays/Intermediate_results/Y_probes_by_study_normalized_all.txt",quote=F,sep="\t",col.names=F,row.names=F,append = T)
  print(a/length(fileswithknownsex))
}



#### Read the obtained results ####
## Read the matrix with the expression values of the Y chromosome probes for all the samples with known gender
expmat<-read.csv2("Microarrays/Intermediate_results/Y_probes_by_study_normalized_all.txt",stringsAsFactors = F,sep="\t",header=F)
expmat<-t(expmat)
samples<-expmat[1,]
## Add numbers
expressionmatrix<-c() ; for(a in 2:length(expmat[,1])){expressionmatrix<-rbind(expressionmatrix,as.numeric(expmat[a,]))}
colnames(expressionmatrix)<-paste("Sample_",1:length(samples),"_____",samples,sep="")
expressionmatrix<-as.data.frame(expressionmatrix)
# rownames(expressionmatrix)<-theprobes
## Extract samples' gender ##
# massi.y.out <-massi_y(expressionmatrix,hgu133plus2yprobes)
# massi_y_plot(massi.y.out)
# massi.select.out <- massi_select(expressionmatrix, hgu133plus2yprobes, threshold=3)
# eset.results1 <-massi_cluster(massi.select.out)
eset.results1 <-massi_cluster(expressionmatrix)
extractedgender1<-eset.results1$massi.results
extractedgender1$ID<-gsub(".+_____","",extractedgender1$ID)

## The known ones
fileswithknownsex<-gsub(".txt","",list.files("Microarrays/Data/Raw_data_sex/"))
allknownsex<-c()
for(a in 1:length(fileswithknownsex)){
  ## Load the known sex
  knownsex<-read.csv2(paste("Microarrays/Data/Raw_data_sex/",fileswithknownsex[a],".txt",sep=""),stringsAsFactors = F,sep="\t")
  allknownsex<-rbind(allknownsex,knownsex)
}
# write.table(allknownsex,"Microarrays/Data/Sex_from_all_samples_with_known_sex.txt",quote=F,sep="\t",row.names = F)
## Add a column with the known sex ##
theknownone<-c()
for(a in 1:length(extractedgender1[,1])){
  # a<-1
  cual<-which(allknownsex$Sample==extractedgender1[a,1])
  if(length(cual)==0){theknownone<-c(theknownone,NA)}
  if(length(cual)>0){theknownone<-c(theknownone,allknownsex$Gender[cual][1])}
}
matched<-cbind(extractedgender1,theknownone)
## Remove the 1s and 0s and the NAs
eliminar<-c(which(is.na(matched$theknownone)),which(matched$theknownone==1),which(matched$theknownone==0))
nmatched<-matched[-eliminar,]
## Percentage of properly classified samples ##
(length(which(nmatched$sex==nmatched$theknownone))/length(nmatched$theknownone))*100 ## 94.09%
malclassificados<-nmatched[which(nmatched$sex!=nmatched$theknownone),]
malclassificados<-malclassificados[order(abs(malclassificados$z_score),decreasing = T),]


#### Create the new "Raw_data_predictedsex" tables with the predicted sex information ####
extractedgender<-extractedgender1[-which(duplicated(extractedgender1)),]
if("Raw_data_predictedsex"%in%list.files("Microarrays/Data/")==FALSE){dir.create("Microarrays/Data/Raw_data_predictedsex")}
datos<-list.files("Microarrays/Data/Raw_data/")
for(a in 1:length(datos)){
  # a<-1
  casos<-list.files(paste("Microarrays/Data/Raw_data/",datos[a],"/Case",sep=""))
  controles<-list.files(paste("Microarrays/Data/Raw_data/",datos[a],"/Control",sep=""))
  sexocasos<-c()
  for(b in 1:length(casos)){
    sexocasos<-rbind(sexocasos,c(casos[b],"Case",extractedgender$sex[which(extractedgender$ID==casos[b])[1]]))
    if(length(which(extractedgender$ID==casos[b]))==0){print(paste(datos[a],casos[b],sep=" - "))}
  }
  sexocontroles<-c()
  for(b in 1:length(controles)){
    sexocontroles<-rbind(sexocontroles,c(controles[b],"Control",extractedgender$sex[which(extractedgender$ID==controles[b])[1]]))
    if(length(which(extractedgender$ID==controles[b]))==0){print(paste(datos[a],a,controles[b],b,sep=" - "))}
  }
  sexo<-rbind(sexocasos,sexocontroles)
  colnames(sexo)<-c("Sample","Category","Gender")
  write.table(sexo,paste("Microarrays/Data/Raw_data_predictedsex/",datos[a],".txt",sep=""),quote=F,sep="\t",row.names=F)
  print(a/length(datos))
}








