## CODE FTO IDENTIFY THE SEX OF SAMPLES AND ELIMINATE LOW QUALITY SAMPLES ##
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
library("genefilter")

#### Identify samples' sex ####
## @@ @@ @@ @@ @ @@ @@ @@ @@ ##
hgu133plus2yprobes <- data.frame(y.probes["affy_hg_u133_plus_2"])
theprobes<-rownames(hgu133plus2yprobes)
## Get the expression of the probes analyzing y chromosome genes ##
if("Y_probes_by_study_normalized_all_samples.txt"%in%list.files("Microarrays/IntermediateData/")==FALSE){
  fileswithknownsex<-list.files("Microarrays/Data/Raw_data/")
  ## Extract the Y chromose probes from the massiR packages and identify the ones for HG U 133 plus2 ##
  hgu133plus2yprobes <- data.frame(y.probes["affy_hg_u133_plus_2"])
  theprobes<-rownames(hgu133plus2yprobes)
  bystudyresults<-c()
  for(a in 1:length(fileswithknownsex)){
    ## Load expression information ##
    control<-list.files(paste("./Microarrays/Data/Raw_data/",fileswithknownsex[a],"/Control",sep=""))
    path_con<-paste("./Microarrays/Data/Raw_data/",fileswithknownsex[a],"/Control/",control,sep="")
    case<-list.files(paste("./Microarrays/Data/Raw_data/",fileswithknownsex[a],"/Case",sep=""))
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
    write.table(tguayes,"Microarrays/IntermediateData/Y_probes_by_study_normalized_all_samples.txt",quote=F,sep="\t",col.names=F,row.names=F,append = T)
    print(paste(fileswithknownsex[a],": ",round((a/length(fileswithknownsex))*100,2),"%",sep=""))
  }
}

## Read the matrix with the expression values of the Y chromosome probes for all the samples
expmat<-read.csv2("Microarrays/IntermediateData/Y_probes_by_study_normalized_all_samples.txt",stringsAsFactors = F,sep="\t",header=F)
expmat<-t(expmat)
samples<-expmat[1,]
## Add numbers
expressionmatrix<-c() ; for(a in 2:length(expmat[,1])){expressionmatrix<-rbind(expressionmatrix,as.numeric(expmat[a,]))}
colnames(expressionmatrix)<-paste("Sample_",1:length(samples),"_____",samples,sep="")
expressionmatrix<-as.data.frame(expressionmatrix)
eset.results1 <-massi_cluster(expressionmatrix)
extractedgender1<-eset.results1$massi.results
extractedgender1$ID<-gsub(".+_____","",extractedgender1$ID)

## The known ones ##
fileswithknownsex<-gsub(".txt","",list.files("Microarrays/IntermediateData/Raw_data_sex/"))
allknownsex<-c()
for(a in 1:length(fileswithknownsex)){
  ## Load the known sex
  knownsex<-read.csv2(paste("Microarrays/IntermediateData/Raw_data_sex/",fileswithknownsex[a],".txt",sep=""),stringsAsFactors = F,sep="\t")
  allknownsex<-rbind(allknownsex,knownsex)
}
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
(length(which(nmatched$sex==nmatched$theknownone))/length(nmatched$theknownone))*100
malclassificados<-nmatched[which(nmatched$sex!=nmatched$theknownone),]
malclassificados<-malclassificados[order(abs(malclassificados$z_score),decreasing = T),]

## Create the new "Raw_data_predictedsex" tables with the predicted sex information ##
extractedgender<-extractedgender1[-which(duplicated(extractedgender1)),]
if("Raw_data_predictedsex"%in%list.files("Microarrays/IntermediateData/")==FALSE){dir.create("Microarrays/IntermediateData/Raw_data_predictedsex")}
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
  write.table(sexo,paste("Microarrays/IntermediateData/Raw_data_predictedsex/",datos[a],".txt",sep=""),quote=F,sep="\t",row.names=F)
  print(paste(round((a/length(datos))*100,2),"%",sep=""))
}


#### Identify and remove low quality samples ####
## @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ ##
## Create a directory to write the table indicating the samples to be removed ##
if("Targets"%in%list.files("Microarrays/Data/")==FALSE){dir.create("Microarrays/Data/Targets")}
## Load each study and identify low quality samples ##
fileswithknownsex<-list.files("Microarrays/IntermediateData/Raw_data_predictedsex/")
fileswithknownsex<-setdiff(fileswithknownsex,list.files("Microarrays/IntermediateData/Targets/"))
for(a in 1:length(fileswithknownsex)){
  # a<-1
  sextab<-read.csv2(paste("Microarrays/IntermediateData/Raw_data_predictedsex/",fileswithknownsex[a],sep=""),stringsAsFactors = F,sep="\t")
  ## Create the targets table ##
  targets<-cbind(paste("Microarrays/Data/Raw_data/",gsub(".txt","",fileswithknownsex[a]),"/",sextab$Category,"/",sextab$Sample,sep=""),sextab,gsub(".+_","",gsub(".txt","",fileswithknownsex[a])))
  colnames(targets)<-c("Path","Sample","Category","Sex","Study")
  ## Read the affy files ##
  affyBatch<-ReadAffy(filenames=targets$Path)
  ## Normalize ##
  expSet<-frma(affyBatch)
  ## Calculate GNUSE median values for each sample, and remove those with a value higher than 1.25
  ## What means that you remove those arrays whose precission is, by mean, a 25% worst that the usual array
  gnuseres<-GNUSE(expSet,type="stats")
  gnuseres<-gnuseres[,targets$Sample]
  removeornot<-rep("no",length(targets$Sample))
  if(length(which(gnuseres[1,]>1.25))>0){removeornot[which(gnuseres[1,]>1.25)]<-"yes"}
  targets<-cbind(targets,removeornot)
  colnames(targets)[6]<-"RemoveThem"
  targets<-cbind(targets,gnuseres[1,])
  colnames(targets)[7]<-"MedianGNUSE"
  write.table(targets,paste("Microarrays/IntermediateData/Targets/",fileswithknownsex[a],sep=""),quote=F,sep="\t",row.names=F)
  print(paste(round((a/length(fileswithknownsex))*100,2),"%",sep=""))
}
files<-list.files("Microarrays/IntermediateData/Targets/")
resumen<-c() ; for(a in 1:length(files)){targets<-read.csv2(paste("Microarrays/IntermediateData/Targets/",files[a],sep=""),stringsAsFactors = F,sep="\t") ; resumen<-rbind(resumen,targets)}
pdf(file="Microarrays/Plots/Histogram_GNUSE.pdf")
  hist(as.numeric(resumen$MedianGNUSE),breaks=100,xlim=c(0.5,2.5),xlab="median GNUSE",main="Histogram of median GNUSE\nmicroarray HG-U-133 plus 2")
dev.off()

## Remove the samples that should not be used ##
## @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ ##
if("FinalTargets"%in%list.files("Microarrays/Data/")==FALSE){dir.create("Microarrays/Data/FinalTargets")}
targetdatasets<-list.files("Microarrays/IntermediateData/Targets/")
for(a in 1:length(targetdatasets)){
  # a<-1
  tt<-read.csv2(paste("Microarrays/IntermediateData/Targets/",targetdatasets[a],sep=""),stringsAsFactors = F,sep="\t")
  if(length(which(tt$RemoveThem=="yes"))>0){tt<-tt[-which(tt$RemoveThem=="yes"),]}
  if(dim(tt)[1]>0){
    write.table(tt,paste("Microarrays/Data/FinalTargets/",targetdatasets[a],sep=""),quote=F,sep="\t",row.names=F)
  }
  if(dim(tt)[1]==0){print(paste("We lost the study:",targetdatasets[a]))}
}


#### Extract the sex of newly added samples ####
## @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ ##
## First load the matrix with the expression of all the samples ##
tabla<-read.csv2("Microarrays/IntermediateData/Y_probes_by_study_normalized_all_samples.txt",stringsAsFactors = F,sep="\t",header=F)
if("Y_probes_by_study_normalized_all_samples_new.txt"%in%list.files("Microarrays/IntermediateData/")==FALSE){
  write.table(tabla,"Microarrays/IntermediateData/Y_probes_by_study_normalized_all_samples_new.txt",quote=F,sep="\t",col.names=F,row.names=F)
}

if(length(setdiff(list.files("Microarrays/Data/Raw_data/"),gsub(".txt","",list.files("Microarrays/IntermediateData/Raw_data_predictedsex/"))))>0){
  fileswithknownsex<-setdiff(list.files("Microarrays/Data/Raw_data/"),gsub(".txt","",list.files("Microarrays/IntermediateData/Raw_data_predictedsex/")))
  ## Extract the Y chromose probes from the massiR packages and identify the ones for HG U 133 plus2 ##
  hgu133plus2yprobes <- data.frame(y.probes["affy_hg_u133_plus_2"])
  theprobes<-rownames(hgu133plus2yprobes)
  bystudyresults<-c()
  for(a in 1:length(fileswithknownsex)){
    ## Load expression information ##
    control<-list.files(paste("./Microarrays/Data/Raw_data/",fileswithknownsex[a],"/Control",sep=""))
    path_con<-paste("./Microarrays/Data/Raw_data/",fileswithknownsex[a],"/Control/",control,sep="")
    case<-list.files(paste("./Microarrays/Data/Raw_data/",fileswithknownsex[a],"/Case",sep=""))
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
    write.table(tguayes,"Microarrays/IntermediateData/Y_probes_by_study_normalized_all_samples_new.txt",quote=F,sep="\t",col.names=F,row.names=F,append = T)
    print(paste(fileswithknownsex[a],": ",round((a/length(fileswithknownsex))*100,2),"%",sep=""))
  }
  
  ## Read the matrix with the expression values of the Y chromosome probes for all the samples with known gender
  expmat<-read.csv2("Microarrays/IntermediateData/Y_probes_by_study_normalized_all_samples_new.txt",stringsAsFactors = F,sep="\t",header=F)
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
  
  #### Create the new "Raw_data_predictedsex" tables with the predicted sex information ####
  extractedgender<-extractedgender1[-which(duplicated(extractedgender1)),]
  for(a in 1:length(fileswithknownsex)){
    # a<-1
    casos<-list.files(paste("Microarrays/Data/Raw_data/",fileswithknownsex[a],"/Case",sep=""))
    controles<-list.files(paste("Microarrays/Data/Raw_data/",fileswithknownsex[a],"/Control",sep=""))
    sexocasos<-c()
    for(b in 1:length(casos)){
      sexocasos<-rbind(sexocasos,c(casos[b],"Case",extractedgender$sex[which(extractedgender$ID==casos[b])[1]]))
      if(length(which(extractedgender$ID==casos[b]))==0){print(paste(fileswithknownsex[a],casos[b],sep=" - "))}
    }
    sexocontroles<-c()
    for(b in 1:length(controles)){
      sexocontroles<-rbind(sexocontroles,c(controles[b],"Control",extractedgender$sex[which(extractedgender$ID==controles[b])[1]]))
      if(length(which(extractedgender$ID==controles[b]))==0){print(paste(fileswithknownsex[a],a,controles[b],b,sep=" - "))}
    }
    sexo<-rbind(sexocasos,sexocontroles)
    colnames(sexo)<-c("Sample","Category","Gender")
    write.table(sexo,paste("Microarrays/IntermediateData/Raw_data_predictedsex/",fileswithknownsex[a],".txt",sep=""),quote=F,sep="\t",row.names=F)
    print(paste(round((a/length(fileswithknownsex))*100,2),"%",sep=""))
  }
  
  #### Identify and remove low quality samples ####
  ## @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ ##
  ## Load each study and identify low quality samples ##
  fileswithknownsex<-list.files("Microarrays/IntermediateData/Raw_data_predictedsex/")
  fileswithknownsex<-setdiff(fileswithknownsex,list.files("Microarrays/IntermediateData/Targets/"))
  for(a in 1:length(fileswithknownsex)){
    # a<-1
    sextab<-read.csv2(paste("Microarrays/IntermediateData/Raw_data_predictedsex/",fileswithknownsex[a],sep=""),stringsAsFactors = F,sep="\t")
    ## Create the targets table ##
    targets<-cbind(paste("Microarrays/Data/Raw_data/",gsub(".txt","",fileswithknownsex[a]),"/",sextab$Category,"/",sextab$Sample,sep=""),sextab,gsub(".+_","",gsub(".txt","",fileswithknownsex[a])))
    colnames(targets)<-c("Path","Sample","Category","Sex","Study")
    ## Read the affy files ##
    affyBatch<-ReadAffy(filenames=targets$Path)
    ## Normalize ##
    expSet<-frma(affyBatch)
    ## Calculate GNUSE median values for each sample, and remove those with a value higher than 1.25
    ## What means that you remove those arrays whose precission is, by mean, a 25% worst that the usual array
    gnuseres<-GNUSE(expSet,type="stats")
    gnuseres<-gnuseres[,targets$Sample]
    removeornot<-rep("no",length(targets$Sample))
    if(length(which(gnuseres[1,]>1.25))>0){removeornot[which(gnuseres[1,]>1.25)]<-"yes"}
    targets<-cbind(targets,removeornot)
    colnames(targets)[6]<-"RemoveThem"
    targets<-cbind(targets,gnuseres[1,])
    colnames(targets)[7]<-"MedianGNUSE"
    write.table(targets,paste("Microarrays/IntermediateData/Targets/",fileswithknownsex[a],sep=""),quote=F,sep="\t",row.names=F)
    print(paste(round((a/length(fileswithknownsex))*100,2),"%",sep=""))
  }
  files<-list.files("Microarrays/IntermediateData/Targets/")
  resumen<-c() ; for(a in 1:length(files)){targets<-read.csv2(paste("Microarrays/IntermediateData/Targets/",files[a],sep=""),stringsAsFactors = F,sep="\t") ; resumen<-rbind(resumen,targets)}
  
  ## Remove the samples that should not be used ##
  ## @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ ##
  for(a in 1:length(fileswithknownsex)){
    # a<-1
    tt<-read.csv2(paste("Microarrays/IntermediateData/Targets/",fileswithknownsex[a],sep=""),stringsAsFactors = F,sep="\t")
    if(length(which(tt$RemoveThem=="yes"))>0){tt<-tt[-which(tt$RemoveThem=="yes"),]}
    if(dim(tt)[1]>0){
      write.table(tt,paste("Microarrays/Data/FinalTargets/",fileswithknownsex[a],sep=""),quote=F,sep="\t",row.names=F)
    }
    if(dim(tt)[1]==0){print(paste("We lost the study:",fileswithknownsex[a]))}
  }
}

#### Create a supplementary dataset for the manuscript indicating the number of samples per disease and study ####
finaltargets<-list.files("Microarrays/Data/FinalTargets/")
datasets<-c()
for(a in 1:length(finaltargets)){
  # a<-1
  tt<-read.csv2(paste("Microarrays/Data/FinalTargets/",finaltargets[a],sep=""),stringsAsFactors = F,sep="\t")
  if(length(which(table(tt$Category)>=3))==2){
    datasets<-rbind(datasets,c(gsub("_.+","",finaltargets[a]),gsub(".txt","",gsub(".+_","",finaltargets[a])),
                               length(intersect(which(tt$Category=="Case"),which(tt$Sex=="female"))),length(intersect(which(tt$Category=="Control"),which(tt$Sex=="female"))),
                               length(intersect(which(tt$Category=="Case"),which(tt$Sex=="male"))),length(intersect(which(tt$Category=="Control"),which(tt$Sex=="male")))))
  }
}
## Load the name of the diseases to create the table ##
disinf<-read.csv2("Microarrays/Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
disname<-disinf$Disease_name ; names(disname)<-disinf$Disease
## Add the correct name of the disease ##
datasets[,1]<-disname[datasets[,1]]
## write the table ##
colnames(datasets)<-c("disease","study","cases (women)","controls (women)","cases (men)","controls (men)")
write.table(datasets,"Submission/Supplementary_Dataset_1_samples_per_study_and_disease.txt",quote=F,sep="\t",row.names=F)

#### Create a supplementary table with information on the categorias and ICD codes for each disease ####
disinf<-read.csv2("Microarrays/Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
write.table(disinf,"Submission/Supplementary_Dataset_2_disease_information_icds_categories_and_colors.txt",quote=F,sep="\t",row.names=F)




