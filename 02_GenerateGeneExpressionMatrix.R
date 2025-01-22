## ## CODE TO GENERATE GENE EXPRESSION MATRIX REMOVING LOWLY EXPRESSED PROBES ##
## Developed by Jon Sanchez-Valle
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es

## Install packages ##
# BiocManager::install(c("affy","frma","limma","hgu133plus2frmavecs","hgu133plus2.db","GEOmetadb","hgu133plus2CellScore","massiR","genefilter"))

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


args = commandArgs(trailingOnly=TRUE)

## @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
#### Gene Expression Matrix removing lowly expressed genes ####
## @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
if(as.numeric(args[1])==1){
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
  
  ## Identify the number of females and males by cases and controls: ##
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
  
  femid<-intersect(which(resultstab[,1]>=3),which(resultstab[,2]>=3))
  malid<-intersect(which(resultstab[,3]>=3),which(resultstab[,4]>=3))
  ## Function to generate gene expression matrices by ICD version and comparison removing lowly expressed genes:
  icdcode<-"ICD9"
  comparison<-"Adjusted"
  # comparison<-"All"
  removepoorqualityandlowexpression<-function(icdcode,comparison,resultstab){
    ## Create the needed directories 
    if(icdcode%in%list.files("Microarrays/")==FALSE){dir.create(paste("Microarrays/",icdcode,sep=""))}
    if(comparison=="All"){
      if("GeneExpressionMatrix"%in%list.files(paste("Microarrays/",icdcode,sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/GeneExpressionMatrix",sep=""))}
      if("Targets"%in%list.files(paste("Microarrays/",icdcode,sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/Targets",sep=""))}
    }
    if(comparison!="All"){
      if("Comparisons"%in%list.files(paste("Microarrays/",icdcode,sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/Comparisons",sep=""))}
      if(comparison%in%list.files(paste("Microarrays/",icdcode,"/Comparisons/",sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/Comparisons/",comparison,sep=""))}
      if("Targets"%in%list.files(paste("Microarrays/",icdcode,"/Comparisons/",comparison,sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/Comparisons/",comparison,"/Targets",sep=""))}
      if("Outputs"%in%list.files(paste("Microarrays/",icdcode,"/Comparisons/",comparison,sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/Comparisons/",comparison,"/Outputs",sep=""))}
      if("GeneExpressionMatrix"%in%list.files(paste("Microarrays/",icdcode,"/Comparisons/",comparison,sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/Comparisons/",comparison,"/GeneExpressionMatrix",sep=""))}
      if("DifferentialExpressions"%in%list.files(paste("Microarrays/",icdcode,"/Comparisons/",comparison,sep=""))==FALSE){dir.create(paste("Microarrays/",icdcode,"/Comparisons/",comparison,"/DifferentialExpressions",sep=""))}
    }
    
    ## Read the table with the ICD9 information ##
    tabla<-read.csv2("Microarrays/Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
    ## Select the desired version of the ICD codes ##
    if(icdcode=="ICD9"){
      if(length(which(tabla$ICD9=="-"))>0){tabla<-tabla[-which(tabla$ICD9=="-"),]}
      nametoicd<-tabla$ICD9 ; names(nametoicd)<-tabla$Disease
    }
    if(icdcode=="ICD10"){
      if(length(which(tabla$ICD10=="-"))>0){tabla<-tabla[-which(tabla$ICD10=="-"),]}
      nametoicd<-tabla$ICD10 ; names(nametoicd)<-tabla$Disease
    }
    ## Which comparison are we carrying out? ##
    ## "All" takes all the samples into consideration, no matter if there are more than 3 samples or not
    if(comparison=="All"){
      compadatasets3<-gsub(".txt","",fileswithknownsex)
    }
    ## The rest will request at least 3 cases and 3 controls of the analysed sex, or of at least one if we are adjusting for sex
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
    ## Remove those diseases or phenotypes with no ICD transformation
    if(length(which(is.na(compaicd)))>0){
      compadatasets<-compadatasets[-which(is.na(compaicd))]
      compadatasets2<-compadatasets3[-which(is.na(compaicd))]
    }
    print(paste("Starting ",icdcode," - ",comparison,":",sep=""))
    ## List of ICDs with enough samples
    compadiseases<-unique(gsub("_.+","",compadatasets))
    setdiff(compadiseases,gsub(".txt","",list.files(paste("Microarrays/",icdcode,"/GeneExpressionMatrix",sep=""))))
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
        if(comparison=="All"){
          targets<-rbind(targets,target[,c(1:3,5,4)])
        }
      }
      colnames(targets)<-c("Path","Sample","Category","Study","Sex")
      if(comparison!="All"){
        ## If the matrix was not generated in the general step avoid this step
        if(length(which(gsub(".txt","",list.files(paste("Microarrays/",icdcode,"/GeneExpressionMatrix",sep="")))==compadiseases[a]))>0){
          write.table(targets,paste("Microarrays/",icdcode,"/Comparisons/",comparison,"/Targets/",compadiseases[a],".txt",sep=""),quote=F,sep="\t",row.names=F)
        }
      }
      if(comparison=="All"){
        write.table(targets,paste("Microarrays/",icdcode,"/Targets/",compadiseases[a],".txt",sep=""),quote=F,sep="\t",row.names=F)
      }
      ## If there are duplicated samples remove them ##
      if(length(which(duplicated(targets$Sample)))>0){
        print(paste(a,"tiene rows duplicadas:",compadiseases[a]))
        targets<-targets[-which(duplicated(targets$Sample)),]
        ## Write targets file again ##
        if(comparison!="All"){
          ## If the matrix was not generated in the general step avoid this step
          if(length(which(gsub(".txt","",list.files(paste("Microarrays/",icdcode,"/GeneExpressionMatrix",sep="")))==compadiseases[a]))>0){
            write.table(targets,paste("Microarrays/",icdcode,"/Comparisons/",comparison,"/Targets/",compadiseases[a],".txt",sep=""),quote=F,sep="\t",row.names=F)
          }
        }
        if(comparison=="All"){
          write.table(targets,paste("Microarrays/",icdcode,"/Targets/",compadiseases[a],".txt",sep=""),quote=F,sep="\t",row.names=F)
        }
      }
      if(comparison=="All"){
        ## Evaluate if controlling for the experiment while having just one experiment gives different results
        affyBatch<-ReadAffy(filenames=targets$Path)
        ## Identify lowly expressed probes ##
        co <- mas5calls(affyBatch)
        co <- assayData(co)[["se.exprs"]] #extract detection p-values
        pvalue.detection.cutoff <- 0.05
        calls.matrix <- co < pvalue.detection.cutoff
        extranumbers<-function(vector){valor<-length(which(vector));return(valor)}
        totalexpressed<-apply(calls.matrix, 1, extranumbers)
        ## Identify the lowly expressed probes ##
        seleccionar<-which(totalexpressed>min(table(targets[,3:4])))
        if(length(seleccionar)>0){
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
          write.table(expression,paste("Microarrays/",icdcode,"/GeneExpressionMatrix/",compadiseases[a],".txt",sep=""),quote=F,sep="\t")
        }
        print(paste(round((a/length(compadiseases))*100,3),"%",sep=""))
      }
      if(comparison!="All"){
        ## If the matrix was not generated in the general step avoid this step
        if(length(which(gsub(".txt","",list.files(paste("Microarrays/",icdcode,"/GeneExpressionMatrix",sep="")))==compadiseases[a]))>0){
          expression<-read.csv(paste("Microarrays/",icdcode,"/GeneExpressionMatrix/",compadiseases[a],".txt",sep=""),stringsAsFactors = F,sep="\t",check.names = F)
          ## Select from the expression matrix the samples selected in the targets object at the beggining of the loop 
          expression<-expression[,targets$Sample]
          ## Check that the amount of samples selected at the beggining in the "targets" is the same as in the expression
          if(dim(expression)[2]==dim(targets)[1]){
            write.table(expression,paste("Microarrays/",icdcode,"/Comparisons/",comparison,"/GeneExpressionMatrix/",compadiseases[a],".txt",sep=""),quote=F,sep="\t")
          }
          if(dim(expression)[2]!=dim(targets)[1]){
            print(a)
          }
        }
        
        print(paste(round((a/length(compadiseases))*100,3),"%",sep=""))
      }
    }
    print("Finished")
  }
  
  ## ICD9 ##
  ## @@@@ ##
  ## General
  removepoorqualityandlowexpression("ICD9","All",resultstab)
  ## Comparisons
  removepoorqualityandlowexpression("ICD9","Women",resultstab)
  removepoorqualityandlowexpression("ICD9","Men",resultstab)
  removepoorqualityandlowexpression("ICD9","Cases",resultstab)
  removepoorqualityandlowexpression("ICD9","Controls",resultstab)
  removepoorqualityandlowexpression("ICD9","Adjusted",resultstab)
  
  ## ICD10 ##
  ## @@ @@ ##
  ## General
  removepoorqualityandlowexpression("ICD10","All",resultstab)
  ## Comparisons
  removepoorqualityandlowexpression("ICD10","Women",resultstab)
  removepoorqualityandlowexpression("ICD10","Men",resultstab)
  removepoorqualityandlowexpression("ICD10","Cases",resultstab)
  removepoorqualityandlowexpression("ICD10","Controls",resultstab)
  removepoorqualityandlowexpression("ICD10","Adjusted",resultstab)
}

## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
#### Number of samples per disease ####
## @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ ##
if(as.numeric(args[1])==4){
  files<-list.files("Microarrays/Data/FinalTargets/")
  femtab<-c() ; maltab<-c() ; stus<-c()
  for(a in 1:length(files)){
    # a<-1
    tabla<-read.csv2(paste("Microarrays/Data/FinalTargets/",files[a],sep=""),stringsAsFactors = F,sep="\t")
    if(length(which(names(table(tabla$Category))=="Case"))>0 && length(which(names(table(tabla$Category))=="Control"))>0){
      if(table(tabla$Category)[which(names(table(tabla$Category))=="Case")]>=3 && table(tabla$Category)[which(names(table(tabla$Category))=="Control")]>=3){
        stus<-c(stus,files[a])
      }
    }
    tablatable<-table(tabla[,3:4])
    sexos<-colnames(tablatable)
    if(length(which(sexos=="female"))>0){
      fem<-tablatable[,which(sexos=="female")]
      if(length(which(fem>3))==2){
        femtab<-rbind(femtab,c(gsub(".txt","",files[a]),as.numeric(fem[which(names(fem)=="Case")]),as.numeric(fem[which(names(fem)=="Control")])))
      }
    }
    if(length(which(sexos=="male"))>0){
      mal<-tablatable[,which(sexos=="male")]
      if(length(which(mal>3))==2){
        maltab<-rbind(maltab,c(gsub(".txt","",files[a]),as.numeric(mal[which(names(mal)=="Case")]),as.numeric(mal[which(names(mal)=="Control")])))
      }
    }
  }
  
  ## Number of studies with at least 3 cases and 3 controls:
  length(stus) ## 211 studies
  ## Number of diseases with at least 3 cases and 3 controls per study:
  length(unique(gsub("_.+","",stus))) # 128 diseases
  
  ## Number of diseases women ##
  length(unique(gsub("_.+","",femtab[,1]))) # 76 diseases 
  ## Number of cases in women ##
  sum(as.numeric(femtab[,2])) # 2,301
  ## Number of controls in women ##
  sum(as.numeric(femtab[,3])) # 1,355
  
  ## Number of diseases men ##
  length(unique(gsub("_.+","",maltab[,1]))) # 84 diseases
  ## Number of cases in men ##
  sum(as.numeric(maltab[,2])) # 3,092
  ## Number of controls in men ##
  sum(as.numeric(maltab[,3])) # 1,928
  
  ## ICD10 ##
  ## @@ @@ ##
  ## Women ##
  wtar<-list.files("Microarrays/ICD10/Comparisons/Women/Targets/")
  wsampls10<-c()
  for(a in 1:length(wtar)){
    tabla<-read.csv2(paste("Microarrays/ICD10/Comparisons/Women/Targets/",wtar[a],sep=""),stringsAsFactors = F,sep="\t")
    wsampls10<-rbind(wsampls10,c(gsub(".txt","",wtar[a]),length(which(tabla$Category=="Case")),length(which(tabla$Category=="Control"))))
  }
  ## Men ##
  mtar<-list.files("Microarrays/ICD10/Comparisons/Men/Targets/")
  msampls10<-c()
  for(a in 1:length(mtar)){
    tabla<-read.csv2(paste("Microarrays/ICD10/Comparisons/Men/Targets/",mtar[a],sep=""),stringsAsFactors = F,sep="\t")
    msampls10<-rbind(msampls10,c(gsub(".txt","",mtar[a]),length(which(tabla$Category=="Case")),length(which(tabla$Category=="Control"))))
  }
  ## Women cases ##
  sum(as.numeric(wsampls10[,2])) # 2,465
  ## Women controls ##
  sum(as.numeric(wsampls10[,3])) # 1,370
  
  ## Men cases ##
  sum(as.numeric(msampls10[,2])) # 3,200
  ## Men controls ##
  sum(as.numeric(msampls10[,3])) # 1,871
  
  ## ICD9 ##
  ## @@@@ ##
  ## Women ##
  wtar<-list.files("Microarrays/ICD9/Comparisons/Women/Targets/")
  wsampls9<-c()
  for(a in 1:length(wtar)){
    tabla<-read.csv2(paste("Microarrays/ICD9/Comparisons/Women/Targets/",wtar[a],sep=""),stringsAsFactors = F,sep="\t")
    wsampls9<-rbind(wsampls9,c(gsub(".txt","",wtar[a]),length(which(tabla$Category=="Case")),length(which(tabla$Category=="Control"))))
  }
  ## Men ##
  mtar<-list.files("Microarrays/ICD9/Comparisons/Men/Targets/")
  msampls9<-c()
  for(a in 1:length(mtar)){
    tabla<-read.csv2(paste("Microarrays/ICD9/Comparisons/Men/Targets/",mtar[a],sep=""),stringsAsFactors = F,sep="\t")
    msampls9<-rbind(msampls9,c(gsub(".txt","",mtar[a]),length(which(tabla$Category=="Case")),length(which(tabla$Category=="Control"))))
  }
  
  ## Women cases ##
  sum(as.numeric(wsampls9[,2])) # 2,465
  ## Women controls ##
  sum(as.numeric(wsampls9[,3])) # 1,339
  
  ## Men cases ##
  sum(as.numeric(msampls9[,2])) # 3,190
  ## Men controls ##
  sum(as.numeric(msampls9[,3])) # 1,810
}







