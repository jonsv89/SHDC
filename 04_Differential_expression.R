## CODE FOR THE DIFFERENTIAL EXPRESSION ANALYSIS FOR EACH DISEASE BY MEN, WOMEN, AND ADJUSTING BY SEX ##
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
library("statmod")


icdcode<-"ICD9"
comparison<-"Cases"
## Function: Conduct differential gene expression analyses in each comparison ##
differentialexpressionanalysis<-function(icdcode,comparison){
  ## We are going to process all the diseases in the same loop, no matter the sex ##
  files<-list.files(paste("Microarrays",icdcode,comparison,"GeneExpressionMatrix/",sep="/"))
  summaryl<-list()
  for(a in 1:length(files)){
    targets<-read.csv2(paste("Microarrays/",icdcode,"/",comparison,"/Targets/",files[a],sep=""),stringsAsFactors = F,sep="\t")
    expression2<-read.csv2(paste("Microarrays/",icdcode,"/",comparison,"/GeneExpressionMatrix/",files[a],sep=""),stringsAsFactors = F,sep="\t")
    expression<-c() ; for(b in 1:length(expression2[1,])){expression<-cbind(expression,as.numeric(expression2[,b]))}
    colnames(expression)<-colnames(expression2) ; rownames(expression)<-rownames(expression2)
    ## Making Factors for the Category and Study ##
    study <- factor(targets$Study)
    if(comparison=="Women" || comparison=="Men"){category <- factor(targets$Category, levels=c("Control","Case"))}
    if(comparison=="Cases" || comparison=="Controls"){category <- factor(targets$Sex, levels=c("male","female"))}
    if(comparison=="Adjusted"){
      category <- factor(targets$Category, levels=c("Control","Case"))
      sex <- factor(targets$Sex)
    }
    if(comparison!="Adjusted"){
      ## We only have one study for the disease ##
      if(length(unique(targets$Study))==1){
        design <- model.matrix(~category)
        fit <- lmFit(expression, design)
        fit2 <- eBayes(fit, trend=TRUE, robust=TRUE)
        results <- decideTests(fit2)
        summaryl[[gsub(".txt","",files[a])]]<-summary(results)
        if(comparison=="Women" || comparison=="Men"){thetoptable<-topTable(fit2, coef="categoryCase",number=length(fit$Amean))}
        if(comparison=="Cases" || comparison=="Controls"){thetoptable<-topTable(fit2, coef="categoryfemale",number=length(fit$Amean))}
        ## save all patient vs. all control differential expression results
        write.table(thetoptable,paste("Microarrays/",icdcode,"/",comparison,"/DifferentialExpressions/",files[a],sep=""),sep="\t",quote=F)
      }
      ## We have more than one study for the disease ##
      ## Compare cases and controls adjusting for differences between studies ##
      if(length(unique(targets$Study))>1){
        design <- model.matrix(~category+study)
        fit <- lmFit(expression, design)
        fit2 <- eBayes(fit, trend=TRUE, robust=TRUE)
        results <- decideTests(fit2)
        summaryl[[gsub(".txt","",files[a])]]<-summary(results)
        if(comparison=="Women" || comparison=="Men"){thetoptable<-topTable(fit2, coef="categoryCase",number=length(fit$Amean))}
        if(comparison=="Cases" || comparison=="Controls"){thetoptable<-topTable(fit2, coef="categoryfemale",number=length(fit$Amean))}
        ## save all patient vs. all control differential expression results
        write.table(thetoptable,paste("Microarrays/",icdcode,"/",comparison,"/DifferentialExpressions/",files[a],sep=""),sep="\t",quote=F)
      }
    }
    if(comparison=="Adjusted"){
      ## We only have one study for the disease ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
      if(length(unique(targets$Study))==1){
        ## We only have one sex ##
        if(length(unique(targets$Sex))==1){
          design <- model.matrix(~category)
          fit <- lmFit(expression, design)
          fit2 <- eBayes(fit, trend=TRUE, robust=TRUE)
          results <- decideTests(fit2)
          summaryl[[gsub(".txt","",files[a])]]<-summary(results)
          thetoptable<-topTable(fit2, coef="categoryCase",number=length(fit$Amean))
          ## save all patient vs. all control differential expression results
          write.table(thetoptable,paste("Microarrays/",icdcode,"/",comparison,"/DifferentialExpressions/",files[a],sep=""),sep="\t",quote=F)
        }
        ## We have both sexs ##
        if(length(unique(targets$Sex))>1){
          design <- model.matrix(~category+sex)
          fit <- lmFit(expression, design)
          fit2 <- eBayes(fit, trend=TRUE, robust=TRUE)
          results <- decideTests(fit2)
          summaryl[[gsub(".txt","",files[a])]]<-summary(results)
          thetoptable<-topTable(fit2, coef="categoryCase",number=length(fit$Amean))
          ## save all patient vs. all control differential expression results
          write.table(thetoptable,paste("Microarrays/",icdcode,"/",comparison,"/DifferentialExpressions/",files[a],sep=""),sep="\t",quote=F)
        }
      }
      
      ## We have more than one study for the disease ##
      ## @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ ##
      ## Compare cases and controls adjusting for differences between studies ##
      if(length(unique(targets$Study))>1){
        if(length(unique(targets$Sex))==1){
          design <- model.matrix(~category+study)
          fit <- lmFit(expression, design)
          fit2 <- eBayes(fit, trend=TRUE, robust=TRUE)
          results <- decideTests(fit2)
          summaryl[[gsub(".txt","",files[a])]]<-summary(results)
          thetoptable<-topTable(fit2, coef="categoryCase",number=length(fit$Amean))
          ## save all patient vs. all control differential expression results
          write.table(thetoptable,paste("Microarrays/",icdcode,"/",comparison,"/DifferentialExpressions/",files[a],sep=""),sep="\t",quote=F)
        }
        if(length(unique(targets$Sex))>1){
          reducido<-targets[-which(duplicated(targets[,c(4,5)])),c(4,5)]
          ## If each study has a different sex, control only for the study
          if(length(which(table(reducido$Study)>1))==0){
            design <- model.matrix(~category+study)
            fit <- lmFit(expression, design)
            fit2 <- eBayes(fit, trend=TRUE, robust=TRUE)
            results <- decideTests(fit2)
            summaryl[[gsub(".txt","",files[a])]]<-summary(results)
            thetoptable<-topTable(fit2, coef="categoryCase",number=length(fit$Amean))
            ## save all patient vs. all control differential expression results
            write.table(thetoptable,paste("Microarrays/",icdcode,"/",comparison,"/DifferentialExpressions/",files[a],sep=""),sep="\t",quote=F)
          }
          ## If not, control for both
          if(length(which(table(reducido$Study)>1))>0){
            design <- model.matrix(~category+study+sex)
            fit <- lmFit(expression, design)
            fit2 <- eBayes(fit, trend=TRUE, robust=TRUE)
            results <- decideTests(fit2)
            summaryl[[gsub(".txt","",files[a])]]<-summary(results)
            thetoptable<-topTable(fit2, coef="categoryCase",number=length(fit$Amean))
            ## save all patient vs. all control differential expression results
            write.table(thetoptable,paste("Microarrays/",icdcode,"/",comparison,"/DifferentialExpressions/",files[a],sep=""),sep="\t",quote=F)
          }
        }
      }
    }
    print(a/length(files))
  }
  
  return(summaryl)
}
## Function: Obtain the number of differentially expressed genes in each comparison ##
numberofsdegs<-function(icdcode,comparison){
  files<-list.files(paste("Microarrays/",icdcode,"/",comparison,"/DifferentialExpressions",sep=""))
  tabla<-c()
  for(a in 1:length(files)){
    degtable<-read.csv2(paste("Microarrays/",icdcode,"/",comparison,"/DifferentialExpressions/",files[a],sep=""),stringsAsFactors = F,sep="\t")
    tabla<-rbind(tabla,c(length(which(as.numeric(degtable$adj.P.Val)<=0.05)),length(intersect(which(as.numeric(degtable$adj.P.Val)<=0.05),which(as.numeric(degtable$logFC)<0))),
                         length(intersect(which(as.numeric(degtable$adj.P.Val)<=0.05),which(as.numeric(degtable$logFC)>0)))))
  }
  rownames(tabla)<-gsub(".txt","",files)
  colnames(tabla)<-c("Total","Down","Up")
  write.table(tabla,paste("Microarrays/",icdcode,"/",comparison,"/Number_DEGs.txt",sep=""),quote=F,sep="\t")
  return(tabla)
}

## Differential expression ##
icd9women<-differentialexpressionanalysis("ICD9","Women")
icd9men<-differentialexpressionanalysis("ICD9","Men")
icd9cases<-differentialexpressionanalysis("ICD9","Cases")
icd9controls<-differentialexpressionanalysis("ICD9","Controls")
icd9adjusted<-differentialexpressionanalysis("ICD9","Adjusted")
icd10women<-differentialexpressionanalysis("ICD10","Women")
icd10men<-differentialexpressionanalysis("ICD10","Men")
icd10cases<-differentialexpressionanalysis("ICD10","Cases")
icd10controls<-differentialexpressionanalysis("ICD10","Controls")
icd10adjusted<-differentialexpressionanalysis("ICD10","Adjusted")

## Number of differentially expressed genes ##
nsdegsicd9women<-numberofsdegs("ICD9","Women")
nsdegsicd9men<-numberofsdegs("ICD9","Men")
nsdegsicd9cases<-numberofsdegs("ICD9","Cases")
nsdegsicd9controls<-numberofsdegs("ICD9","Controls")
nsdegsicd9adjusted<-numberofsdegs("ICD9","Adjusted")
nsdegsicd10women<-numberofsdegs("ICD10","Women")
nsdegsicd10men<-numberofsdegs("ICD10","Men")
nsdegsicd10cases<-numberofsdegs("ICD10","Cases")
nsdegsicd10controls<-numberofsdegs("ICD10","Controls")
nsdegsicd10adjusted<-numberofsdegs("ICD10","Adjusted")








