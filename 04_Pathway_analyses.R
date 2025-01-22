## CODE FOR ANALYZING THE SIGNIFICANTLY ENRICHED PATHWAYS IDENTIFIED IN THE DIFFERENT COMPARISONS ##
## Developed by Jon Sanchez-Valle
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es

## Load libraries ##
library("ggplot2")
library("gplots")
library("forestplot")
library("dplyr")
if(!require(dendextend)) install.packages("dendextend")
library("dendextend")
library("pvclust")


## Create a ranked list of genes for each disease and sex to conduct enrichment analyses (GSEA) ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
comparison<-"Women"
icdcode<-"ICD10"
experiment<-"Microarrays"
rankgenes<-function(comparison,icdcode,experiment){
  if("Ranks"%in%list.files(paste(experiment,"/",icdcode,"/Comparisons/",comparison,sep=""))==FALSE){dir.create(paste(experiment,"/",icdcode,"/Comparisons/",comparison,"/Ranks",sep=""))}
  if("Pathways"%in%list.files(paste(experiment,"/",icdcode,"/Comparisons/",comparison,sep=""))==FALSE){dir.create(paste(experiment,"/",icdcode,"/Comparisons/",comparison,"/Pathways",sep=""))}
  files<-list.files(paste(experiment,"/",icdcode,"/Comparisons/",comparison,"/DifferentialExpressions/",sep=""))
  genes<-c()
  for(a in 1:length(files)){
    tabla<-read.csv2(paste(experiment,"/",icdcode,"/Comparisons/",comparison,"/DifferentialExpressions/",files[a],sep=""),stringsAsFactors = F,sep="\t")
    genes<-rbind(genes,
                 c(length(intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))>0))),
                   length(intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))<0))),
                   length(tabla[,1])))
    tabla<-cbind(rownames(tabla),tabla[,3])
    tabli<-tabla[order(as.numeric(as.character(tabla[,2])),decreasing = F),]
    write.table(tabli,paste(experiment,"/",icdcode,"/Comparisons/",comparison,"/Ranks/",gsub("-","__",gsub(".txt",".rnk",files[a])),sep=""),quote=F,sep="\t",row.names=F,col.names=F)
  }
  rownames(genes)<-gsub(".txt","",files)
  colnames(genes)<-c("up","down","total")
  write.table(genes,paste(experiment,"/",icdcode,"/Comparisons/",comparison,"/Outputs/Number_sDEGs_and_genes.txt",sep=""),quote=F,sep="\t")
  return(genes)
}

## Adjusting for sex ##
adj9net<-rankgenes("Adjusted","ICD9","Microarrays")
adj10net<-rankgenes("Adjusted","ICD10","Microarrays")

## Men and Women separately ##
women9net<-rankgenes("Women","ICD9","Microarrays")
men9net<-rankgenes("Men","ICD9","Microarrays")
women10net<-rankgenes("Women","ICD10","Microarrays")
men10net<-rankgenes("Men","ICD10","Microarrays")

## Women vs. Men in cases and controls ##
cases9net<-rankgenes("Cases","ICD9","Microarrays")
controls9net<-rankgenes("Controls","ICD9","Microarrays")
cases10net<-rankgenes("Cases","ICD10","Microarrays")
controls10net<-rankgenes("Controls","ICD10","Microarrays")


#### Run GSEA ####
#### This code is to generate the myrun file for running GSEA in the marenostrum ####
## module load java/16.0.1 (this is the version needed for running GSEA on the mare)

# therun9<-"GSEA_4.3.2/gsea-cli.sh GSEAPreranked -gmx GSEAfiles/gmt/c2.cp.kegg.v2023.1.Hs.symbols.gmt,GSEAfiles/gmt/c2.cp.reactome.v2023.1.Hs.symbols.gmt,GSEAfiles/gmt/c5.go.bp.v2023.1.Hs.symbols.gmt -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnd_seed timestamp -rnk ./Microarrays/ICD9/Comparisons/Controls/Ranks/996.rnk -scoring_scheme weighted -rpt_label 996 -chip ./GSEAfiles/chip/Human_Gene_Symbol_with_Remapping_MSigDB.v2023.1.Hs.chip -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -set_max 500 -set_min 15 -zip_report false -out ./Microarrays/ICD9/Comparisons/Controls/Pathways"
# therun10<-"GSEA_4.3.2/gsea-cli.sh GSEAPreranked -gmx GSEAfiles/gmt/c2.cp.kegg.v2023.1.Hs.symbols.gmt,GSEAfiles/gmt/c2.cp.reactome.v2023.1.Hs.symbols.gmt,GSEAfiles/gmt/c5.go.bp.v2023.1.Hs.symbols.gmt -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnd_seed timestamp -rnk ./Microarrays/ICD10/Comparisons/Controls/Ranks/996.rnk -scoring_scheme weighted -rpt_label 996 -chip ./GSEAfiles/chip/Human_Gene_Symbol_with_Remapping_MSigDB.v2023.1.Hs.chip -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -set_max 500 -set_min 15 -zip_report false -out ./Microarrays/ICD10/Comparisons/Controls/Pathways"
# comparaciones<-c("Adjusted","Cases","Controls","Women","Men")
# theruns<-c()
# for(z in 1:length(comparaciones)){
#   ## ICD9 ##
#   ficheros9<-list.files(paste("Microarrays/ICD9/Comparisons/",comparaciones[z],"/Ranks/",sep=""))
#   theruns9<-c() ; for(a in 1:length(ficheros9)){theruns9<-c(theruns9,gsub("Controls",comparaciones[z],gsub(996,gsub(".rnk","",ficheros9[a]),therun9)))}
#   ## ICD10 ##
#   ficheros10<-list.files(paste("Microarrays/ICD10/Comparisons/",comparaciones[z],"/Ranks/",sep=""))
#   theruns10<-c() ; for(a in 1:length(ficheros10)){theruns10<-c(theruns10,gsub("Controls",comparaciones[z],gsub(996,gsub(".rnk","",ficheros10[a]),therun10)))}
#   theruns<-c(theruns,theruns9,theruns10)
#   print(paste(round((z/length(comparaciones))*100,2),"%",sep=""))
# }
# write.table(theruns,"my_run_GSEA.txt",quote=F,sep="\t",row.names=F,col.names=F)

## Run it on the mare nostrum ##

#### Create a table for each condition ####
comparison<-"Women"
icdcode<-"ICD9"
experiment<-"Microarrays"
pathwaytables<-function(comparison,icdcode,experiment){
  if("PathwayTables"%in%list.files(paste(experiment,"/",icdcode,"/Comparisons/",comparison,sep=""))==FALSE){dir.create(paste(paste(experiment,"/",icdcode,"/Comparisons/",comparison,sep=""),"/PathwayTables",sep=""))}
  files<-list.files(paste(paste(experiment,"/",icdcode,"/Comparisons/",comparison,sep=""),"/Pathways",sep=""))
  pathnames<-c() ; pathlist<-c() ; allpathnames<-c() ; allpathlist<-list()
  for(a in 1:length(files)){
    # a<-1
    files2<-list.files(paste(paste(experiment,"/",icdcode,"/Comparisons/",comparison,sep=""),"/Pathways/",files[a],sep=""))
    negfiles<-files2[intersect(grep("gsea_report_for_na_neg",files2),grep(".tsv",files2))]
    neg<-read.csv2(paste(paste(experiment,"/",icdcode,"/Comparisons/",comparison,sep=""),"/Pathways/",files[a],"/",negfiles,sep=""),stringsAsFactors = F,sep="\t")
    posfiles<-files2[intersect(grep("gsea_report_for_na_pos",files2),grep(".tsv",files2))]
    pos<-read.csv2(paste(paste(experiment,"/",icdcode,"/Comparisons/",comparison,sep=""),"/Pathways/",files[a],"/",posfiles,sep=""),stringsAsFactors = F,sep="\t")
    tabla<-rbind(pos[,c(1,4,6,8,9)],neg[,c(1,4,6,8,9)])
    ## Remove those lines where we see a "---" in the NES ##
    if(length(grep("---",tabla$NES))>0){tabla<-tabla[-grep("---",tabla$NES),]}
    tabla<-tabla[order(abs(as.numeric(as.character(tabla$NES))),decreasing = T),]
    ## Save the table with pathways and their associated NES for each ICD code separately ##
    write.table(tabla,paste(paste(experiment,"/",icdcode,"/Comparisons/",comparison,sep=""),"/PathwayTables/",gsub(".GseaPre.+","",files[a]),".txt",sep=""),quote=F,sep="\t",row.names = F)
  }
  print("All the tables have been saved")
  files<-list.files(paste(paste(experiment,"/",icdcode,"/Comparisons/",comparison,sep=""),"/PathwayTables",sep=""))
  pathnames<-c() ; pathlist<-c() ; allpathnames<-c()
  for(a in 1:length(files)){
    # a<-1
    tabla<-read.csv2(paste(experiment,"/",icdcode,"/Comparisons/",comparison,"/PathwayTables/",files[a],sep=""),stringsAsFactors = F,sep="\t")
    tablan<-cbind(tabla$NAME,0,tabla$NES)
    if(length(which(tabla$FDR.q.val<=0.05))>0){
      if(length(intersect(which(tabla$FDR.q.val<=0.05),which(tabla$NES>0)))>0){
        tablan[intersect(which(tabla$FDR.q.val<=0.05),which(tabla$NES>0)),2]<-1
      }
      if(length(intersect(which(tabla$FDR.q.val<=0.05),which(tabla$NES<0)))>0){
        tablan[intersect(which(tabla$FDR.q.val<=0.05),which(tabla$NES<0)),2]<-(-1)
      }
    }
    tablan<-as.data.table(tablan)
    setkey(tablan,V1)
    pathlist[[gsub(".txt","",files[a])]]<-tablan
    ## Vector with the pathways significantly altered in at least one disease ##
    pathnames<-c(pathnames,tablan[which(tablan$V2!=0)]$V1)
  }
  ## Create the matrixes for the heatmap ##
  paths<-unique(pathnames)
  ## Binary and NES matrices ##
  binpathmat<-matrix(nrow=length(paths),ncol=length(pathlist),0)
  rownames(binpathmat)<-paths ; colnames(binpathmat)<-names(pathlist)
  nespathmat<-binpathmat
  for(a in 1:length(names(pathlist))){
    # a<-1
    binpathmat[intersect(pathlist[[a]]$V1,rownames(nespathmat)),names(pathlist)[a]]<-pathlist[[a]][intersect(pathlist[[a]]$V1,rownames(nespathmat))]$V2
    nespathmat[intersect(pathlist[[a]]$V1,rownames(nespathmat)),names(pathlist)[a]]<-pathlist[[a]][intersect(pathlist[[a]]$V1,rownames(nespathmat))]$V3
  }
  print("Matrices for the heatmap have been created")
  write.table(nespathmat,paste(experiment,"/",icdcode,"/Comparisons/",comparison,"/Outputs/Enrichments_nes.txt",sep=""),quote=F,sep="\t")
  write.table(binpathmat,paste(experiment,"/",icdcode,"/Comparisons/",comparison,"/Outputs/Enrichments_binarized.txt",sep=""),quote=F,sep="\t")
  ## Save the results ##
  resultados<-list("nes"=nespathmat,"binarized"=binpathmat)
  return(resultados)
}

## Adjusting for sex ##
adj9net<-pathwaytables("Adjusted","ICD9","Microarrays")
adj10net<-pathwaytables("Adjusted","ICD10","Microarrays")

## Men and Women separately ##
women9net<-pathwaytables("Women","ICD9","Microarrays")
men9net<-pathwaytables("Men","ICD9","Microarrays")
women10net<-pathwaytables("Women","ICD10","Microarrays")
men10net<-pathwaytables("Men","ICD10","Microarrays")

## Women vs. Men in cases and controls ##
cases9net<-pathwaytables("Cases","ICD9","Microarrays")
controls9net<-pathwaytables("Controls","ICD9","Microarrays")
cases10net<-pathwaytables("Cases","ICD10","Microarrays")
controls10net<-pathwaytables("Controls","ICD10","Microarrays")