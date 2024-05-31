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

#### Microarrays ####
## @@ @@ @ @ @@ @@ ##

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

#### Prepare the plots and dendrograms ####
## Colors for the ICD9 diseases ##
info<-read.csv2("Microarrays/Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
sindupcol<-info[-(which(duplicated(info[,c(3,7)]))),c(3,7)]
colorcat9<-sindupcol[,2] ; names(colorcat9)<-sindupcol[,1]
## ICD10
sindupcol<-info[-(which(duplicated(info[,c(8,12)]))),c(8,12)]
colorcat10<-sindupcol[,2] ; names(colorcat10)<-sindupcol[,1]
## Colors for the pathways ##
patcatcol<-read.csv2("Microarrays/Data/final_colors_reactome_category.txt",stringsAsFactors = F,sep="\t")
reaccolors<-patcatcol[,2] ; names(reaccolors)<-patcatcol[,1]
## Identify the reactome parents ##
reactomeparents<-read.csv2("Microarrays/Data/GSEAReactome_ParentNames.txt",stringsAsFactors = F,sep="\t")
reactomeparents<-reactomeparents[order(reactomeparents[,3]),]
factreactomeparents<-reactomeparents[,3] ; names(factreactomeparents)<-reactomeparents[,1]

## Function to extract generate the table selecting only the pathways of interest, in this case REACTOME pathways ##
filtro<-"REACTOME" ; comparison<-"Women" ; icdcode<-"ICD10" ; experiment<-"Microarrays" ; nombre<-"Enrichments_allnes.txt"
prepare_tables<-function(comparison,icdcode,experiment,nombre,reactomeparents,filtro){
  tabla<-read.csv2(paste(experiment,"/",icdcode,"/Comparisons/",comparison,"/Outputs/",nombre,sep=""),stringsAsFactors = F,sep="\t")
  if(icdcode=="ICD9"){colnames(tabla)<-gsub("X","",colnames(tabla))}
  innumbers<-c();for(a in 1:length(tabla[1,])){innumbers<-cbind(innumbers,as.numeric(tabla[,a]))} 
  colnames(innumbers)<-colnames(tabla) ; rownames(innumbers)<-rownames(tabla)
  if(filtro=="REACTOME"){
    innumbers<-innumbers[grep("^REACTOME",rownames(innumbers)),]
    indices<-c();for(a in 1:length(rownames(innumbers))){
      indices<-c(indices,which(reactomeparents[,1]==rownames(innumbers)[a])[1])
    }
    innumbers2<-innumbers[order(indices),]
    categorias<-as.character(factreactomeparents[rownames(innumbers2)])
    ## In case we do not have the category for a specific pathway ## 
    if(length(which(is.na(categorias)))>0){
      innumbers2<-innumbers2[-which(is.na(categorias)),]
      categorias<-categorias[-which(is.na(categorias))]
    }
    unicats<-unique(categorias)
    separador<-c();for(a in 1:length(unicats)){separador<-c(separador,max(which(categorias==unicats[a])))}
  }
  listas<-list("table"=innumbers2,"separator"=separador)
  return(listas)
}

#### Compare dendrograms ####
icdcode<-"ICD10"
experiment<-"Microarrays"
nombre<-"women_vs_men_clustering"
comparedendograms<-function(icdcode,experiment,nombre,firstcluster,secondcluster){
  comfirstcluster<-firstcluster[,intersect(colnames(firstcluster),colnames(secondcluster))]
  comsecondcluster<-secondcluster[,intersect(colnames(firstcluster),colnames(secondcluster))]
  ## Are there NAs? in this case remove them ##
  ## In the first comparison ##
  haynasfirst<-c()
  for(a in 1:length(comfirstcluster[,1])){
    haynasfirst<-c(haynasfirst,length(which(is.na(comfirstcluster[a,]))))
  }
  comfirstcluster_2<-comfirstcluster[which(haynasfirst==0),]
  ## In the second comparison ##
  haynassecond<-c()
  for(a in 1:length(comsecondcluster[,1])){
    haynassecond<-c(haynassecond,length(which(is.na(comsecondcluster[a,]))))
  }
  comsecondcluster_2<-comsecondcluster[which(haynassecond==0),]
  
  ## Calculate pvalues by bootstrapping ##
  set.seed(1)
  ## Males ##
  clustfirst <- pvclust(comfirstcluster_2, method.dist="euclidean",
                    method.hclust="ward.D2", nboot = 1000)
  ## Females ##
  clustsecond <- pvclust(comsecondcluster_2, method.dist="euclidean",
                    method.hclust="ward.D2", nboot = 1000)
  
  pdf(file = paste("Microarrays/Plots/ClustersWithPvalues_",icdcode,".pdf",sep=""))
    plot(clustfirst, hang = -1, cex = 0.5,main="Cluster dendrogram with p-values (%)\nWomen")
    pvrect(clustfirst)
    plot(clustsecond, hang = -1, cex = 0.5,main="Cluster dendrogram with p-values (%)\nMen")
    pvrect(clustsecond)
  dev.off()
  
  ## Distances and clusters ##
  ## Males ##
  firstclusterdist<-dist(t(comfirstcluster_2),method = "euclidean")
  firstclusterclust<-hclust(firstclusterdist,method="ward.D2")
  firstclusterdend<-as.dendrogram(firstclusterclust)
  ## Females ##
  secondclusterdist<-dist(t(comsecondcluster_2),method = "euclidean")
  secondclusterclust<-hclust(secondclusterdist,method="ward.D2")
  secondclusterdend<-as.dendrogram(secondclusterclust)
  
  ## Correlation between distances ##
  correlation<-cor.test(firstclusterdist,secondclusterdist)
  
  ## Differences between dendrograms ##
  dend_diff(firstclusterdend,secondclusterdend)
  pdf(file=paste(experiment,"/Plots/",nombre,"_",icdcode,".pdf",sep=""))
    tanglegram(dend_diff(firstclusterdend,secondclusterdend))
  dev.off()
  ## The lower the better ##
  similarity<-dend_diff(firstclusterdend,secondclusterdend) %>% entanglement
  resultado<-list("cor"=correlation,"sim"=similarity,"womenclustlab"=labels(firstclusterdend),"menclustlab"=labels(secondclusterdend),
                  "clustwompval"=clustfirst,"clustmenpval"=clustsecond)
  return(resultado)
}

#### Plot the heatmaps to be combined with the tanglegram above ####
## ICD9 ##
## @@@@ ##
lallnes9<-prepare_tables("Women","ICD9","Microarrays","Enrichments_nes.txt",reactomeparents,"REACTOME")
firstcluster<-lallnes9$table
lallnes9<-prepare_tables("Men","ICD9","Microarrays","Enrichments_nes.txt",reactomeparents,"REACTOME")
secondcluster<-lallnes9$table
icd9comp<-comparedendograms("ICD9","Microarrays","women_vs_men_clustering",firstcluster,secondcluster)
## Get the binarized tables in the needed shape ##
wbin9<-prepare_tables("Women","ICD9","Microarrays","Enrichments_binarized.txt",reactomeparents,"REACTOME")
mbin9<-prepare_tables("Men","ICD9","Microarrays","Enrichments_binarized.txt",reactomeparents,"REACTOME")
## Put the diseases ordered properly to add the to the tanglegram ##
w9<-wbin9$table[,icd9comp$womenclustlab[length(icd9comp$womenclustlab):1]]
m9<-mbin9$table[,icd9comp$menclustlab[length(icd9comp$menclustlab):1]]
## Transpose the matrixes ##
tw9<-t(w9)
tm9<-t(m9)
## Get the colors of the pathways and the diseases ##
wcolorcategory9<-as.character(colorcat9[rownames(tw9)])
wcolorpathways9<-as.character(reaccolors[as.character(factreactomeparents[colnames(tw9)])])
mcolorcategory9<-as.character(colorcat9[rownames(tm9)])
mcolorpathways9<-as.character(reaccolors[as.character(factreactomeparents[colnames(tm9)])])
## Plot the heatmaps ##
pdf(file = "Microarrays/Plots/ICD9_heatmaps_for_tanglegram_addition_shareddiseases.pdf")
  heatmap.2(tw9,key=F,
            colsep = wbin9$separator,sepcolor = "gray",
            col=colorpanel(100, "#405191","#FFFFFF","#C81E17"),
            dendrogram = "none",Colv=FALSE,Rowv=FALSE,labCol = "",cexRow = 0.7,
            scale="none", margins=c(5,5), trace="none",colRow = wcolorcategory9,
            ColSideColors = wcolorpathways9)
  heatmap.2(tm9,key=F,
            colsep = mbin9$separator,sepcolor = "gray",
            col=colorpanel(100, "#405191","#FFFFFF","#C81E17"),
            dendrogram = "none",Colv=FALSE,Rowv=FALSE,labCol = "",cexRow = 0.7,
            scale="none", margins=c(5,5), trace="none",colRow = mcolorcategory9,
            ColSideColors = mcolorpathways9)
dev.off()

## Create the legend for the heatmaps ##
## @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ ##
patcatcol<-read.csv2("Microarrays/Data/final_colors_reactome_category.txt",stringsAsFactors = F,sep="\t")
reaccolors<-patcatcol[,2] ; names(reaccolors)<-patcatcol[,1]
matriz<-matrix(ncol=2,nrow=length(reaccolors),0)
rownames(matriz)<-names(reaccolors)
## The plot ##
pdf(file="Microarrays/Plots/ColorPathways.pdf")
  heatmap.2(matriz,RowSideColors = as.character(reaccolors),key=F,
            scale="none",margins=c(2,20), trace="none",dendrogram = "none",Rowv=FALSE)
dev.off()
## Colors for the ICD9 diseases ##
info<-read.csv2("Microarrays/Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
sindupcol<-info[-(which(duplicated(info[,c(6,7)]))),c(6,7)]
colorcat9<-sindupcol[,2] ; names(colorcat9)<-sindupcol[,1]
matriz<-matrix(ncol=2,nrow=length(colorcat9),0)
rownames(matriz)<-names(colorcat9)
## The plot ##
pdf(file="Microarrays/Plots/ColorCategoriesICD9.pdf",width = 20,height = 10)
  heatmap.2(matriz,RowSideColors = as.character(colorcat9),key=F,
            scale="none",margins=c(2,100), trace="none",dendrogram = "none",Rowv=FALSE,cexRow = 3)
dev.off()


## ICD10 ##
## @@ @@ ##
lallnes10<-prepare_tables("Women","ICD10","Microarrays","Enrichments_nes.txt",reactomeparents,"REACTOME")
firstcluster<-lallnes10$table
lallnes10<-prepare_tables("Men","ICD10","Microarrays","Enrichments_nes.txt",reactomeparents,"REACTOME")
secondcluster<-lallnes10$table
icd10comp<-comparedendograms("ICD10","Microarrays","women_vs_men_clustering",firstcluster,secondcluster)
## Get the binarized tables in the needed shape ##
wbin10<-prepare_tables("Women","ICD10","Microarrays","Enrichments_binarized.txt",reactomeparents,"REACTOME")
mbin10<-prepare_tables("Men","ICD10","Microarrays","Enrichments_binarized.txt",reactomeparents,"REACTOME")
## Put the diseases ordered properly to add the to the tanglegram ##
w10<-wbin10$table[,icd10comp$womenclustlab[length(icd10comp$womenclustlab):1]]
m10<-mbin10$table[,icd10comp$menclustlab[length(icd10comp$menclustlab):1]]
## Transpose the matrixes ##
tw10<-t(w10)
tm10<-t(m10)
## Get the colors of the pathways and the diseases ##
wcolorcategory10<-as.character(colorcat10[rownames(tw10)])
wcolorpathways10<-as.character(reaccolors[as.character(factreactomeparents[colnames(tw10)])])
mcolorcategory10<-as.character(colorcat10[rownames(tm10)])
mcolorpathways10<-as.character(reaccolors[as.character(factreactomeparents[colnames(tm10)])])
## Plot the heatmaps ##
pdf(file = "Microarrays/Plots/ICD10_heatmaps_for_tanglegram_addition_shareddiseases.pdf")
  heatmap.2(tw10,key=F,
            colsep = wbin10$separator,sepcolor = "gray",
            col=colorpanel(100, "#405191","#FFFFFF","#C81E17"),
            dendrogram = "none",Colv=FALSE,Rowv=FALSE,labCol = "",cexRow = 0.7,
            scale="none", margins=c(5,5), trace="none",colRow = wcolorcategory10,
            ColSideColors = wcolorpathways10)
  heatmap.2(tm10,key=F,
            colsep = mbin10$separator,sepcolor = "gray",
            col=colorpanel(100, "#405191","#FFFFFF","#C81E17"),
            dendrogram = "none",Colv=FALSE,Rowv=FALSE,labCol = "",cexRow = 0.7,
            scale="none", margins=c(5,5), trace="none",colRow = mcolorcategory10,
            ColSideColors = mcolorpathways10)
dev.off()
## Colors for the ICD10 diseases ##
info<-read.csv2("Microarrays/Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
sindupcol<-info[-(which(duplicated(info[,c(11,12)]))),c(11,12)]
colorcat10<-sindupcol[,2] ; names(colorcat10)<-sindupcol[,1]
matriz<-matrix(ncol=2,nrow=length(colorcat10),0)
rownames(matriz)<-names(colorcat10)
## The plot ##
pdf(file="Microarrays/Plots/ColorCategoriesICD10.pdf",width = 20,height = 10)
  heatmap.2(matriz,RowSideColors = as.character(colorcat10),key=F,
            scale="none",margins=c(2,100), trace="none",dendrogram = "none",Rowv=FALSE,cexRow = 2.5)
dev.off()


