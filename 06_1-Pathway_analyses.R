## CODE FOR ANALYZING THE SIGNIFICANTLY ENRICHED PATHWAYS IDENTIFIED IN THE DIFFERENT COMPARISONS ##
## Developed by Jon Sanchez-Valle
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es

## Load libraries ##
library("dendextend")
library("ggplot2")
library("gplots")
library("forestplot")
library("dplyr")


## Create a ranked list of genes for each disease and sex to conduct enrichment analyses (GSEA) ##
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
comparison<-"Women"
icdcode<-"ICD9"
experiment<-"Microarrays"
rankgenes<-function(comparison,icdcode,experiment){
  if("Ranks"%in%list.files(paste(experiment,"/",icdcode,"/",comparison,sep=""))==FALSE){dir.create(paste(experiment,"/",icdcode,"/",comparison,"/Ranks",sep=""))}
  if("Pathways"%in%list.files(paste(experiment,"/",icdcode,"/",comparison,sep=""))==FALSE){dir.create(paste(experiment,"/",icdcode,"/",comparison,"/Pathways",sep=""))}
  files<-list.files(paste(experiment,"/",icdcode,"/",comparison,"/DifferentialExpressions/",sep=""))
  genes<-c()
  for(a in 1:length(files)){
    tabla<-read.csv2(paste(experiment,"/",icdcode,"/",comparison,"/DifferentialExpressions/",files[a],sep=""),stringsAsFactors = F,sep="\t")
    genes<-rbind(genes,
                 c(length(intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))>0))),
                   length(intersect(which(as.numeric(as.character(tabla$adj.P.Val))<=0.05),which(as.numeric(as.character(tabla$logFC))<0))),
                   length(tabla[,1])))
    tabla<-cbind(rownames(tabla),tabla[,3])
    tabli<-tabla[order(as.numeric(as.character(tabla[,2])),decreasing = F),]
    write.table(tabli,paste(experiment,"/",icdcode,"/",comparison,"/Ranks/",gsub("-","__",gsub(".txt",".rnk",files[a])),sep=""),quote=F,sep="\t",row.names=F,col.names=F)
  }
  rownames(genes)<-gsub(".txt","",files)
  colnames(genes)<-c("up","down","total")
  write.table(genes,paste(experiment,"/",icdcode,"/",comparison,"/Number_sDEGs_and_genes.txt",sep=""),quote=F,sep="\t")
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
# ## Run, on the SexAssociatedComorbidities directory, the following scripts for the enrichment analysis ##
# ## Adjusted ##
# for a in `ls SexAssociatedComorbidities/Microarrays/ICD10/Adjusted/Ranks/`; do java -Xmx2000m -cp ./gsea.jar xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.reactome.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v7.0.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ./SexAssociatedComorbidities/Microarrays/ICD10/Adjusted/Ranks/$a -scoring_scheme weighted -rpt_label ${a/%????} -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./SexAssociatedComorbidities/Microarrays/ICD10/Adjusted/Pathways -gui false; done
# for a in `ls SexAssociatedComorbidities/Microarrays/ICD10/Adjusted/Ranks/`; do java -Xmx2000m -cp ./gsea.jar xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.reactome.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v7.0.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ./SexAssociatedComorbidities/Microarrays/ICD10/Adjusted/Ranks/$a -scoring_scheme weighted -rpt_label ${a/%????} -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./SexAssociatedComorbidities/Microarrays/ICD10/Adjusted/Pathways -gui false; done
# ## Women ##
# for a in `ls SexAssociatedComorbidities/Microarrays/ICD10/Women/Ranks/`; do java -Xmx2000m -cp ./gsea.jar xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.reactome.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v7.0.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ./SexAssociatedComorbidities/Microarrays/ICD10/Women/Ranks/$a -scoring_scheme weighted -rpt_label ${a/%????} -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./SexAssociatedComorbidities/Microarrays/ICD10/Women/Pathways -gui false; done
# for a in `ls SexAssociatedComorbidities/Microarrays/ICD10/Women/Ranks/`; do java -Xmx2000m -cp ./gsea.jar xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.reactome.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v7.0.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ./SexAssociatedComorbidities/Microarrays/ICD10/Women/Ranks/$a -scoring_scheme weighted -rpt_label ${a/%????} -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./SexAssociatedComorbidities/Microarrays/ICD10/Women/Pathways -gui false; done
# ## Men ##
# for a in `ls SexAssociatedComorbidities/Microarrays/ICD10/Men/Ranks/`; do java -Xmx2000m -cp ./gsea.jar xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.reactome.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v7.0.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ./SexAssociatedComorbidities/Microarrays/ICD10/Men/Ranks/$a -scoring_scheme weighted -rpt_label ${a/%????} -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./SexAssociatedComorbidities/Microarrays/ICD10/Men/Pathways -gui false; done
# for a in `ls SexAssociatedComorbidities/Microarrays/ICD10/Men/Ranks/`; do java -Xmx2000m -cp ./gsea.jar xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.reactome.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v7.0.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ./SexAssociatedComorbidities/Microarrays/ICD10/Men/Ranks/$a -scoring_scheme weighted -rpt_label ${a/%????} -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./SexAssociatedComorbidities/Microarrays/ICD10/Men/Pathways -gui false; done
# ## Cases ##
# for a in `ls SexAssociatedComorbidities/Microarrays/ICD10/Cases/Ranks/`; do java -Xmx2000m -cp ./gsea.jar xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.reactome.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v7.0.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ./SexAssociatedComorbidities/Microarrays/ICD10/Cases/Ranks/$a -scoring_scheme weighted -rpt_label ${a/%????} -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./SexAssociatedComorbidities/Microarrays/ICD10/Cases/Pathways -gui false; done
# for a in `ls SexAssociatedComorbidities/Microarrays/ICD10/Cases/Ranks/`; do java -Xmx2000m -cp ./gsea.jar xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.reactome.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v7.0.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ./SexAssociatedComorbidities/Microarrays/ICD10/Cases/Ranks/$a -scoring_scheme weighted -rpt_label ${a/%????} -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./SexAssociatedComorbidities/Microarrays/ICD10/Cases/Pathways -gui false; done
# ## Controls ##
# for a in `ls SexAssociatedComorbidities/Microarrays/ICD10/Controls/Ranks/`; do java -Xmx2000m -cp ./gsea.jar xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.reactome.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v7.0.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ./SexAssociatedComorbidities/Microarrays/ICD10/Controls/Ranks/$a -scoring_scheme weighted -rpt_label ${a/%????} -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./SexAssociatedComorbidities/Microarrays/ICD10/Controls/Pathways -gui false; done
# for a in `ls SexAssociatedComorbidities/Microarrays/ICD10/Controls/Ranks/`; do java -Xmx2000m -cp ./gsea.jar xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.reactome.v7.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.bp.v7.0.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ./SexAssociatedComorbidities/Microarrays/ICD10/Controls/Ranks/$a -scoring_scheme weighted -rpt_label ${a/%????} -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./SexAssociatedComorbidities/Microarrays/ICD10/Controls/Pathways -gui false; done


#### Create a table for each condition ####
comparison<-"Women"
icdcode<-"ICD9"
experiment<-"Microarrays"
pathwaytables<-function(comparison,icdcode,experiment){
  if("PathwayTables"%in%list.files(paste(experiment,"/",icdcode,"/",comparison,sep=""))==FALSE){dir.create(paste(paste(experiment,"/",icdcode,"/",comparison,sep=""),"/PathwayTables",sep=""))}
  files<-list.files(paste(paste(experiment,"/",icdcode,"/",comparison,sep=""),"/Pathways",sep=""))
  pathnames<-c() ; pathlist<-c() ; allpathnames<-c() ; allpathlist<-list()
  for(a in 1:length(files)){
    # a<-1
    files2<-list.files(paste(paste(experiment,"/",icdcode,"/",comparison,sep=""),"/Pathways/",files[a],sep=""))
    negfiles<-files2[intersect(grep("gsea_report_for_na_neg",files2),grep(".xls",files2))]
    neg<-read.csv2(paste(paste(experiment,"/",icdcode,"/",comparison,sep=""),"/Pathways/",files[a],"/",negfiles,sep=""),stringsAsFactors = F,sep="\t")
    posfiles<-files2[intersect(grep("gsea_report_for_na_pos",files2),grep(".xls",files2))]
    pos<-read.csv2(paste(paste(experiment,"/",icdcode,"/",comparison,sep=""),"/Pathways/",files[a],"/",posfiles,sep=""),stringsAsFactors = F,sep="\t")
    tabla<-rbind(pos[,c(1,4,6,8,9)],neg[,c(1,4,6,8,9)])
    tabla<-tabla[order(abs(as.numeric(as.character(tabla$NES))),decreasing = T),]
    write.table(tabla,paste(paste(experiment,"/",icdcode,"/",comparison,sep=""),"/PathwayTables/",gsub(".GseaPre.+","",files[a]),".txt",sep=""),quote=F,sep="\t",row.names = F)
    upi<-as.numeric(as.character(tabla$NES[which(as.numeric(as.character(tabla$NES))>0)]))
    names(upi)<-as.character(tabla$NAME[which(as.numeric(as.character(tabla$NES))>0)])
    downi<-as.numeric(as.character(tabla$NES[which(as.numeric(as.character(tabla$NES))<0)]))
    names(downi)<-as.character(tabla$NAME[which(as.numeric(as.character(tabla$NES))<0)])
    allpathlist[[gsub(".GseaPre.+","",files[a])]]$up<-upi
    allpathlist[[gsub(".GseaPre.+","",files[a])]]$down<-downi
    allpathnames<-c(allpathnames,tabla$NAME)
    ## Significantly enriched pathways ##
    signif<-which(as.numeric(as.character(tabla$FDR.q.val))<=0.05)
    if(length(signif)>0){
      pathnames<-c(pathnames,tabla$NAME[signif])
      indupreg<-which(as.numeric(as.character(tabla[signif,3]))>0)
      if(length(indupreg)>0){
        ups<-as.numeric(as.character(tabla[signif,3][indupreg]))
        names(ups)<-as.character(tabla[signif,1][indupreg])
        pathlist[[gsub(".GseaPre.+","",files[a])]]$up<-ups
      }
      inddownreg<-which(as.numeric(as.character(tabla[signif,3]))<0)
      if(length(inddownreg)>0){
        downs<-as.numeric(as.character(tabla[signif,3][inddownreg]))
        names(downs)<-as.character(tabla[signif,1][inddownreg])
        pathlist[[gsub(".GseaPre.+","",files[a])]]$down<-downs
      }
    }
  }
  print("All the tables have been saved")
  ## Create the matrixes for the heatmap ##
  paths<-unique(pathnames)
  allpaths<-unique(allpathnames)
  ## Binary and NES matrices ##
  binpathmat<-matrix(nrow=length(paths),ncol=length(pathlist),0)
  rownames(binpathmat)<-paths ; colnames(binpathmat)<-names(pathlist)
  nespathmat<-binpathmat
  for(a in 1:length(names(pathlist))){
    # a<-1
    if(length(which(names(pathlist[[a]])=="up"))>0){
      nespathmat[names(pathlist[[a]]$up),names(pathlist)[a]]<-as.numeric(pathlist[[a]]$up)
      binpathmat[names(pathlist[[a]]$up),names(pathlist)[a]]<-1
    }
    if(length(which(names(pathlist[[a]])=="down"))>0){
      nespathmat[names(pathlist[[a]]$down),names(pathlist)[a]]<-as.numeric(pathlist[[a]]$down)
      binpathmat[names(pathlist[[a]]$down),names(pathlist)[a]]<-(-1)
    }
  }
  ## New, include all the NES and not only significant (15/02/2022) ##
  allnespathmat<-matrix(nrow=length(allpaths),ncol=length(allpathlist),NA)
  rownames(allnespathmat)<-allpaths ; colnames(allnespathmat)<-names(allpathlist)
  for(a in 1:length(names(allpathlist))){
    if(length(which(names(allpathlist[[a]])=="up"))>0){
      nombresup<-names(allpathlist[[a]]$up)
      allnespathmat[nombresup,names(allpathlist)[a]]<-as.numeric(allpathlist[[a]]$up[nombresup])
    }
    if(length(which(names(allpathlist[[a]])=="down"))>0){
      nombresdown<-names(allpathlist[[a]]$down)
      allnespathmat[nombresdown,names(allpathlist)[a]]<-as.numeric(allpathlist[[a]]$down[nombresdown])
    }
  }
  print("Matrices for the heatmap have been created")
  write.table(nespathmat,paste(experiment,"/",icdcode,"/",comparison,"/Enrichments_nes.txt",sep=""),quote=F,sep="\t")
  write.table(allnespathmat,paste(experiment,"/",icdcode,"/",comparison,"/Enrichments_allnes.txt",sep=""),quote=F,sep="\t")
  write.table(binpathmat,paste(experiment,"/",icdcode,"/",comparison,"/Enrichments_binarized.txt",sep=""),quote=F,sep="\t")
  ## Save the results ##
  resultados<-list("nes"=nespathmat,"binarized"=binpathmat,"allnes"=allnespathmat)
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
## ICD9
info<-read.csv2("Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
sindupcol<-info[-(which(duplicated(info[,c(2,5)]))),c(2,5)]
colorcat9<-sindupcol[,2] ; names(colorcat9)<-sindupcol[,1]
## ICD10
info<-read.csv2("Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
sindupcol<-info[-(which(duplicated(info[,c(6,5)]))),c(6,5)]
colorcat10<-sindupcol[,2] ; names(colorcat10)<-sindupcol[,1]
## Colors for the pathways ##
patcatcol<-read.csv2("Data/final_colors_reactome_category.txt",stringsAsFactors = F,sep="\t")
reaccolors<-patcatcol[,2] ; names(reaccolors)<-patcatcol[,1]
## Identify the reactome parents ##
reactomeparents<-read.csv2("Data/Reactome_parents.txt",stringsAsFactors = F,sep="\t")
reactomeparents<-reactomeparents[order(reactomeparents[,2]),]
factreactomeparents<-reactomeparents[,2] ; names(factreactomeparents)<-reactomeparents[,1]

## Function to extract the tables ##
filtro<-"REACTOME"
comparison<-"Women"
icdcode<-"ICD9"
experiment<-"Microarrays"
nombre<-"Enrichments_allnes.txt"
prepare_tables<-function(comparison,icdcode,experiment,nombre,reactomeparents,filtro){
  tabla<-read.csv2(paste(experiment,"/",icdcode,"/",comparison,"/",nombre,sep=""),stringsAsFactors = F,sep="\t")
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

#### Plot clustering ####
filtro<-"REACTOME"
comparison<-"Men"
icdcode<-"ICD9"
experiment<-"Microarrays"
nombre<-"Enrichments_allnes.txt"
nombre2<-"allnes"
plotclusterings<-function(comparison,icdcode,experiment,nombre,reactomeparents,filtro,nombre2,colorcat){
  lallnes<-prepare_tables(comparison,icdcode,experiment,nombre,reactomeparents,filtro)
  femallnes<-lallnes$table
  pdf(file = paste(experiment,"/Plots/",icdcode,"_",comparison,"_",nombre2,".pdf",sep=""))
    ## Gender differences in cases ##
    if(icdcode=="ICD9"){
      colorcategory<-as.character(colorcat[gsub("ICD9__","ICD9_",colnames(femallnes))])
      colnames(femallnes)<-gsub("ICD9__","",colnames(femallnes))
    }
    colorpathways<-as.character(reaccolors[as.character(factreactomeparents[rownames(femallnes)])])
    heatmap.2(femallnes,
              distfun = function(x) dist(x, method="euclidean"),
              hclustfun = function(x) hclust(x, method="ward.D2"),
              rowsep = lallnes$separator,sepcolor = "gray",key=F,
              col=colorpanel(100, "#405191","#FFFFFF","#C81E17"),
              dendrogram = "col",Rowv=FALSE,labRow = "",cexCol = 0.7,
              scale="none", margins=c(5,5), trace="none",colCol = colorcategory,
              RowSideColors = colorpathways)
  dev.off()
}

plotclusterings("Women","ICD9","Microarrays","Enrichments_binarized.txt",reactomeparents,"REACTOME","binarized",colorcat9)
plotclusterings("Men","ICD9","Microarrays","Enrichments_binarized.txt",reactomeparents,"REACTOME","binarized",colorcat9)
plotclusterings("Women","ICD10","Microarrays","Enrichments_binarized.txt",reactomeparents,"REACTOME","binarized",colorcat10)
plotclusterings("Men","ICD10","Microarrays","Enrichments_binarized.txt",reactomeparents,"REACTOME","binarized",colorcat10)

#### Compare dendrograms ####
icdcode<-"ICD9"
experiment<-"Microarrays"
nombre<-"women_vs_men_clustering"
comparedendograms<-function(icdcode,experiment,nombre,firstcluster,secondcluster){
  comfirstcluster<-firstcluster[,intersect(colnames(firstcluster),colnames(secondcluster))]
  comsecondcluster<-secondcluster[,intersect(colnames(firstcluster),colnames(secondcluster))]
  if(icdcode=="ICD9"){
    colnames(comfirstcluster)<-gsub("ICD9__","",colnames(comfirstcluster))
    colnames(comsecondcluster)<-gsub("ICD9__","",colnames(comsecondcluster))
  }
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
  resultado<-list("cor"=correlation,"sim"=similarity)
  return(resultado)
}

## ICD9 ##
lallnes9<-prepare_tables("Women","ICD9","Microarrays","Enrichments_allnes.txt",reactomeparents,"REACTOME")
firstcluster<-lallnes9$table
lallnes9<-prepare_tables("Men","ICD9","Microarrays","Enrichments_allnes.txt",reactomeparents,"REACTOME")
secondcluster<-lallnes9$table
comparedendograms("ICD9","Microarrays","women_vs_men_clustering",firstcluster,secondcluster)

## ICD10 ##
lallnes10<-prepare_tables("Women","ICD10","Microarrays","Enrichments_allnes.txt",reactomeparents,"REACTOME")
firstcluster<-lallnes10$table
lallnes10<-prepare_tables("Men","ICD10","Microarrays","Enrichments_allnes.txt",reactomeparents,"REACTOME")
secondcluster<-lallnes10$table
comparedendograms("ICD10","Microarrays","women_vs_men_clustering",firstcluster,secondcluster)




















