# A catalog of sex differences in human diseases
The goal of this project is to study the similarities and differences in gene expression between men and women across human diseases to better understand comorbidity differences between sexes.

To this end, we have downloaded publicly available microarray and RNAseq gene expression datasets from the Gene Expression Omnibus (GEO), ArrayExpress, and the GREIN platform (http://www.ilincs.org/apps/grein/?gse=). The first two provide raw gene expression data (microarrays), while the second one provides the transcript counts (https://www.nature.com/articles/s41598-019-43935-8).

Since not all the studies register samples' sex, for those who doesn't, in the case of microarray data we extracted them using the MassiR package (https://academic.oup.com/bioinformatics/article/30/14/2084/2390865), which used the probes analyzing Y chromosome genes to identify samples' sex.
