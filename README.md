# Sex-specific transcriptome similarity networks elucidate comorbidity relationships
The goal of this project is to study the similarities and differences in gene expression between men and women across human diseases to better understand comorbidity differences between sexes.

To this end, we have downloaded publicly available microarray gene expression datasets from the Gene Expression Omnibus (GEO) and ArrayExpress, and the GREIN platform (http://www.ilincs.org/apps/grein/?gse=).

Since not all the studies register samples' sex, for those who doesn't, in the case of microarray data we extracted them using the MassiR package (https://academic.oup.com/bioinformatics/article/30/14/2084/2390865), which used the probes analyzing Y chromosome genes to identify samples' sex.
