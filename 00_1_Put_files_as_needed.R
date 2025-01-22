## CODE TO ADAPT THE FILES AS NEEDED FOR THE ANALYSIS ##
## Developed by Jon Sanchez-Valle
## Barcelona Supercomputing Center
## Life Science Department
## Computational Biology Group
## Email: jon.sanchez@bsc.es

## Adaptations ##
info<-read.csv2("Microarrays/Data/Disease_information.txt",stringsAsFactors = F,sep="\t")

#### ICD9 names ####
icd9codes<-unique(info$ICD9)

icd9namesgiven<-c("Acne","Diseases of white blood cells","Atopic dermatitis","Acute myeloid leukemia","Acute myocardial infarction",
             "Adrenocortical carcinoma","Lymphoid leukemia","Other phenotypes","Alcoholic hepatitis","Alopecia",
             "Alzheimer's disease","Amyotrophic lateral sclerosis","ANCA crescentic glomerulonephritis","Oral diseases","Ventricular cardiomyopathy",
             "Asthma","Brain tumors","Atrial fibrillation","Autism","Axial spondyloarthropathy",
             "Lung cancer","Behcet's disease","Bladder cancer","Breast cancer","Campylobacter jenuni infection",
             "Stroke","Congenital anomalies","Cervix cancer","Choroideremia","COPD",
             "Chronic rhinosinusitis","Unspecified disorders of metabolism","Colon cancer","Crohns disease","Sarcoidosis",
             "Dengue","Diffuse diseases of connective tissue","Endometriosis","Esophagus cancer","Other neoplasms",
             "Fatal familian insomnia","Lymphoma","Gastric cancer","Viral hepatitis","Liver cancer",
             "Human papillomavirus","Huntington's disease","Progeria","Other alveolar and parietoalveolar pneumonopathy","Cystitis",
             "Irritable bowel syndrome","Ischemia","Chromosomal anomalies (turner / klinefelter)","Leigh syndrome","Leishmaniasis",
             "Mood disorders (depression / bipolar disorder)","Pleura cancer","MGUS","Multiple myeloma","Multiple sclerosis",
             "Tuberculosis (mycobacterium)","Myositis","Nasopharingeal cancer","Oral cavity cancer","Osteoarthritis",
             "Osteosarcoma","Ovary cancer","Pancreas cancer","Parkinson's disease","Polycystic ovary syndrome",
             "Peripheral arterial disease","Polen allergy","Primary myelofibrosis","Prostate cancer","Psoriasis",
             "Kidney cancer","Rosacea","Schizophrenia","Seborrheic keratosis","Septic shock",
             "Setleis syndrome","Hereditary hemolytic anemias (Sickle / Thalassemia)","Smoking","Sotos syndrome","Spina bifida",
             "Thrombocytopenia","Thyroid cancer","Muscular dystrophy","Trachoma","Diabetes mellitus",
             "Ulcerative enterocolitis","Uremia","Vitiligo","Lung allograft infection by aspergillus","Lung allograft dysfunction",
             "IgA nephropathy","Melanoma","Melanoma metastasis","Polymyalgia rheumatica")

icd9namesreal<-c("Diseases of sebaceous glands","Diseases of white blood cells","Atopic dermatitis and related conditions","Myeloid leukemia acute","Acute myocardial infarction",
                  "Malignant neoplasm of other endocrine glands and related structures","Lymphoid leukemia","Other phenotypes","Chronic liver disease and cirrhosis","Diseases of hair and hair follicles",
                 "Other cerebral degenerations","Anterior horn cell disease","Acute glomerulonephritis","Diseases of the oral soft tissues excluding lesions specific for gingiva and tongue","Cardiomyopathy",
                 "Asthma","Malingnant neoplasm of brain","Cardiac dysrhythmias","Pervasive developmental disorders","Spondylosis and allied disorders",
                 "Malignant neoplasm of trachea bronchus and lung","Other and unspecified infectious and parasitic diseases","Malignant neoplasm of bladder","Malignant neoplasm of female breast","Intestinal infections due to other organisms",
                 "Occlusion of cerebral arteries","Other and unspecified congenital anomalies","Malignant neoplasm of cervix uteri","Chorioretinal inflammations scars and other disorders of choroid","Chronic airway obstruction, not elsewhere classified",
                 "Chronic sinusitis","Other and unspecified disorders of metabolism","Malignant neoplasm of colon","Regional enteritis","Sarcoidosis",
                 "Dengue","Diffuse diseases of connective tissue","Endometriosis","Malignant neoplasm of esophagus","Neoplasm of uncertain behavior of other and unspecified sites and tissues",
                 "Slow virus infection and prion diseases of central nervous system","Other malignant neoplasms of lymphoid and histiocytic tissue","Malignant neoplasm of stomach","Viral hepatitis","Malignant neoplasm of liver and intrahepatic bile ducts",
                 "Viral and chlamydial infection in conditions classified elsewhere and of unspecified site","Other extrapyramidal disease and abnormal movement disorders","Other endocrine disorders","Other alveolar and parietoalveolar pneumonopathy","Cystitis",
                 "Functional digestive disorders not elsewhere classified","Other forms of chronic ischemic heart disease","Chromosomal anomalies","Cerebral degenerations usually manifest in childhood","Leishmaniasis",
                 "Episodic mood disorders","Malignant neoplasm of pleura","Disorders of plasma protein metabolism","Multiple myeloma and immunoproliferative neoplasms","Multiple sclerosis",
                 "Diseases due to other mycobacteria","Other disorders of soft tissues","Malignant neoplasm of nasopharynx","Malignant neoplasm of other and ill-defined sites within the lip oral cavity and pharynx","Osteoarthrosis and allied disorders",
                 "Malignant neoplasm of bone and articular cartilage","Malignant neoplasm of ovary and other uterine adnexa","Malignant neoplasm of pancreas","Parkinson's disease","Ovarian dysfunction",
                 "Other peripheral vascular disease","Allergic rhinitis","Other diseases of blood and blood-forming organs","Malignant neoplasm of prostate","Psoriasis and similar disorders",
                 "Malignant neoplasm of kidney and other and unspecified urinary organs","Erythematous conditions","Schizophrenic disorders","Other dermatoses","Symptoms involving cardiovascular system",
                 "Congenital anomalies of the integument","Hereditary hemolytic anemias","Nondependent abuse of drugs","Disorders of the pituitary gland and its hypothalamic control","Spina bifida",
                 "Purpura and other hemorrhagic conditions","Malignant neoplasm of thyroid gland","Muscular dystrophies and other myopathies","Hemarthrosis, unspecified foot","Diabetes mellitus",
                 "Ulcerative enterocolitis","Renal failure, unspecified","Other disorders of skin and subcutaneous tissue","Pneumonia in infectious diseases classified elsewhere","Complications peculiar to certain specified procedures",
                 "Nephritis and nephropathy not specified as acute or chronic","Malignant melanoma of skin","Secondary malignant neoplasm of other specified sites","Polymyalgia rheumatica")


info[which(info$ICD9==icd9codes[99]),1:2]

length(icd9namesgiven)
length(icd9namesreal)

names(icd9namesgiven)<-icd9codes
names(icd9namesreal)<-icd9codes

info<-cbind(info,as.character(icd9namesgiven[info$ICD9]),as.character(icd9namesreal[info$ICD9]))
colnames(info)[7:8]<-c("ICD9_our_names","ICD9_name")

write.table(info,"Microarrays/Data/Disease_information.txt",quote = F,sep="\t",row.names = F)

#### ICD10 names ####
icd10codes<-unique(info$ICD10)

info[which(info$ICD10==icd10codes[3]),c(1,6)]

icd10namesgiven<-c("Acne","Hemophagocytic lymphohistiocytosis","Atopic dermatitis","Myeloid leukemia","Acute myocardial infarction",
                  "Adrenal cancer","Lymphoid leukemia","","Alcoholic liver disease","Alopecia",
                  "Alzheimer's disease","Amyotrophic lateral sclerosis","ANCA crescentic glomerulonephritis","Aphthous stomatitis","Cardiomyopathy",
                  "Asthma","Brain cancer","Atrial fibrillation","Autism","Autosomal dominant monocytopenia",
                  "Axial Spondyloarthropathy","Lung cancer","Other systemic involvement of connective tissue","Bladder cancer","Breast cancer",
                  "Campylobacter Jenuni Infection","Cerebral infarction","Other specified congenital malformation syndromes affecting multiple systems","Cervix cancer","Choroideremia",
                  "Chronic Obstructive Pulmonary Disease","Chronic sinusitis","Other and unspecified metabolic disorders","Colon cancer","Crohn's disease",
                  "Sarcoidosis","Dengue fever","Dermatomyositis","Endometriosis","Esophageal cancer",
                  "Other neoplasms of uncertain behavior of lymphoid, hematopoietic and related tissue","Fatal Familial Insomnia","Follicular lymphoma","Gastric cancer","Hepatitis C",
                  "Hepatitis B","Liver cancer","human papillomavirus","Huntington's disease","Hutchinson Gilford Progeria",
                  "Other interstitial pulmonary diseases","Cystitis","Irritable bowel syndrome","Ischemic Heart","ISCU Myopathy",
                  "JMML","Job's syndrome","Klinefelter's syndrome","Leigh syndrome","Leishmaniasis",
                  "Major depression","Mesothelioma","Myeloma","Multiple sclerosis","Mycobacterium Tuberculosis",
                  "Myelodysplastic syndromes","Myositis","Nasopharingeal cancer","Oral cavity cancer","Oral dysplasia",
                  "Osteoarthritis","Osteosarcoma","Ovary cancer","Pancreatic cancer","Parkinson's disease",
                  "PCOS","Peripheral arterial disease","Polen allergy","Polycythemia vera","Primary myelofibrosis",
                  "Prostate cancer","Psoriasis","Kidney cancer","Rhabdoid cancer","Rosacea",
                  "Schizophrenia","Seborrheic keratosis","Septic Shock","Setleis syndrome","Sickle-cell anemia",
                  "Systemic lupus erythematosus","Smoker","Spina bifida","T-cell lymphoma","Thalassemia",
                  "Thrombocytopenia","Thyroid cancer","Tibial muscular dystrophy","Trachoma","Turner's syndrome",
                  "Type 1 diabetes mellitus","Type 2 diabetes mellitus","Ulcerative colitis","Uremia","Vitiligo",
                  "Aspergillus colonization of lung allograft","Bipolar disorder","Chronic lung allograft dysfunction","Fallopian tube epithelium carcinoma","IgA nephropathy",
                  "Melanoma","Melanoma metastasis","Respiratory syncytial virus")

icd10namesreal<-c("Acne","Other specified diseases with participation of lymphoreticular and reticulohistiocytic tissue","Atopic dermatitis","Myeloid leukemia","Acute myocardial infarction",
                 "Malignant neoplasm of adrenal gland","Lymphoid leukemia","","Alcoholic liver disease","Other nonscarring hair loss",
                 "Alzheimer's disease","Spinal muscular atrophy and related syndromes","Rapidly progressive nephritic syndrome","Stomatitis and related lesions","Cardiomyopathy",
                 "Asthma","Malignant neoplasm of brain","Atrial fibrillation and flutter","Pervasive developmental disorders","Other disorders of white blood cells",
                 "Spondylosis","Malignant neoplasm of bronchus and lung","Other systemic involvement of connective tissue","Malignant neoplasm of bladder","Malignant neoplasm of breast",
                 "Other bacterial intestinal infections","Cerebral infarction","Other specified congenital malformation syndromes affecting multiple systems","Malignant neoplasm of cervix uteri","Other disorders of choroid",
                 "Other chronic obstructive pulmonary disease","Chronic sinusitis","Other and unspecified metabolic disorders","Malignant neoplasm of colon","Crohn's disease",
                 "Sarcoidosis","Dengue fever","Dermatopolymyositis","Endometriosis","Malignant neoplasm of esophagus",
                 "Other neoplasms of uncertain behavior of lymphoid, hematopoietic and related tissue","Atypical virus infections of central nervous system","Follicular lymphoma","Malignant neoplasm of stomach","Other acute viral hepatitis",
                 "Unspecified viral hepatitis","Malignant neoplasm of liver and intrahepatic bile ducts","Other predominantly sexually transmitted diseases, not elsewhere classified","Huntington's disease","Other endocrine disorders",
                 "Other interstitial pulmonary diseases","Cystitis","Irritable bowel syndrome","Chronic ischemic heart disease","Other and unspecified myopathies",
                 "Monocytic leukemia","Immunodeficiency associated with other major defects","Other sex chromosome abnormalities, male phenotype, not elsewhere classified","Other degenerative diseases of nervous system, not elsewhere classified","Leishmaniasis",
                 "Major depressive disorder, recurrent","Malignant neoplasm of heart, mediastinum and pleura","Multiple myeloma and malignant plasma cell neoplasms","Multiple sclerosis","Respiratory tuberculosis",
                 "Myelodysplastic syndromes","Myositis","Malignant neoplasm of nasopharynx","Malignant neoplasm of other and ill-defined sites in the lip, oral cavity and pharynx","Other diseases of lip and oral mucosa",
                 "Other and unspecified osteoarthritis","Malignant neoplasm of bone and articular cartilage of other and unspecified sites","Malignant neoplasm of ovary","Malignant neoplasm of pancreas","Parkinson's disease",
                 "Ovarian dysfunction","Other peripheral vascular diseases","Vasomotor and allergic rhinitis","Polycythemia vera","Other and unspecified diseases of blood and blood-forming organs",
                 "Malignant neoplasm of prostate","Psoriasis","Malignant neoplasm of kidney, except renal pelvis","Malignant neoplasm of other connective and soft tissue","Rosacea",
                 "Schizophrenia","Seborrheic keratosis","Symptoms and signs specifically associated with systemic inflammation and infection","Other congenital malformations of skin","Sickle-cell disorders",
                 "Systemic lupus erythematosus (SLE)","Nicotine dependence","Spina bifida","Mature T/NK-cell lymphomas","Thalassemia",
                 "Purpura and other hemorrhagic conditions","Malignant neoplasm of thyroid gland","Primary disorders of muscles","Trachoma","Turner's syndrome",
                 "Type 1 diabetes mellitus","Type 2 diabetes mellitus","Ulcerative colitis","Unspecified kidney failure","Vitiligo",
                 "Aspergillosis","Bipolar disorder","Complications of transplanted organs and tissue","Malignant neoplasm of other and unspecified female genital organs","Recurrent and persistent hematuria",
                 "Malignant melanoma of skin","Secondary malignant neoplasm of other and unspecified sites","Viral agents as the cause of diseases classified elsewhere")


length(icd10namesgiven)
length(icd10namesreal)

names(icd10namesgiven)<-icd10codes
names(icd10namesreal)<-icd10codes

info<-cbind(info,as.character(icd10namesgiven[info$ICD10]),as.character(icd10namesreal[info$ICD10]))
colnames(info)[9:10]<-c("ICD10_our_names","ICD10_name")

write.table(info,"Microarrays/Data/Disease_information.txt",quote=F,sep="\t",row.names=F)

## Add the categories based on ICD10 ##
info<-read.csv2("Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
category<-info[,c(2,6,4,5)]
if(length(which(duplicated(category)))>0){category<-category[-which(duplicated(category)),]}
category<-category[order(category$ICD10),]

newcats<-c("Unclassified diseases and risk factors",rep("Certain infectious and parasitic diseases",11),rep("Neoplasms",33),
  rep("Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism",8),
  rep("Endocrine, nutritional and metabolic diseases",5),rep("Mental and behavioural disorders",5),
  rep("Diseases of the nervous system",8),"Diseases of the eye and adnexa",rep("Diseases of the circulatory system",6),
  rep("Diseases of the respiratory system",5),rep("Diseases of the digestive system",6),rep("Diseases of the skin and subcutaneous tissue",7),
  rep("Diseases of the musculoskeletal system and connective tissue",8),rep("Diseases of the genitourinary system",5),
  rep("Congenital malformations, deformations and chromosomal abnormalities",6),
  "Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified","Injury, poisoning and certain other consequences of external causes")
names(newcats)<-category$ICD10

newcolors<-c("#97F52C",rep("#4EF9EA",11),rep("#570892",33),rep("#758F55",8),rep("#C55F1B",5),rep("#03A412",5),
           rep("#75C0D9",8),"#FFC9DE",rep("#0BF57B",6),rep("#FE1659",5),rep("#82DD7A",6),rep("#5045F3",7),
           rep("#1E5448",8),rep("#42223B",5),rep("#272399",6),"#D1EB04","#DCD41E")
names(newcolors)<-category$ICD10

info2<-cbind(info,as.character(newcats[info$ICD10]),as.character(newcolors[info$ICD10]))

colnames(info2)[c(4,5,11,12)]<-c("ICD9_category","ICD9_colors","ICD10_category","ICD10_colors")

info2<-info2[,c(1,3,2,7:8,4:5,6,9:10,11:12)]
write.table(info2,"Microarrays/Data/Disease_information.txt",quote=F,sep="\t",row.names=F)


#### Prepare Reactome tables ####
## Identify the reactome parents and select their genes ##
reactomeparents<-read.csv2("Microarrays/Data/Reactome_parents.txt",stringsAsFactors = F,sep="\t")
reactomeparents<-reactomeparents[order(reactomeparents[,2]),]
factreactomeparents<-reactomeparents[,2] ; names(factreactomeparents)<-reactomeparents[,1]
## Get the genes associated to each pathway ##
gseareac<-read.csv2("Microarrays/Data/c2.cp.reactome.v7.0.symbols.gmt",stringsAsFactors = F)
pathwaygenelist<-list()
for(a in 1:length(gseareac[,1])){
  # a<-1
  pathgenevector<-strsplit(gseareac[a,1],"\t")[[1]]
  pathwaygenelist[[pathgenevector[1]]]<-pathgenevector[3:length(pathgenevector)]
}
## Get the genes associated to each reactome parent ##
theparents<-unique(reactomeparents$Parent)
reactomeparentgenelist<-list()
genesperparent<-c()
for(a in 1:length(theparents)){
  # a<-1
  thepathways<-reactomeparents$GSEA_pathway[which(reactomeparents$Parent==theparents[a])]
  thegenes<-c();for(t in 1:length(thepathways)){thegenes<-c(thegenes,pathwaygenelist[[thepathways[t]]])}
  reactomeparentgenelist[[theparents[a]]]<-unique(thegenes)
  genesperparent<-rbind(genesperparent,c(theparents[a],length(unique(thegenes))))
}
genesperparent<-genesperparent[order(as.numeric(genesperparent[,2]),decreasing = T),]

## Reactome child-parent relation ##
parentchild<-read.csv2("Microarrays/Data/Reactome/ReactomePathwaysRelation.txt",stringsAsFactors = F,sep="\t",header = F)
parentchild<-parentchild[intersect(grep("HSA",parentchild$V1),grep("HSA",parentchild$V2)),]
reacnameid<-read.csv2("Microarrays/Data/Reactome/ReactomePathways.txt",stringsAsFactors = F,sep="\t",header=F)
reacnameid<-reacnameid[grep("HSA",reacnameid$V1),]
namereacnameid<-reacnameid[,2] ; names(namereacnameid)<-reacnameid[,1]

welostl<-list()
secondparents<-list() ; secondparentnumbers<-c()
for(a in 1:length(genesperparent[,1])){
  # a<-1
  theid<-reacnameid$V1[which(reacnameid$V2==genesperparent[a,1])]
  theids2<-parentchild$V2[which(parentchild$V1==theid)]
  thenames2<-c() ; for(b in 1:length(theids2)){thenames2<-c(thenames2,reacnameid$V2[which(reacnameid$V1==theids2[b])])}
  if(length(theids2)!=length(thenames2)){print(paste(a,"different numbers of ids and names"))}
  thenames2_1<-paste("REACTOME_",toupper(gsub("-","_",gsub(" ","_",thenames2))),sep="")
  if(length(intersect(thenames2_1,names(pathwaygenelist)))!=length(thenames2)){welostl[[genesperparent[a,1]]]<-setdiff(thenames2_1,names(pathwaygenelist))}
  for(b in 1:length(thenames2_1)){
    # b<-1
    cual<-which(names(pathwaygenelist)==thenames2_1[b])
    if(length(cual)>0){
      secondparents[[thenames2[b]]]<-pathwaygenelist[[cual]]
      secondparentnumbers<-rbind(secondparentnumbers,c(genesperparent[a,1],theids2[b],thenames2[b],thenames2_1[b],length(pathwaygenelist[[cual]])))
    }
  }
}
secondnames<-secondparentnumbers[,1] ; names(secondnames)<-secondparentnumbers[,3]
## For each pathway we are going to obtain the entire list of sons iteratively ##
secondparentschilds<-c()
for(a in 1:length(secondparentnumbers[,1])){
  allchildren<-c()
  # a<-1
  childs<-parentchild$V2[which(parentchild$V1==secondparentnumbers[a,2])]
  while(length(childs)>0){
    allchildren<-c(allchildren,childs)
    childs2<-c(); for(b in 1:length(childs)){childs2<-c(childs2,parentchild$V2[which(parentchild$V1==childs[b])])}
    childs<-unique(childs2)
  }
  allchildren<-c(secondparentnumbers[a,2],allchildren)
  secondparentschilds<-rbind(secondparentschilds,cbind(secondparentnumbers[a,1],secondparentnumbers[a,2],secondparentnumbers[a,3],allchildren))
}
if(dim(secondparentschilds)[2]==4){
  secondparentschilds<-cbind(secondparentschilds,as.character(namereacnameid[secondparentschilds[,4]]),paste("REACTOME_",toupper(gsub("-","_",gsub(" ","_",as.character(namereacnameid[secondparentschilds[,4]])))),sep=""))
}

segundoshijos<-unique(secondparentschilds[,3])
secondparentgenes<-c()
for(a in 1:length(segundoshijos)){
  # a<-1
  thepathways<-secondparentschilds[which(secondparentschilds[,3]==segundoshijos[a]),6]
  interseccion<-intersect(thepathways,names(pathwaygenelist))
  if(length(interseccion)>0){
    allgenes<-c()
    for(b in 1:length(interseccion)){
      allgenes<-c(allgenes,pathwaygenelist[[interseccion[b]]])
    }
    allgenes<-unique(allgenes)
    secondparentgenes<-rbind(secondparentgenes,cbind(as.character(secondnames[segundoshijos[a]]),segundoshijos[a],allgenes))
  }
}
colnames(secondparentgenes)<-c("Parent","Second","Genes")
write.table(secondparentgenes,"Microarrays/Data/Reactome_second_parents.txt",quote = F,sep="\t",row.names=F)


#### Just a table with the colors by disease - For Iker ####
## ICD10 ##
disinf<-read.csv2("Microarrays/Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
if(dim(disinf)[2]==12){disinf<-disinf[,c(8,11,12)]}
if(length(which(duplicated(disinf)))>0){disinf<-disinf[-which(duplicated(disinf)),]}
if(length(which(disinf$ICD10=="-"))>0){disinf<-disinf[-which(disinf$ICD10=="-"),]}
disinf<-disinf[order(disinf$ICD10,decreasing = F),]
write.table(disinf,"Microarrays/Data/ICD10_colors.txt",quote=F,sep="\t",row.names = F)

## ICD9 ##
disinf<-read.csv2("Microarrays/Data/Disease_information.txt",stringsAsFactors = F,sep="\t")
if(dim(disinf)[2]==12){disinf<-disinf[,c(3,6,7)]}
if(length(which(duplicated(disinf)))>0){disinf<-disinf[-which(duplicated(disinf)),]}
if(length(which(disinf$ICD10=="-"))>0){disinf<-disinf[-which(disinf$ICD10=="-"),]}
disinf<-disinf[order(disinf$ICD10,decreasing = F),]
write.table(disinf,"Microarrays/Data/ICD9_colors.txt",quote=F,sep="\t",row.names = F)

## Select only those interactions from IID that have been experimentally validated ##
tt<-fread("DrugAnalysis/human_annotated_PPIs.txt",stringsAsFactors = F)
iid<-tt[grep("exp",tt$evidence_type)]
if(dim(iid)[2]==319){iid<-iid[,c(1:10)]}
write.table(iid,"DrugAnalysis/PPIN.txt",quote=F,sep="\t",row.names = F)

## New colors for the Reactome categories ##
reaccol<-fread("Microarrays/Data/final_colors_reactome_category.txt",stringsAsFactors = F,sep="\t")

newcols<-c("#641E16","#E74C3C","#D98880","#512E5F","#8E44AD",
           "#C39BD3","#154360","#3498DB","#7FB3D5","#0E6251",
           "#16A085","#76D7C4","#145A32","#2ECC71","#7DCEA0",
           "#7D6608","#F39C12","#F7DC6F","#784212","#D35400",
           "#F0B27A","#7B7D7D","#BDC3C7","#F4F6F7","#4D5656",
           "#7F8C8D","#BFC9CA","#1B2631","#2C3E50","#85929E")

reaccol$Color<-newcols
write.table(reaccol,"Microarrays/Data/new_final_colors_reactome_category.txt",quote=F,sep="\t",row.names=F)




