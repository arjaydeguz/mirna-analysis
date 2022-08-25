if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install("miRNAtap")
BiocManager::install("miRNAtap.db")
BiocManager::install("topGO")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("annotate")
BiocManager::install("AnnotationDbi")
###################################################

#Make sure to install and load all of these packages

library(miRNAtap)
library(org.Hs.eg.db)
library(annotate)
library("AnnotationDbi")
library(Biobase)
library(miRNAtap.db)
library(topGO)

###################################################
#input your selected miRNA as the object "mir"

mir = 'miR-342-3p' 

predictions = getPredictedTargets(mir, species = 'hsa',
                                  method = 'geom', min_src = 2)
head(predictions)


############################################

#Getting the list of gene targets
list<-row.names(predictions)
targets<-getSYMBOL(list, data='org.Hs.eg')
targets<-as.data.frame(targets)
names(targets)<-c("SYMBOL")
row.names(targets)<-NULL
targets<-na.omit(targets)
head(targets)
dim(targets)

#export list of targets as csv
#change name according to miRNA


write.csv(targets,"Gene Targets of mir-342-3p.csv")

#note: list is already ordered according to signficance (based on whether gene-miRNA interaction was reported in 5 miRNA databases that was checked by miRNAtap)
#for your thesis: you can select the top 5 genes, and report them in a table with their known function. and then postulate how they are associated with the outcome that tested signficant for that miRNA

###################################################
#GO enrichment analysis

rankedGenes = predictions[,"rank_product"]
selection = function(x)TRUE
allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes =NULL, mapping='org.Hs.eg.db', ID='entrez')
GOdata = new('topGOdata', ontology = 'BP', allGenes = rankedGenes, annot = annFUN.GO2genes, GO2genes = allGO2genes, geneSel = selection, nodeSize=10)

results.ks = runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks

allRes = GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 50)
allRes[,c('GO.ID','Term','KS')]

head(allRes)

#export Gene Ontology results as csv
#change name according to miRNA

write.csv(allRes, "Gene Ontology Results for mir-342-3p.csv")

#note:list of GO terms is ranked according to p values using the 2 sample Kolmogorov Smirnov test (read up on KS test in Enrichment Analysis)
#for your thesis: only report the GO terms with a p value less than 0.05 and hypothesize significance in outcome that that tested significant for the miRNA

###################################################
#input your selected miRNA as the object "mir"

mir = 'miR-122-3p' 

predictions = getPredictedTargets(mir, species = 'hsa',
                                  method = 'geom', min_src = 2)
head(predictions)


############################################

#Getting the list of gene targets
list<-row.names(predictions)
targets<-getSYMBOL(list, data='org.Hs.eg')
targets<-as.data.frame(targets)
names(targets)<-c("SYMBOL")
row.names(targets)<-NULL
targets<-na.omit(targets)
head(targets)
dim(targets)

#export list of targets as csv
#change name according to miRNA


write.csv(targets,"Gene Targets of mir-122-3p.csv")

#note: list is already ordered according to signficance (based on whether gene-miRNA interaction was reported in 5 miRNA databases that was checked by miRNAtap)
#for your thesis: you can select the top 5 genes, and report them in a table with their known function. and then postulate how they are associated with the outcome that tested signficant for that miRNA

###################################################
#GO enrichment analysis

rankedGenes = predictions[,"rank_product"]
selection = function(x)TRUE
allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes =NULL, mapping='org.Hs.eg.db', ID='entrez')
GOdata = new('topGOdata', ontology = 'BP', allGenes = rankedGenes, annot = annFUN.GO2genes, GO2genes = allGO2genes, geneSel = selection, nodeSize=10)

results.ks = runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks

allRes = GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 50)
allRes[,c('GO.ID','Term','KS')]

head(allRes)

#export Gene Ontology results as csv
#change name according to miRNA

write.csv(allRes, "Gene Ontology Results for mir-122-3p.csv")

#note:list of GO terms is ranked according to p values using the 2 sample Kolmogorov Smirnov test (read up on KS test in Enrichment Analysis)
#for your thesis: only report the GO terms with a p value less than 0.05 and hypothesize significance in outcome that that tested significant for the miRNA

###################################################

#input your selected miRNA as the object "mir"

mir = 'miR-16-1-3p' 

predictions = getPredictedTargets(mir, species = 'hsa',
                                  method = 'geom', min_src = 2)
head(predictions)


############################################

#Getting the list of gene targets
list<-row.names(predictions)
targets<-getSYMBOL(list, data='org.Hs.eg')
targets<-as.data.frame(targets)
names(targets)<-c("SYMBOL")
row.names(targets)<-NULL
targets<-na.omit(targets)
head(targets)
dim(targets)

#export list of targets as csv
#change name according to miRNA


write.csv(targets,"Gene Targets of miR-16-1-3p.csv")

#note: list is already ordered according to signficance (based on whether gene-miRNA interaction was reported in 5 miRNA databases that was checked by miRNAtap)
#for your thesis: you can select the top 5 genes, and report them in a table with their known function. and then postulate how they are associated with the outcome that tested signficant for that miRNA

###################################################
#GO enrichment analysis

rankedGenes = predictions[,"rank_product"]
selection = function(x)TRUE
allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes =NULL, mapping='org.Hs.eg.db', ID='entrez')
GOdata = new('topGOdata', ontology = 'BP', allGenes = rankedGenes, annot = annFUN.GO2genes, GO2genes = allGO2genes, geneSel = selection, nodeSize=10)

results.ks = runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks

allRes = GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 50)
allRes[,c('GO.ID','Term','KS')]

head(allRes)

#export Gene Ontology results as csv
#change name according to miRNA

write.csv(allRes, "Gene Ontology Results for miR-16-1-3p.csv")

#note:list of GO terms is ranked according to p values using the 2 sample Kolmogorov Smirnov test (read up on KS test in Enrichment Analysis)
#for your thesis: only report the GO terms with a p value less than 0.05 and hypothesize significance in outcome that that tested significant for the miRNA
###################################################



#input your selected miRNA as the object "mir"

mir = 'miR-223-5p' 

predictions = getPredictedTargets(mir, species = 'hsa',
                                  method = 'geom', min_src = 2)
head(predictions)


############################################

#Getting the list of gene targets
list<-row.names(predictions)
targets<-getSYMBOL(list, data='org.Hs.eg')
targets<-as.data.frame(targets)
names(targets)<-c("SYMBOL")
row.names(targets)<-NULL
targets<-na.omit(targets)
head(targets)
dim(targets)

#export list of targets as csv
#change name according to miRNA


write.csv(targets,"Gene Targets of miR-223-5p.csv")

#note: list is already ordered according to signficance (based on whether gene-miRNA interaction was reported in 5 miRNA databases that was checked by miRNAtap)
#for your thesis: you can select the top 5 genes, and report them in a table with their known function. and then postulate how they are associated with the outcome that tested signficant for that miRNA

###################################################
#GO enrichment analysis

rankedGenes = predictions[,"rank_product"]
selection = function(x)TRUE
allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes =NULL, mapping='org.Hs.eg.db', ID='entrez')
GOdata = new('topGOdata', ontology = 'BP', allGenes = rankedGenes, annot = annFUN.GO2genes, GO2genes = allGO2genes, geneSel = selection, nodeSize=10)

results.ks = runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks

allRes = GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 50)
allRes[,c('GO.ID','Term','KS')]

head(allRes)

#export Gene Ontology results as csv
#change name according to miRNA

write.csv(allRes, "Gene Ontology Results for miR-223-5p.csv")

#note:list of GO terms is ranked according to p values using the 2 sample Kolmogorov Smirnov test (read up on KS test in Enrichment Analysis)
#for your thesis: only report the GO terms with a p value less than 0.05 and hypothesize significance in outcome that that tested significant for the miRNA
###################################################

#input your selected miRNA as the object "mir"

mir = 'miR-297' 

predictions = getPredictedTargets(mir, species = 'hsa',
                                  method = 'geom', min_src = 2)
head(predictions)


############################################

#Getting the list of gene targets
list<-row.names(predictions)
targets<-getSYMBOL(list, data='org.Hs.eg')
targets<-as.data.frame(targets)
names(targets)<-c("SYMBOL")
row.names(targets)<-NULL
targets<-na.omit(targets)
head(targets)
dim(targets)

#export list of targets as csv
#change name according to miRNA


write.csv(targets,"Gene Targets of miR-297.csv")

#note: list is already ordered according to signficance (based on whether gene-miRNA interaction was reported in 5 miRNA databases that was checked by miRNAtap)
#for your thesis: you can select the top 5 genes, and report them in a table with their known function. and then postulate how they are associated with the outcome that tested signficant for that miRNA

###################################################
#GO enrichment analysis

rankedGenes = predictions[,"rank_product"]
selection = function(x)TRUE
allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes =NULL, mapping='org.Hs.eg.db', ID='entrez')
GOdata = new('topGOdata', ontology = 'BP', allGenes = rankedGenes, annot = annFUN.GO2genes, GO2genes = allGO2genes, geneSel = selection, nodeSize=10)

results.ks = runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks

allRes = GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 50)
allRes[,c('GO.ID','Term','KS')]

head(allRes)

#export Gene Ontology results as csv
#change name according to miRNA

write.csv(allRes, "Gene Ontology Results for miR-297.csv")

#note:list of GO terms is ranked according to p values using the 2 sample Kolmogorov Smirnov test (read up on KS test in Enrichment Analysis)
#for your thesis: only report the GO terms with a p value less than 0.05 and hypothesize significance in outcome that that tested significant for the miRNA
###################################################





#input your selected miRNA as the object "mir"

mir = 'miR-4772-3p' 

predictions = getPredictedTargets(mir, species = 'hsa',
                                  method = 'geom', min_src = 2)
head(predictions)


############################################

#Getting the list of gene targets
list<-row.names(predictions)
targets<-getSYMBOL(list, data='org.Hs.eg')
targets<-as.data.frame(targets)
names(targets)<-c("SYMBOL")
row.names(targets)<-NULL
targets<-na.omit(targets)
head(targets)
dim(targets)

#export list of targets as csv
#change name according to miRNA


write.csv(targets,"Gene Targets of miR-4772-3p.csv")

#note: list is already ordered according to signficance (based on whether gene-miRNA interaction was reported in 5 miRNA databases that was checked by miRNAtap)
#for your thesis: you can select the top 5 genes, and report them in a table with their known function. and then postulate how they are associated with the outcome that tested signficant for that miRNA

###################################################
#GO enrichment analysis

rankedGenes = predictions[,"rank_product"]
selection = function(x)TRUE
allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes =NULL, mapping='org.Hs.eg.db', ID='entrez')
GOdata = new('topGOdata', ontology = 'BP', allGenes = rankedGenes, annot = annFUN.GO2genes, GO2genes = allGO2genes, geneSel = selection, nodeSize=10)

results.ks = runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks

allRes = GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 50)
allRes[,c('GO.ID','Term','KS')]

head(allRes)

#export Gene Ontology results as csv
#change name according to miRNA

write.csv(allRes, "Gene Ontology Results for miR-4772-3p.csv")

#note:list of GO terms is ranked according to p values using the 2 sample Kolmogorov Smirnov test (read up on KS test in Enrichment Analysis)
#for your thesis: only report the GO terms with a p value less than 0.05 and hypothesize significance in outcome that that tested significant for the miRNA
###################################################






#input your selected miRNA as the object "mir"

mir = 'miR-574-5p' 

predictions = getPredictedTargets(mir, species = 'hsa',
                                  method = 'geom', min_src = 2)
head(predictions)


############################################

#Getting the list of gene targets
list<-row.names(predictions)
targets<-getSYMBOL(list, data='org.Hs.eg')
targets<-as.data.frame(targets)
names(targets)<-c("SYMBOL")
row.names(targets)<-NULL
targets<-na.omit(targets)
head(targets)
dim(targets)

#export list of targets as csv
#change name according to miRNA


write.csv(targets,"Gene Targets of miR-574-5p.csv")

#note: list is already ordered according to signficance (based on whether gene-miRNA interaction was reported in 5 miRNA databases that was checked by miRNAtap)
#for your thesis: you can select the top 5 genes, and report them in a table with their known function. and then postulate how they are associated with the outcome that tested signficant for that miRNA

###################################################
#GO enrichment analysis

rankedGenes = predictions[,"rank_product"]
selection = function(x)TRUE
allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes =NULL, mapping='org.Hs.eg.db', ID='entrez')
GOdata = new('topGOdata', ontology = 'BP', allGenes = rankedGenes, annot = annFUN.GO2genes, GO2genes = allGO2genes, geneSel = selection, nodeSize=10)

results.ks = runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks

allRes = GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 50)
allRes[,c('GO.ID','Term','KS')]

head(allRes)

#export Gene Ontology results as csv
#change name according to miRNA

write.csv(allRes, "Gene Ontology Results for miR-574-5p.csv")

#note:list of GO terms is ranked according to p values using the 2 sample Kolmogorov Smirnov test (read up on KS test in Enrichment Analysis)
#for your thesis: only report the GO terms with a p value less than 0.05 and hypothesize significance in outcome that that tested significant for the miRNA
###################################################





