
######################packages########################
library(vegan)
library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(SIAMCAT)

##########################remove batch effect using ConQuR##########################
####################Risk based removal###############
library(ConQuR)
library(doParallel) 
library(reshape2)
library("LiblineaR")

L6dataPublic<-read.csv("Risk Model Data/Before_ConQuR_60public_count.csv",header=TRUE)
L6data74<-read.csv("Risk Model Data/Before_ConQuR_74samples_count.csv",header=TRUE)
L6data34<-read.csv("Risk Model Data/Before_ConQuR_34samples_count.csv",header=TRUE)
batch<-read.csv("Risk Model Data/16Slist_all_batch3.csv",header=TRUE)

dataPublic<-melt(L6dataPublic)
data74<-melt(L6data74)
data34<-melt(L6data34)

dataALL<-rbind(dataPublic,data74,data34)
colnames(dataALL)<-c("SampleID","Genus","Count")
dataALL<-dataALL[dataALL$Count>2,]
dataALL<-merge(batch,dataALL,by="SampleID",all=FALSE)
dataALL_matrix<-tapply(dataALL$Count,list(dataALL$SampleID,dataALL$Genus),median)

dataALL_matrix_WithGroup<-cbind(dataALL_matrix,batch,by="SampleID")
dataALL_matrix_WithGroup[is.na(dataALL_matrix_WithGroup)]<-0

write.csv(dataALL_matrix_WithGroup,"data_matrix_WithGroup.csv")
#dataALL_matrix_WithGroup<-read.csv("dataALL_matrix_WithGroup.csv",header=TRUE)
#rownames(dataALL_matrix_WithGroup)<-dataALL_matrix_WithGroup$SampleID
#dataALL_matrix_WithGroup<-dataALL_matrix_WithGroup[,-1]

taxa = dataALL_matrix_WithGroup[, 1:1377]
batchid = as.factor(dataALL_matrix_WithGroup[, 'Country'])
covar = dataALL_matrix_WithGroup[, c("Group","Gender","Age")]

options(warn=-1)
taxa_corrected1 = ConQuR(tax_tab=taxa, batchid=batchid, covariates=covar, batch_ref="1")

taxa_corrected1[146:150, 1:5]

options(warn=-1)
taxa_corrected2 = ConQuR(tax_tab=taxa, batchid=batchid, covariates=covar, batch_ref="1",
                         logistic_lasso=T, quantile_type="lasso", interplt=T)
taxa_corrected2[146:150, 1:5]

par(mfrow=c(2, 3))

##########check##################
pdf("Remove_batch_effect_byCountry.pdf") 
Plot_PCoA(TAX=taxa, factor=batchid, main="Before Correction, Bray-Curtis")
Plot_PCoA(TAX=taxa_corrected1, factor=batchid, main="ConQuR (Default), Bray-Curtis")
Plot_PCoA(TAX=taxa_corrected2, factor=batchid, main="ConQuR (Penalized), Bray-Curtis")

Plot_PCoA(TAX=taxa, factor=batchid, dissimilarity="Aitch", main="Before Correction, Aitchison")
Plot_PCoA(TAX=taxa_corrected1, factor=batchid, dissimilarity="Aitch", main="ConQuR (Default), Aitchison")
Plot_PCoA(TAX=taxa_corrected2, factor=batchid, dissimilarity="Aitch", main="ConQuR (Penalized), Aitchison")
dev.off()

library(coda.base)
PERMANOVA_R2(TAX=taxa, batchid=batchid, covariates=covar, key_index=1)
PERMANOVA_R2(TAX=taxa_corrected1, batchid=batchid, covariates=covar, key_index=1)
PERMANOVA_R2(TAX=taxa_corrected2, batchid=batchid, covariates=covar, key_index=1)

sbp = covar[, 'Group']
taxa_result = list(taxa, taxa_corrected1, taxa_corrected2)

pred_rmse = matrix(ncol=3, nrow=5)
colnames(pred_rmse) = c("Original", "ConQuR (Default)", "ConQuR (Penalized)")

for (ii in 1:3){
  pred_rmse[, ii] = RF_Pred_Regression(TAX=taxa_result[[ii]], variable=sbp)$rmse_across_fold
}

par(mfrow=c(1,1))
boxplot(pred_rmse, main="RMSE of Predicting Risk")

#write.csv(t(taxa_corrected2),"taxa_ALL_corrected2.csv")
#############using nornalized data to predict#########
taxa_corrected2_table<-melt(taxa_corrected2)
colnames(taxa_corrected2_table)<-c("index","genus","count")

sum<-aggregate(taxa_corrected2_table$count,by=list(taxa_corrected2_table$index),sum)

colnames(sum)<-c("index","sum")

table<-merge(taxa_corrected2_table,sum,by="index")

table$RA<-table$count/table$sum

to_compare<-tapply(table$RA,list(table$genus,table$index),median)
to_compare[is.na(to_compare)]<-0
write.csv(to_compare,"Risk_batcheffect_removal_relative_abundance.csv")

#############################Risk Cross-validation model###########################
meta.and.features <- read.lefse("Risk Model Data/108_trainSet_Risk.tsv",
                                rows.meta = 1:3, row.samples = 4)

meta <- meta.and.features$meta
feat <- meta.and.features$feat

siamcat  <- siamcat(feat=feat, meta=meta,
                    label='Risk', case='High')


#######60个samples tset##########
meta.and.features_test <- read.lefse("Risk Model Data/60_testSet_Risk.tsv",
                                     rows.meta = 1:3 ,row.samples = 4)

meta_test <- meta.and.features_test$meta
feat_test <- meta.and.features_test$feat

siamcat_test  <- siamcat(feat=feat_test, meta=meta_test,
                         label='Risk', case='High')

siamcat <- filter.features(
  siamcat,
  filter.method = 'abundance',
  cutoff = 0.0001,
  rm.unmapped = TRUE,
  verbose=2
)
siamcat <- normalize.features(
  siamcat,
  norm.method = "log.unit",
  norm.param=list(log.n0=1e-06,n.p=2, norm.margin=1),
  verbose = 2
)        
siamcat <-  create.data.split(
  siamcat,
  num.folds = 7,
  num.resample = 5
)
siamcat <- train.model(
  siamcat,
  method = "lasso"
)
siamcat <- make.predictions(siamcat)
siamcat <-  evaluate.predictions(siamcat)

##### prediction#######
siamcat_test <- normalize.features(siamcat_test,
                                   norm.param=norm_params(siamcat),
                                   feature.type='original',
                                   verbose = 2)

siamcat_test <- make.predictions(
  siamcat = siamcat,
  siamcat.holdout = siamcat_test,
  normalize.holdout =FALSE);

siamcat_test <- evaluate.predictions(siamcat_test)


model.evaluation.plot('108samples-train'=siamcat,
                      '60samples-public-test'=siamcat_test,
                      colours=c('orange', 'dodgerblue4'),
                      fn.plot = 'Risk_Cross-validation_model_2cohorts.pdf')



#########74 train samples and test on 2 cohorts############

meta.and.features <- read.lefse("Risk Model Data/LEfSeRisk_74_train.tsv",
                                rows.meta = 1:3, row.samples = 4)

meta <- meta.and.features$meta
feat <- meta.and.features$feat

siamcat  <- siamcat(feat=feat, meta=meta,
                    label='Risk', case='High')

#######34 samples test (cohort2)###########
#######60个samples test (Public data)###########
meta.and.features_test <- read.lefse("Risk Model Data/60_testSet_Risk.tsv",
                                     rows.meta = 1:3 ,row.samples = 4)

meta.and.features_test.T <- read.lefse("Risk Model Data/LEfSeRisk_34_test.tsv",
                                       rows.meta = 1:3, row.samples = 4)
meta_test <- meta.and.features_test$meta
feat_test <- meta.and.features_test$feat

siamcat_test  <- siamcat(feat=feat_test, meta=meta_test,
                         label='Risk', case='High')

meta_test.T <- meta.and.features_test.T$meta
feat_test.T <- meta.and.features_test.T$feat

siamcat_test.T  <- siamcat(feat=feat_test.T, meta=meta_test.T,
                           label='Risk', case='High')


siamcat <- filter.features(
  siamcat,
  filter.method = 'abundance',
  cutoff = 0.0001,
  rm.unmapped = TRUE,
  verbose=2
)
siamcat <- normalize.features(
  siamcat,
  norm.method = "log.unit",
  norm.param=list(log.n0=1e-06,n.p=2, norm.margin=1),
  verbose = 2
)        
siamcat <-  create.data.split(
  siamcat,
  num.folds = 7,
  num.resample = 5
)
siamcat <- train.model(
  siamcat,
  method = "lasso"
)
siamcat <- make.predictions(siamcat)
siamcat <-  evaluate.predictions(siamcat)

##### prediction#######
siamcat_test <- normalize.features(siamcat_test,
                                   norm.param=norm_params(siamcat),
                                   feature.type='original',
                                   verbose = 2)

siamcat_test.T <- normalize.features(siamcat_test.T,
                                     norm.param=norm_params(siamcat),
                                     feature.type='original',
                                     verbose = 2)

siamcat_test <- make.predictions(
  siamcat = siamcat,
  siamcat.holdout = siamcat_test,
  normalize.holdout =FALSE);

siamcat_test.T <- make.predictions(
  siamcat = siamcat,
  siamcat.holdout = siamcat_test.T,
  normalize.holdout = FALSE);

siamcat_test <- evaluate.predictions(siamcat_test)
siamcat_test.T <- evaluate.predictions(siamcat_test.T)

model.evaluation.plot('74samples-train'=siamcat,
                      '60samples-public-test'=siamcat_test,
                      '34samples-test'=siamcat_test.T,
                      colours=c('orange', 'dodgerblue4',"skyblue"),
                      fn.plot = 'Risk_Cross-validation model_3cohorts.pdf')

####################Liver metastasis based removal###################

L6data74<-read.csv("Metastasis Model Data/Before_ConQuR_74samples_count.csv",header=TRUE)
L6data34<-read.csv("Metastasis Model Data/Before_ConQuR_34samples_count.csv",header=TRUE)
batch<-read.csv("Metastasis Model Data/16Slist_union_batch2.csv",header=TRUE)


data74<-melt(L6data74)
data34<-melt(L6data34)

dataALL<-rbind(data74,data34)
colnames(dataALL)<-c("SampleID","Genus","Count")
dataALL<-dataALL[dataALL$Count>2,]
dataALL<-merge(batch,dataALL,by="SampleID",all=FALSE)
dataALL_matrix<-tapply(dataALL$Count,list(dataALL$SampleID,dataALL$Genus),median)

dataALL_matrix_WithGroup<-cbind(dataALL_matrix,batch,by="SampleID")
dataALL_matrix_WithGroup[is.na(dataALL_matrix_WithGroup)]<-0

taxa = dataALL_matrix_WithGroup[, 1:1190]
batchid = as.factor(dataALL_matrix_WithGroup[, 'batchid'])
covar = dataALL_matrix_WithGroup[, c("Risk","Gender","Age","Tumor_size","Mitotic_index","Adjuvant_imatinib","Necrosis","Liver_metastasis")]

options(warn=-1)
taxa_corrected1 = ConQuR(tax_tab=taxa, batchid=batchid, covariates=covar, batch_ref="1")

#taxa_corrected1[146:150, 1:5]

options(warn=-1)
taxa_corrected2 = ConQuR(tax_tab=taxa, batchid=batchid, covariates=covar, batch_ref="1",
                         logistic_lasso=T, quantile_type="lasso", interplt=T)
#taxa_corrected2[146:150, 1:5]

par(mfrow=c(2, 3))

##########check##################
pdf("Remove_batch_effect_34+74.pdf") 
 Plot_PCoA(TAX=taxa, factor=batchid, main="Before Correction, Bray-Curtis")
 Plot_PCoA(TAX=taxa_corrected1, factor=batchid, main="ConQuR (Default), Bray-Curtis")
 Plot_PCoA(TAX=taxa_corrected2, factor=batchid, main="ConQuR (Penalized), Bray-Curtis")

 Plot_PCoA(TAX=taxa, factor=batchid, dissimilarity="Aitch", main="Before Correction, Aitchison")
 Plot_PCoA(TAX=taxa_corrected1, factor=batchid, dissimilarity="Aitch", main="ConQuR (Default), Aitchison")
 Plot_PCoA(TAX=taxa_corrected2, factor=batchid, dissimilarity="Aitch", main="ConQuR (Penalized), Aitchison")
 dev.off()
 
library(coda.base)
PERMANOVA_R2(TAX=taxa, batchid=batchid, covariates=covar, key_index=8)
PERMANOVA_R2(TAX=taxa_corrected1, batchid=batchid, covariates=covar, key_index=8)
PERMANOVA_R2(TAX=taxa_corrected2, batchid=batchid, covariates=covar, key_index=8)

sbp = covar[, 'Liver_metastasis']
taxa_result = list(taxa, taxa_corrected1, taxa_corrected2)

pred_rmse = matrix(ncol=3, nrow=5)
colnames(pred_rmse) = c("Original", "ConQuR (Default)", "ConQuR (Penalized)")

for (ii in 1:3){
  pred_rmse[, ii] = RF_Pred_Regression(TAX=taxa_result[[ii]], variable=sbp)$rmse_across_fold
}

par(mfrow=c(1,1))
boxplot(pred_rmse, main="RMSE of Predicting Liver_metastasis")

#############using nornalized data to predict#########
taxa_corrected2_table<-melt(taxa_corrected2)
colnames(taxa_corrected2_table)<-c("index","genus","count")

sum<-aggregate(taxa_corrected2_table$count,by=list(taxa_corrected2_table$index),sum)

colnames(sum)<-c("index","sum")

table<-merge(taxa_corrected2_table,sum,by="index")

table$RA<-table$count/table$sum

to_compare<-tapply(table$RA,list(table$genus,table$index),median)
to_compare[is.na(to_compare)]<-0
write.csv(to_compare,"LM_batcheffect_removal_relative_abundance.csv")

#########non-AI group predicts AI group############
#########31-nonAI smaples train######
meta.and.features <- read.lefse("Metastasis Model Data/LEfSe_31_nonAI_all_v2.tsv",
                                rows.meta = 1:8, row.samples = 9)
meta <- meta.and.features$meta
feat <- meta.and.features$feat

siamcat  <- siamcat(feat=feat, meta=meta,
                    label='Liver_metastasis', case='Yes')


#######77-AI samples test###########
meta.and.features_test_AI <- read.lefse("Metastasis Model Data/LEfSe_77_AI_all_v2.tsv",
                                     rows.meta = 1:8, row.samples = 9)

meta_test_AI <- meta.and.features_test_AI$meta
feat_test_AI <- meta.and.features_test_AI$feat

siamcat_test_AI  <- siamcat(feat=feat_test_AI, meta=meta_test_AI,
                         label='Liver_metastasis', case='Yes')

siamcat <- filter.features(
  siamcat,
  filter.method = 'abundance',
  cutoff = 0.0001,
  rm.unmapped = TRUE,
  verbose=2
)
siamcat <- normalize.features(
  siamcat,
  norm.method = "log.unit",
  norm.param=list(log.n0=1e-05, n.p=2, norm.margin=1),
  #norm.param = list(log.n0 = 1e-05, sd.min.q = 0.1),
  verbose = 2
)        
siamcat <-  create.data.split(
  siamcat,
  num.folds = 5,
  num.resample = 3
)
siamcat <- train.model(
  siamcat,
  method = "lasso"
)
siamcat <- make.predictions(siamcat)
siamcat <-  evaluate.predictions(siamcat)

#####predict#######

siamcat_test_AI <- normalize.features(siamcat_test_AI,
                                   norm.param=norm_params(siamcat),
                                   feature.type='original',
                                   verbose = 2)

siamcat_test_AI <- make.predictions(
  siamcat = siamcat,
  siamcat.holdout = siamcat_test_AI,
  normalize.holdout = TRUE);

siamcat_test_AI <- evaluate.predictions(siamcat_test_AI)

model.evaluation.plot('31samples-nonAI-train'=siamcat,
                      '77samples-AI-test'=siamcat_test_AI,
                      colours=c("purple","darkgreen"),
                      fn.plot = 'nonAI_model_evaluation_selection_logunit.pdf')

siamcat <- check.associations(siamcat, log.n0 = 1e-06, alpha = 0.05)
association.plot(siamcat, sort.by = 'fc', fn.plot = 'nonAI_association_selection.pdf',
                 panels = c('fc', 'prevalence', 'auroc'))
####################nonAI feature selection model (Optional)#########################
library(SIAMCAT)
meta.and.features <- read.lefse("Metastasis Model Data/LEfSe_31_nonAI_all_v2.tsv",
                                rows.meta = 1:8, row.samples = 9)
meta <- meta.and.features$meta
feat <- meta.and.features$feat

label <- create.label(meta=meta, label="Liver_metastasis", case = "Yes")

siamcat <- siamcat(feat=feat, label=label, meta=meta)

sc.obj <- filter.features(siamcat,
                          filter.method = 'abundance',
                          cutoff = 0.001)

sc.obj <- check.associations(sc.obj, log.n0 = 1e-06, alpha = 0.05)
association.plot(sc.obj, sort.by = 'fc', fn.plot = '31association_plots_Live_metastasis.pdf',
                 panels = c('fc', 'prevalence', 'auroc'))

check.confounders(sc.obj, fn.plot = '31confounder_plotsLive_metastasis.pdf',
                  meta.in = NULL, feature.type = 'filtered')
sc.obj <- normalize.features(sc.obj, norm.method = "log.unit",
                             norm.param = list(log.n0 = 1e-05, n.p = 2,norm.margin = 1))
sc.obj <-  create.data.split(sc.obj, num.folds = 5, num.resample = 3)
sc.obj <- train.model(sc.obj, method = "lasso")
model_type(sc.obj)
models <- models(sc.obj)
models[[1]]$model
sc.obj <- make.predictions(sc.obj)
pred_matrix <- pred_matrix(sc.obj)
sc.obj <-  evaluate.predictions(sc.obj)
model.evaluation.plot(sc.obj,fn.plot = '31evaluationLive_metastasis.pdf')
#model.interpretation.plot(sc.obj, fn.plot = '31interpretationLive_metastasis.pdf',
#                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')


######################using cohort 1 to predict cohort 2########################
meta.and.features <- read.lefse("Metastasis Model Data/74_trainSet_LM.tsv",
                                rows.meta = 1:8, row.samples = 9)
meta <- meta.and.features$meta
feat <- meta.and.features$feat

siamcat  <- siamcat(feat=feat, meta=meta,
                    label='Liver_metastasis', case='Yes')

#######cohort2-test###########

meta.and.features_test <- read.lefse("Metastasis Model Data/34_testSet_LM.tsv",
                                     rows.meta = 1:8, row.samples = 9)

meta_test <- meta.and.features_test$meta
feat_test <- meta.and.features_test$feat


siamcat_test  <- siamcat(feat=feat_test, meta=meta_test,
                            label='Liver_metastasis', case='Yes')

siamcat <- filter.features(
  siamcat,
  filter.method = 'abundance',
  cutoff = 0.0001,
  rm.unmapped = TRUE,
  verbose=2
)
siamcat <- normalize.features(
  siamcat,
  norm.method = "log.unit",
  norm.param=list(log.n0=1e-05, n.p=2, norm.margin=1),
  #norm.param = list(log.n0 = 1e-05, sd.min.q = 0.1),
  verbose = 2
)        
siamcat <-  create.data.split(
  siamcat,
  num.folds = 5,
  num.resample = 3
)
siamcat <- train.model(
  siamcat,
  method = "lasso"
)
siamcat <- make.predictions(siamcat)
siamcat <-  evaluate.predictions(siamcat)

######预测#######
siamcat_test <- normalize.features(siamcat_test,
                                   norm.param=norm_params(siamcat),
                                   feature.type='original',
                                   verbose = 2)


siamcat_test <- make.predictions(
  siamcat = siamcat,
  siamcat.holdout = siamcat_test,
  normalize.holdout = TRUE);

siamcat_test <- evaluate.predictions(siamcat_test)


model.evaluation.plot('Cohort 1 train'=siamcat,
                      'Cohort 2 test'=siamcat_test,
                      colours=c("#3E5CC5", "#E64E00"),
                      fn.plot = 'LM_Cohort 1 predict cohort 2_model_.pdf')

################################feature selection model##############################
#######LM feature selection model########

meta.and.features <- read.lefse("Metastasis Model Data/LEfSe_74+34_0.1%.tsv",
                                rows.meta = 1:8, row.samples = 9)
meta <- meta.and.features$meta
feat <- meta.and.features$feat

label <- create.label(meta=meta, label="Liver_metastasis", case = "Yes")

siamcat <- siamcat(feat=feat, label=label, meta=meta)

sc.obj <- filter.features(siamcat,
                          filter.method = 'abundance',
                          cutoff = 0.001)

sc.obj <- check.associations(sc.obj, log.n0 = 1e-06, alpha = 0.05)
association.plot(sc.obj, sort.by = 'fc', fn.plot = 'LM_association_plots.pdf',
                 panels = c('fc', 'prevalence', 'auroc'))

check.confounders(sc.obj, fn.plot = 'LM_confounder_plots.pdf',
                  meta.in = NULL, feature.type = 'filtered')
sc.obj <- normalize.features(sc.obj, norm.method = "log.unit",
                             norm.param = list(log.n0 = 1e-05, n.p = 2,norm.margin = 1))
sc.obj <-  create.data.split(sc.obj, num.folds = 5, num.resample = 3)
sc.obj <- train.model(sc.obj, method = "lasso")
model_type(sc.obj)
models <- models(sc.obj)
models[[1]]$model
sc.obj <- make.predictions(sc.obj)
pred_matrix <- pred_matrix(sc.obj)
sc.obj <-  evaluate.predictions(sc.obj)
model.evaluation.plot(sc.obj,fn.plot = 'LM_evaluation.pdf')
model.interpretation.plot(sc.obj, fn.plot = 'LM_interpretation.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')

#####Risk feature selection model######
meta.and.features <- read.lefse("Risk Model Data/LEfSeRisk_168.tsv",
                                rows.meta = 1:3, row.samples = 4)
meta <- meta.and.features$meta
feat <- meta.and.features$feat

label <- create.label(meta=meta, label="Risk", case = "High")

siamcat <- siamcat(feat=feat, label=label, meta=meta)

sc.obj <- filter.features(siamcat,
                          filter.method = 'abundance',
                          cutoff = 0.0001)

sc.obj <- check.associations(sc.obj, log.n0 = 1e-05, alpha = 0.05)
association.plot(sc.obj, sort.by = 'fc', fn.plot = 'Risk_association_plots.pdf',
                 panels = c('fc', 'prevalence', 'auroc'))

check.confounders(sc.obj, fn.plot = 'LM_confounder_plots.pdf',
                  meta.in = NULL, feature.type = 'filtered')
sc.obj <- normalize.features(sc.obj, norm.method = "log.unit",
                             norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1))
sc.obj <-  create.data.split(sc.obj, num.folds = 7, num.resample = 5)
sc.obj <- train.model(sc.obj, method = "lasso")
model_type(sc.obj)
models <- models(sc.obj)
models[[1]]$model
sc.obj <- make.predictions(sc.obj)
pred_matrix <- pred_matrix(sc.obj)
sc.obj <-  evaluate.predictions(sc.obj)
model.evaluation.plot(sc.obj,fn.plot = 'Risk_evaluation.pdf')
model.interpretation.plot(sc.obj, fn.plot = 'Risk_interpretation.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')

