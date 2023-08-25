setwd("D:/GIST 16s/GIST microbiome")

library("survival")
library("survminer")
library(reshape2)
library(lefser)
library(ggplot2)
library(ade4)
library(RColorBrewer)
library(vegan)
library(rstatix)
library(ggpubr)
library(ggprism)
library(dplyr)
library(SIAMCAT)
library("plotrix")

###################病例信息统计##########################
Final_description_patient_info<-read.csv("Final_patient_info.csv",header=TRUE)
#####densityplot比较两组数据的肝转移发生时间#####

library(plyr)
library(ggpubr)
mu <- ddply(Final_description_patient_info, "Liver_metastasis", summarise, grp.mean=mean(Time))
head(mu)
densityplot_time<-ggplot(Final_description_patient_info, aes(x=Time,fill=Liver_metastasis)) + 
                  scale_fill_manual(values=c("#999999","#E69F00"))+
                  geom_density(alpha=0.4)+
                  #geom_vline(data=mu, aes(xintercept=grp.mean, color=Liver_metastasis),linetype="dashed")+
                  scale_color_manual(values=c( "#999999","#E69F00"))+
                  theme_classic()+
                  facet_wrap(~Batch,ncol=1)

ggsave(filename = "densityplot_time2.pdf", plot =densityplot_time, width = 6, height = 6)

#######pcoa展示数据分布######
####肝转移#####
library(vegan)
library(ape)
matrix_union<-read.csv("taxa_corrected2_Union_matrix.csv",header=TRUE)
rownames(matrix_union)<-matrix_union$X
matrix_union<-matrix_union[,-1]

group_union<-Final_description_patient_info
#rownames(group_union)<-group_union$Specimen_number

otu <- t(matrix_union)
otu.distance <- vegdist(otu)
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)

pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)

p <- ggplot(pc12,aes(x=V1, y=V2))+
  geom_point(size=3)+
  theme_bw()
p

colnames(group_union)<-c("samples","Gender","Age","Tumor_size","Mitotic_index","Risk", "Adjuvant_imatinib", "Gene_test","PLR","NLR","FIB","PNI","Presence_of.necrosis","Liver_metastasis","Time", "Batch")
df <- merge(pc12,group_union,by="samples")

color=c("#999999","#E69F00")
dune.div <- adonis2(otu ~ Liver_metastasis, data = unique(df[c(1,16)]), permutations = 999, method="bray")

dune_adonis <- paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)


p1<-ggplot(data=df,aes(x=V1,y=V2,
                       color=Liver_metastasis,shape=Batch))+
  theme_bw()+
  geom_point(size=1.8)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  #geom_text(aes(label=samples, y=V2+0.02,x=V1+0.02, vjust=0),size=3.5)+
  guides(color=guide_legend(title="Liver_metastasis"))+
  labs(x=paste0("PC1 ",pc[1],"%"),
       y=paste0("PC2 ",pc[2],"%"),title=dune_adonis)+
  scale_color_manual(values = color) +
  scale_fill_manual(values = c("#999999","#E69F00"))+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())
  #stat_ellipse(data=df,
  #             #geom = "polygon",
  #             level=0.95,
  #             linetype = 2,
  #             size=0.5,
  #             aes(fill=Liver_metastasis),
  #             alpha=0.2)


ggsave(p1, file="Live.metastasis_pcoa.pdf", width=6, height=5,limitsize = FALSE)
#ggsave(p1, file="Adjuvant_imatinib_pcoa.pdf", width=6, height=5,limitsize = FALSE)

#####风险评估#####
group_all<-read.csv("16Slist_all_batch.csv",header=TRUE)
matrix_all<-read.csv("taxa_corrected2_ALL_B1_matrix_2group.csv",header=TRUE)
rownames(matrix_all)<-matrix_all$X
matrix_all<-matrix_all[,-1]

otu <- t(matrix_all)
otu.distance <- vegdist(otu)
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)

pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)

p <- ggplot(pc12,aes(x=V1, y=V2))+
  geom_point(size=3)+theme_bw()
p

colnames(group_all)<-c("samples","batch","Gender","Age","Risk","Country", "batchid")
df <- merge(pc12,group_all,by="samples")

color=c("#E64E00","#65B48E")
dune.div <- adonis2(otu ~ Risk, data = unique(df[c(1,7)]), permutations = 999, method="bray")

dune_adonis <- paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)


p1<-ggplot(data=df,aes(x=V1,y=V2,
                       color=Risk,shape=Country))+
  theme_bw()+
  geom_point(size=1.8)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  #geom_text(aes(label=samples, y=V2+0.02,x=V1+0.02, vjust=0),size=3.5)+
  guides(color=guide_legend(title="Liver metastasis"))+
  labs(x=paste0("PC1 ",pc[1],"%"),
       y=paste0("PC2 ",pc[2],"%"),title=dune_adonis)+
  scale_color_manual(values = color) +
  scale_fill_manual(values = c("#E64E00","#65B48E"))+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())
#stat_ellipse(data=df,
#             #geom = "polygon",
#             level=0.95,
#             linetype = 2,
#             size=0.5,
#             aes(fill=Liver_metastasis),
#             alpha=0.2)


ggsave(p1, file="Risk_pcoa.pdf", width=6, height=5,limitsize = FALSE)

##############shanno diversity#############

#####肝转移分组#####
sample_id<-final_patient_cox$Specimen_number
group<-group_union$Liver_metastasis
#matrix_union<-read.csv("taxa_corrected2_Union_matrix.csv",header=TRUE)
union_list_9.30_0.0001<-union_list_9.30[union_list_9.30$RA>0.0001,]
union_list_9.30_0.0001_inhibit<-union_list_9.30_0.0001[union_list_9.30_0.0001$group_LM=="inhibitLM",]
union_list_9.30_0.0001_promote<-union_list_9.30_0.0001[union_list_9.30_0.0001$group_LM=="promoteLM",]


matrix_union2<-tapply(union_list_9.30_0.0001_promote$RA,list(union_list_9.30_0.0001_promote$genus,union_list_9.30_0.0001_promote$samples),median)
matrix_union2[is.na(matrix_union2)]<-0
#matrix_withGroup<-melt(matrix_union)
#rownames(matrix_union)<-matrix_union$X
#matrix_union<-matrix_union[,-1]
data_norm=t(matrix_union2)
data_group=data.frame(sample_id, group)

data_norm_shannon=diversity(data_norm, "shannon")
data_ggplot=data.frame(data_norm_shannon)
data_ggplot=data.frame(data_ggplot, data_group["group"])

Liver_boxplot<-ggplot(data_ggplot, aes(x=group, y=data_norm_shannon, fill=group))+
  geom_boxplot(outlier.shape=16,
               outlier.size=2)+
  scale_fill_manual(values = c("#999999","#E69F00"))+
  labs(title="Alpha diversity", x="Liver metastasis", y="Shannon Index")+
  #geom_jitter(color="black", size=1, alpha=0.9)+
  theme(plot.title=element_text(hjust=0.5), legend.title=element_blank())
####p-value#######
#p_adj <- get_adj_p(data_ggplot, .col = "data_norm_shannon", .grp = "group",method = "wilcox.test", p.adjust.method = "BH")
#p_adj

my_comparisons <- list( c("Yes","No") )
Liver_boxplot<-Liver_boxplot+
  #stat_compare_means()+
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+
  theme_bw()

ggsave(Liver_boxplot, file="Liver_shannon.pdf", width=3, height=5)

#####Risk分组#####
sample_id<-group_all$samples
group<-group_all$Risk
#matrix_union<-read.csv("taxa_corrected2_Union_matrix.csv",header=TRUE)
#matrix_withGroup<-melt(matrix_union)

data_norm=t(matrix_all)
data_group=data.frame(sample_id, group)

data_norm_shannon=diversity(data_norm, "shannon")
data_ggplot=data.frame(data_norm_shannon)
data_ggplot=data.frame(data_ggplot, data_group["group"])

Liver_boxplot<-ggplot(data_ggplot, aes(x=group, y=data_norm_shannon, fill=group))+
  geom_boxplot(outlier.shape=16,
               outlier.size=2)+
  scale_fill_manual(values = c("#E64E00","#65B48E"))+
  labs(title="Alpha diversity", x="Risk", y="Shannon Index")+
  #geom_jitter(color="black", size=1, alpha=0.9)+
  #facet_wrap(~ Country)+
  theme(plot.title=element_text(hjust=0.5), legend.title=element_blank())
####p-value#######
#p_adj <- get_adj_p(data_ggplot, .col = "data_norm_shannon", .grp = "group",method = "wilcox.test", p.adjust.method = "BH")
#p_adj

my_comparisons <- list( c("High","Other") )
Liver_boxplot<-Liver_boxplot+
  #stat_compare_means()+
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+
  theme_bw()

ggsave(Liver_boxplot, file="Risk_shannon.pdf", width=3, height=5)
################################临床标量统计####################
Final_description_patient_info$count<-1
Gender_count<-tapply(Final_description_patient_info$count,list(Final_description_patient_info$Liver_metastasis,Final_description_patient_info$Gender),sum)
Age_count<-tapply(Final_description_patient_info$count,list(Final_description_patient_info$Liver_metastasis,Final_description_patient_info$Age),sum)
Tumor_size_count<-tapply(Final_description_patient_info$count,list(Final_description_patient_info$Liver_metastasis,Final_description_patient_info$Tumor_size),sum)
Mitotic_index_count<-tapply(Final_description_patient_info$count,list(Final_description_patient_info$Liver_metastasis,Final_description_patient_info$Mitotic_index),sum)
Risk_count<-tapply(Final_description_patient_info$count,list(Final_description_patient_info$Liver_metastasis,Final_description_patient_info$Risk),sum)
Adjuvant_imatinib_count<-tapply(Final_description_patient_info$count,list(Final_description_patient_info$Liver_metastasis,Final_description_patient_info$Adjuvant_imatinib),sum)
Gene_test_count<-tapply(Final_description_patient_info$count,list(Final_description_patient_info$Liver_metastasis,Final_description_patient_info$Gene_test),sum)
Presence_of.necrosis_count<-tapply(Final_description_patient_info$count,list(Final_description_patient_info$Liver_metastasis,Final_description_patient_info$Presence_of.necrosis),sum)

table_bind<-cbind(Gender_count,Age_count,Tumor_size_count,Mitotic_index_count,Risk_count,Adjuvant_imatinib_count,Gene_test_count,Presence_of.necrosis_count)
table_bind<-as.data.frame(table_bind)
write.csv(table_bind,"table_coundition_sum.csv")

####在相对丰度表里筛选出26个 feature#############
####计算correlation###
select_genus<-read.csv("ALL_features_true.csv",header=TRUE)
matrix_union<-read.csv("union108_genus_all_matrix.csv",header=TRUE)
select_ab<-merge(select_genus,matrix_union,by="X",all=FALSE)
select_ab_list<-melt(select_ab)
select_ab_list$value<-select_ab_list$value*100
select_ab_list<-select_ab_list[select_ab_list$value>0.01,]
colnames(select_ab_list)<-c("genus","Specimen_number","RA")

select_table<-merge(select_ab_list,Final_description_patient_info,by="Specimen_number",all=FALSE)

select_table$Risk[which(select_table$Risk =='High')] <- '3'
select_table$Risk[which(select_table$Risk =='Low')] <- '2'
select_table$Risk[which(select_table$Risk =='Intermediate')] <- '2'

select_table$Adjuvant_imatinib[which(select_table$Adjuvant_imatinib =='Yes')] <- '2'
select_table$Adjuvant_imatinib[which(select_table$Adjuvant_imatinib =='No')] <- '1'

select_table$Liver_metastasis[which(select_table$Liver_metastasis =='Yes')] <- '2'
select_table$Liver_metastasis[which(select_table$Liver_metastasis =='No')] <- '1'

select_table_1<-select_table[select_table$genus=="Staphylococcus",]
res <- cor.test(select_table_1$RA, as.numeric(select_table_1$Risk), 
                method = "spearman")
res2 <- cor.test(as.numeric(select_table_1$Adjuvant_imatinib), as.numeric(select_table_1$Liver_metastasis), 
                method = "spearman")
res
res2

#write.csv(ca,"surval_influencedby_im.csv")
#####画气泡图#####
correlation_plot<-read.csv("spearman_correlation_genus.csv",header=TRUE)

pp <- ggplot(correlation_plot,aes(x=rho,y=genus))+
      geom_point(aes(color=-log10(P.value),size=propotion))+
      scale_colour_gradient(low="grey90",high="red")+
      facet_wrap( ~ Condition)+
      #geom_vline(xintercept=1.30103)+
      theme_bw()

ggsave(filename = "correlation_plot.pdf", plot = pp, width = 6, height = 4)

#######################cox回归生存分析####################
Final_patient_info<-read.csv("Final_patient_info_numeric.csv",header=TRUE)

#####批量单因素分析######
covariates <- c("Age","Gender","Tumor_size","Mitotic_index","Adjuvant_imatinib","Risk","Presence_of_necrosis","Gene_test")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(Time, Liver_metastasis)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = Final_patient_info)})
univ_models
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)","wald.test", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
class(univ_results)
res <- t(as.data.frame(univ_results, check.names = FALSE))
w<-as.data.frame(res)
write.csv(w,"Surv.csv")

#####单因素分析######
fit<- survfit(Surv(Time, Liver_metastasis) ~ Age, data = Final_patient_info)
fit<- survfit(Surv(Time, Liver_metastasis) ~ Gender, data = Final_patient_info)
fit<- survfit(Surv(Time, Liver_metastasis) ~ Tumor_size, data = Final_patient_info)
fit<- survfit(Surv(Time, Liver_metastasis) ~ Mitotic_index, data = Final_patient_info)
fit<- survfit(Surv(Time, Liver_metastasis) ~ Adjuvant_imatinib, data = Final_patient_info)
fit<- survfit(Surv(Time, Liver_metastasis) ~ Risk, data = Final_patient_info)
fit<- survfit(Surv(Time, Liver_metastasis) ~ Presence_of_necrosis, data = Final_patient_info)
fit<- survfit(Surv(Time, Liver_metastasis) ~ Gene_test, data = Final_patient_info)
######高危和中危的病人伊马提尼对肝转移的影响####
Final_patient_info_High<-Final_patient_info[Final_patient_info$Risk==c(2,3),]
fit<- survfit(Surv(Time, Liver_metastasis) ~ Adjuvant_imatinib, data = Final_patient_info_High)

c<-ggsurvplot(fit, data = Final_patient_info_High,
           surv.median.line = "hv", # Add medians survival
           # Change legends: title & labels
           legend.title = "Adjuvant_imatinib",
           legend.labs=c("No","Yes"),
           # Add p-value and tervals
           pval = TRUE,
           conf.int = TRUE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF","#5DC863),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#E7B800", "#2E9FDF"),
           #,"#5DC863"
           ggtheme = theme_bw() 
)

ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)}

g_to_save <- ggsave_workaround(c)

ggsave(filename = "survplot_HighandInter_Adjuvant_imatinib.pdf", plot = g_to_save, width = 6, height = 6)

##############多因素cox 回归分析######################
library("survival")
library("survminer")
fit<- survfit(Surv(Time, Liver_metastasis) ~ Risk+Adjuvant_imatinib, data = Final_patient_info)
summary(fit)
f<-coxph(Surv(Time, Liver_metastasis) ~ Risk+Adjuvant_imatinib+Gender+Age+Tumor_size+Mitotic_index+Presence_of_necrosis+Empedobacter+Shewanella, data = final_patient_cox)
#f<-coxph(Surv(Time, Liver_metastasis) ~ Risk+Gender+Age+Tumor_size+Mitotic_index+Gene_test+Presence_of_necrosis, data = final_patient_cox)
#sum.surv<-summary(f)
#c_index<-sum.surv$concordance

#ggsurvplot(fit, 
           #color = c("#E7B800", "#2E9FDF","#5DC863"),
#           ggtheme = theme_minimal())

HR<-ggforest(f, final_patient_cox)
ggsave(filename = "survplot_2marker_condition.pdf", plot = HR, width = 6, height = 4)
##########加入基因的cox回归分析##############
#union_list_9.30_2<-union_list_9.30[union_list_9.30$RA>0.001,]
#union_matrix<-tapply(union_list_9.30$RA*100,list(union_list_9.30$samples,union_list_9.30$genus),median)
#union_matrix[is.na(union_matrix)]<-0
#write.csv(union_matrix,"21genus_matrix_0.001.csv")

final_patient_cox<-read.csv("Final_patient_info_cox0.0001.csv",header=TRUE)
#f<-coxph(Surv(Time, Liver_metastasis) ~ Risk+Brevundimonas+Coprococcus+Chryseobacterium+Dialister+Dietzia+Dorea+Empedobacter+Escherichia.Shigella+Klebsiella+Massilia+Noviherbaspirillum+Pseudomonas+Shewanella+Sphingobacterium+Staphylococcus+Stenotrophomonas+Subdoligranulum+Subgroup_7+TRA3.20+Turicibacter, data = final_patient_cox)
f<-coxph(Surv(Time, Liver_metastasis) ~ Risk+
                                        Adjuvant_imatinib
                                        #+Gender+Age
                                        #+Gene_test
                                        #+Tumor_size+Mitotic_index+Tumor_site
                                        #+Presence_of_necrosis
                                        #+Empedobacter
                                        #+Shewanella
                                        #+Escherichia.Shigella
                                        #+Listeria
                                        #+Brevundimonas+Dialister+Dietzia+Dorea+Escherichia.Shigella+Klebsiella+Massilia+Noviherbaspirillum+Pseudomonas+Sphingobacterium+Staphylococcus+Stenotrophomonas+Subgroup_7+TRA3.20+Turicibacter
                                        #+shannon_21genus
                                         +shannon_promote
                                         +shannon_inhibit
                                        ,data = final_patient_cox)

HR_C<-ggforest(f,final_patient_cox)

ggsave(filename = "survplot_Risk+shannon_nonrisk.pdf", plot = HR_C, width = 7, height = 6)

###########c index#############
cI<-read.csv("C_index.csv",header=TRUE)
cI <- cI[order(cI$Concordance_Index), ]
cI$Sample <- factor(cI$Sample,levels=cI$Sample)

c_index<-ggplot(cI, aes(x=Sample, y=Concordance_Index, group=1)) + 
  #geom_point(color="red", size=3)
  geom_line(color="red")+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x =element_text(angle =60))
ggsave(filename = "c_index.pdf", plot = c_index, width = 4, height = 4)
###############nomogram##################
#library(glmnet)
library(rms)
library(VIM)
library(survival)
library(ggfortify)
library(rpart)

final_patient_cox<-read.csv("Final_patient_info_cox0.0001.csv",header=TRUE)
f2<-aareg(Surv(Time, Liver_metastasis) ~ Risk
         +Adjuvant_imatinib
         #+Gender+Age
         #+Gene_test
         #+Tumor_size+Mitotic_index+Tumor_site
         #+Presence_of_necrosis
         #+Brevundimonas+Coprococcus+Chryseobacterium+Dialister+Dietzia+Dorea+Empedobacter+Escherichia.Shigella+Klebsiella+Massilia+Noviherbaspirillum+Pseudomonas+Shewanella+Sphingobacterium+Staphylococcus+Stenotrophomonas+Subdoligranulum+Subgroup_7+TRA3.20+Turicibacter
         #+Empedobacter
         #+Shewanella
         #+Brevundimonas+Dialister+Dietzia+Dorea+Escherichia.Shigella+Klebsiella+Massilia+Noviherbaspirillum+Pseudomonas+Sphingobacterium+Staphylococcus+Stenotrophomonas+Subgroup_7+TRA3.20+Turicibacter
         #+shannon_21genus
         +shannon_promote
         +shannon_inhibit
         #+simpson_21genus
          ,data = final_patient_cox,dfbeta=TRUE)
#cox.zph(f2)
summary(f2)
a<-autoplot(f2)+theme_bw()
ggsave(filename = "covariates change over time.pdf", plot = a, width = 8, height = 5)

f1<-psm(Surv(Time, Liver_metastasis) ~ Risk
          +Adjuvant_imatinib
          #+Empedobacter
          #+Shewanella
          +shannon_promote
          +shannon_inhibit
          #+shannon_21genus
          ,data = final_patient_cox)
surv <- Survival(f1) # 建立生存函数

surv1 <- function(x)surv(1*24,lp=x)
surv2 <- function(x)surv(1*36,lp=x) 
surv3 <- function(x)surv(1*48,lp=x)
surv4 <- function(x)surv(1*66,lp=x)

dd<-datadist(final_patient_cox) 
options(datadist='dd') 

pdf("nomogram.pdf",width=8, height=6) 
b<-plot(nomogram(f1,
              fun=list(surv1,surv2,surv3,surv4),
              lp= F,
              funlabel=c('24-Month survival','36-Month survival','48-Month survival','66-Month survival'),
              maxscale=100,
              fun.at=c('1','0.9','0.85','0.80','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0')
              ),
     xfrac=.45)
#maxscale 参数指定最高分数，一般设置为100或者10分
#fun.at 设置生存率的刻度
#xfrac 设置数值轴与最左边标签的距离，可以调节下数值观察下图片变化情况
ggsave(filename = "nomogram.pdf", plot = b, width = 8, height = 5)
dev.off()

####################原始数据筛选########################

union108samples<-read.csv("taxa_corrected2_Union_matrix.csv",header=TRUE)
union108samples_list<-melt(union108samples)
union108samples_list<-union108samples_list[union108samples_list$value>0.001,]
union108samples_matrix<-tapply(union108samples_list$value,list(union108samples_list$X,union108samples_list$variable),median)
union108samples_matrix[is.na(union108samples_matrix)]<-0
write.csv(union108samples_matrix,"union108samples_matrix.csv")

public60samples<-read.csv("LEfSeRisk_public.csv",header=TRUE)
public60samples_list<-melt(public60samples)
public60samples_list<-public60samples_list[public60samples_list$value>0.001,]
public60samples_matrix<-tapply(public60samples_list$value,list(public60samples_list$ID,public60samples_list$variable),median)
public60samples_matrix[is.na(public60samples_matrix)]<-0
write.csv(public60samples_matrix,"public60samples_matrix.csv")

#######################肝转移分组#############################

library(SIAMCAT)
#meta.and.features <- read.lefse("LEfSe_74+34.tsv",
#                                rows.meta = 1:8, row.samples = 9)
meta.and.features <- read.lefse("LEfSe_LM_108.tsv",
                                rows.meta = 1:8, row.samples = 9)
meta <- meta.and.features$meta
feat <- meta.and.features$feat

label <- create.label(meta=meta, label="Liver_metastasis", case = "Yes")

siamcat <- siamcat(feat=feat, label=label, meta=meta)

sc.obj <- filter.features(siamcat,
                          filter.method = 'abundance',
                          cutoff = 0.001)

sc.obj <- check.associations(sc.obj, log.n0 = 1e-05, alpha = 0.05)
association.plot(sc.obj, sort.by = 'fc', fn.plot = 'V2_108association_plots_Live_metastasis_cotoff0.1.pdf',
                 panels = c('fc', 'prevalence', 'auroc'))

check.confounders(sc.obj, fn.plot = 'V2_108confounder_plotsLive_metastasis_cotoff0.1.pdf',
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
model.evaluation.plot(sc.obj,fn.plot = 'V2_108evaluationLive_metastasis_cotoff0.1.pdf')
model.interpretation.plot(sc.obj, fn.plot = 'V2_108interpretationLive_metastasis_cotoff0.1.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')

##################median relative feature weight################################
library(effsize)
feature_weights<-feature_weights(sc.obj)
feature_weights<-feature_weights[order(feature_weights$median.rel.weight),]
#write.csv(feature_weights,"feature_weights0.922.csv")
feature_weights_1<-feature_weights[feature_weights$median.rel.weight>0,]
feature_weights_1$group<-"Positive"
feature_weights_2<-feature_weights[feature_weights$median.rel.weight<0,]
feature_weights_2$group<-"Negative"
feature_weights<-rbind(feature_weights_2,feature_weights_1)
write.csv(feature_weights,"median_weight_feature_weights0.930_fullname.csv")
feature_weights<-read.csv("median_weight_feature_weights0.930.csv",header=TRUE)
#feature_weights_table<-melt(feature_weights[c(1:7)])
feature_weights <- feature_weights[order(feature_weights$median.rel.weight), ]

feature_weights$genus <- factor(feature_weights$genus,levels=feature_weights$genus)
feature_weights$percentage<-round(feature_weights$percentage,2)

##################################cliff delta################################
cliffD<-read.csv("selected_genus9.30cliffD_AIGroupP.csv",header=TRUE)
#cliffD<-read.csv("ALL_features_cliff.csv",header=TRUE)
#cliffD$group<-
cliffD2<-cliffD[cliffD$group=="AI",]
D2<-aggregate(cliffD2$CliffDelta,list(cliffD2$genus),mean)
cliffD2 <- D2[order(-D2$x),]
cliffD$genus <- factor(cliffD$genus,levels=cliffD2$Group.1)
#cliffD$group <-factor(cliffD$group,ordered=TRUE,levels=unique(cliffD2$group))

p<-ggplot(data = cliffD, mapping = aes(x =genus, y = CliffDelta,color=group))+ 
  geom_point(stat='identity',aes(color=group,size=Prevelence))+
  scale_color_manual(values=c("#FF0000", "#663366"))+
  geom_hline(aes(yintercept=0),colour="#999999",linetype=1)+
  geom_hline(aes(yintercept=0.33),colour="#999999",linetype="dashed")+
  geom_hline(aes(yintercept=-0.33),colour="#999999",linetype="dashed")+
  #scale_color_brewer(palette="Set1")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  #geom_text(color="white", size=2) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggsave(p, filename="ALLfeature_weights_CilffD_AI_dot.pdf", height=4, width=7)

##################  cliff.delta effsize  ###############
#LM matrix##
matrix_union<-read.csv("union108_genus_all_matrix.csv",header=TRUE)
#Risk matrix##
#matrix_union<-read.csv("all168samples_matrix_genus.csv",header=TRUE)
union_list<-melt(matrix_union)
#union_list<-union_list[union_list$value>0.0001,]
gene_list<-read.csv("selected_genus9.30.csv",header=TRUE)
colnames(gene_list)<-c("genus","group_LM")
#gene_list<-read.csv("compare_enrichment_socre_2group.csv",header=TRUE)
#gene_list<-read.csv("ALL_features_true.csv",header=TRUE)

colnames(union_list)<-c("genus","samples","RA")

union_list_9.30<-merge(gene_list,union_list,by="genus",all=FALSE)
#union_list_9.30$RA<-union_list_9.30$RA*100
group_union<-read.csv("sample_table_risk.csv",header=TRUE)
colnames(group_union)<-c("samples","Gender","Age","Tumor_size","Mitotic_index","Risk", "Adjuvant_imatinib", "Gene_test","PLR","NLR","FIB","PNI","Presence_of.necrosis","Liver_metastasis","Time", "Batch","group")
union_table<-merge(group_union,union_list_9.30,by="samples",all=FALSE)
#union_table<-union_table[union_table$RA>0,]

library(effsize)
#union_table.0<-union_table[union_table$Risk=="High",]
#union_table.1<-union_table[union_table$Risk==c("Low","Intermediate"),]

#union_table.2<-union_table.0[union_table.0$genus=="Nonomuraea",]
#union_table.3<-union_table.0[union_table.1$genus=="Nonomuraea",]

#cliff.delta(union_table.2$RA,union_table.3$RA,conf.level=.95)

union_table.0<-union_table[union_table$Adjuvant_imatinib=="Yes",]
union_table.1<-union_table[union_table$Adjuvant_imatinib=="No",]

union_table.2<-union_table.0[union_table.0$Liver_metastasis=="Yes",]
union_table.3<-union_table.0[union_table.0$Liver_metastasis=="No",]

union_table.4<-union_table.2[union_table.2$genus=="Coprococcus",]
union_table.5<-union_table.3[union_table.3$genus=="Coprococcus",]

cliff.delta(union_table.4$RA,union_table.5$RA,conf.level=.95)


#########################Risk分组###############################################

#####协和数据##########
meta.and.featuresR <- read.lefse("LEfSeRisk_168.tsv",
                                rows.meta = 1:3, row.samples = 4)

#meta.and.featuresR <- read.lefse("LEfSeRisk_168_0.1.tsv",
#                                 rows.meta = 1:3, row.samples = 4)

metaR <- meta.and.featuresR$meta
featR <- meta.and.featuresR$feat

labelR <- create.label(meta=metaR, label="Risk", case = "High")

siamcatR <- siamcat(feat=featR, label=labelR, meta=metaR)

sc.objR <- filter.features(siamcatR,
                          filter.method = 'abundance',
                          #rm.unmapped = TRUE,
                          cutoff = 0.0001)

sc.objR <- check.associations(sc.objR, log.n0 = 1e-05, alpha = 0.1)
association.plot(sc.objR, sort.by = 'fc', fn.plot = 'V2_168association_plots_Risk_cotoff0.1.pdf', max.show = 1000,
                 panels = c('fc', 'prevalence', 'auroc'))

check.confounders(sc.objR, fn.plot = 'V2_168confounder_plots_Risk_cotoff0.1.pdf',
                  meta.in = NULL, feature.type = 'filtered')
sc.objR <- normalize.features(sc.objR, norm.method = "log.unit",
                             norm.param = list(log.n0 = 1e-05, n.p = 2,norm.margin = 1))
sc.objR <-  create.data.split(sc.objR, num.folds = 5, num.resample = 3)
sc.objR <- train.model(sc.objR, method = "lasso")
model_type(sc.objR)
modelsR <- models(sc.objR)
modelsR[[1]]$model
sc.objR <- make.predictions(sc.objR)
pred_matrixR <- pred_matrix(sc.objR)
sc.objR <-  evaluate.predictions(sc.objR)
model.evaluation.plot(sc.objR,fn.plot = 'V2_168evaluation_Risk_cotoff0.1.PR.pdf',verbose=2)
model.interpretation.plot(sc.objR, fn.plot = 'V2_168interpretation_Risk_cotoff0.1.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type='zscore')

###############Risk calculation###############
feature_weightsR<-feature_weights(sc.objR)
feature_weightsR<-feature_weightsR[order(feature_weightsR$median.rel.weight),]
#write.csv(feature_weightsR,"feature_weights0.936.csv")
feature_weights_1R<-feature_weightsR[feature_weightsR$median.rel.weight>0,]
feature_weights_1R$group<-"Positive"
feature_weights_2R<-feature_weightsR[feature_weightsR$median.rel.weight<0,]
feature_weights_2R$group<-"Negative"
feature_weightsR<-rbind(feature_weights_2R,feature_weights_1R)
#write.csv(feature_weightsR,"median_weight_feature_weights0.936_fullname.csv")
##################venn plot############
lefse_matrix<-read.csv("LEfSeRisk_168_matrix.csv",header=TRUE)
lefse_list<-melt(lefse_matrix)
lefse_list<-lefse_list[lefse_list$value>0.001,]
lefse_matrix<-tapply(lefse_list$value,list(lefse_list$ID,lefse_list$variable),median)
lefse_matrix[is.na(lefse_matrix)]<-0
write.csv(lefse_matrix,"lefse_matrix_risk.csv")

#######公共数据(abandon)########
meta.and.features <- read.lefse("LEfSeRisk_public.tsv",
                                rows.meta = 1:3, row.samples = 4)
meta <- meta.and.features$meta
feat <- meta.and.features$feat

label <- create.label(meta=meta, label="Risk", case = "High")

siamcat <- siamcat(feat=feat, label=label, meta=meta)

sc.obj <- filter.features(siamcat,
                          filter.method = 'abundance',
                          cutoff = 0.0001)

sc.obj <- check.associations(sc.obj, log.n0 = 1e-05, alpha = 0.1)
association.plot(sc.obj, sort.by = 'fc', fn.plot = 'V2_Public_association_plots_Risk_public_cotoff0.1.pdf',
                 panels = c('fc', 'prevalence', 'auroc'))

check.confounders(sc.obj, fn.plot = 'V2_Public_confounder_plots_Risk_public_cotoff0.1.pdf',
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
model.evaluation.plot(sc.obj,fn.plot = 'V2_Public_evaluation_Risk_public_cotoff0.1.pdf')
model.interpretation.plot(sc.obj, fn.plot = 'V2_Public_interpretation_Risk_public_cotoff0.1.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')








##############heatmap for selected genus union important genus#####################
library(ComplexHeatmap)

genus_matrix<-read.csv("union108_genus_all_matrix.csv",header=TRUE)

#genus_list<-read.csv("selected_genus.csv",header=TRUE)
genus_list<-read.csv("selected_abundance_genus.csv",header=TRUE)
#genus_list<-read.csv("selected_vary_important_genus.csv",header=TRUE)

genus_matrix<-merge(genus_list,genus_matrix,by="X",all=FALSE)
#rownames(genus_matrix)<-genus_matrix$X
#genus_matrix<-genus_matrix[,-1]
genus_all<-melt(genus_matrix)
genus_all$value<-genus_all$value*100

colnames(genus_all)<-c("Genus","SampleID","Abundance")

groupinfo<-read.csv("16Slist_union_batch.csv",header=TRUE)

genus_group<-merge(groupinfo,genus_all,by="SampleID",all=FALSE)
genus_group<-genus_group[genus_group$Abundance>0.01,]
 
rownames(genus_matrix)<-genus_matrix$X
genus_matrix<-genus_matrix[,-1]

group<-groupinfo[c(1,7,10)]
rownames(group)<-group$SampleID
group<-group[,-1]


p<-pheatmap(log2(genus_matrix+0.00001),annotation_col=group,cluster_rows =TRUE,
            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
            border_color = "grey",
            cellwidth = 10, cellheight = 10,
            gaps_col = 73,
            cluster_cols = FALSE)    

p

ggsave(p, filename="heatmap_20_metagenomic_features_genus.pdf", height=10, width=20)
  
#####boxplot#####

b <- ggplot(genus_group, aes(x=Liver_metastasis, y=log2(Abundance),fill=Liver_metastasis)) + 
   facet_grid(Risk ~ Genus)+
  scale_fill_manual(values=c("#999999","#E69F00"))+
  #facet_wrap(~ Genus,scales = "free_y")+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5 ,fill = "black")+
  stat_compare_means()+
  theme_bw()
  
b
ggsave(b, filename="boxplot_12_abundance_genus.pdf", height=8, width=23)

####plot Shewanella##########
genus_matrix<-read.csv("union108_genus_all_matrix.csv",header=TRUE)
groupinfo<-read.csv("16Slist_union_batch.csv",header=TRUE)
genus_all<-melt(genus_matrix)
genus_all$value<-genus_all$value*100
colnames(genus_all)<-c("Genus","SampleID","Abundance")
genus_group<-merge(groupinfo,genus_all,by="SampleID",all=FALSE)

genus_group<-genus_group[genus_group$Abundance>0.01,]
genus_Shewanella<-genus_group[genus_group$Genus=="Brevundimonas",]

plot_Shewanella<-ggplot(genus_Shewanella, aes(x=Liver_metastasis, y=log2(Abundance),fill=Liver_metastasis)) + 
  facet_grid(Risk ~ Adjuvant_imatinib)+
  scale_fill_manual(values=c("#999999","#E69F00"))+
  labs(title="Brevundimonas")+
  #facet_wrap(~ Adjuvant_imatinib,scales = "free_y")+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5 ,fill = "black")+
  stat_compare_means()+
  theme_bw()

ggsave(plot_Shewanella, filename="boxplot_plot_Brevundimonas_genus.pdf", height=8, width=6)

#####################PCoA four can be modified by Aimtinib genus########################################
matrix_union<-read.csv("union108_genus_all_matrix.csv",header=TRUE)
#matrix_union<-read.csv("taxa_corrected2_Union_matrix.csv",header=TRUE)

#gene_list<-read.csv("selected_genus.csv",header=TRUE)
union_list<-melt(matrix_union)
#union_list<-union_list[union_list$value>0.0001,]
#union_list<-merge(gene_list,union_list,by="X")

colnames(union_list)<-c("genus","samples",'value')
#union_list<-union_list[union_list$genus==c("Stenotrophomonas","Shewanella","Pseudomonas","Brevundimonas"),]

matrix_union<-tapply(union_list$value,list(union_list$genus,union_list$samples),median)
matrix_union[is.na(matrix_union)]<-0

#rownames(matrix_union)<-matrix_union$X
#matrix_union<-matrix_union[,-1]

group_union<-Final_description_patient_info
#rownames(group_union)<-group_union$Specimen_number

otu <- t(matrix_union)
otu.distance <- vegdist(otu)
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)

pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)

p <- ggplot(pc12,aes(x=V1, y=V2))+
  geom_point(size=3)+theme_bw()
p

colnames(group_union)<-c("samples","Gender","Age","Tumor_size","Mitotic_index","Risk", "Adjuvant_imatinib", "Gene_test","PLR","NLR","FIB","PNI","Presence_of.necrosis","Liver_metastasis","Time", "Batch","group")
df <- merge(pc12,group_union,by="samples")

color=c("#56B4E9","#E69F00","#993333","#003399")
dune.div <- adonis2(otu ~ group, data = unique(df[c(1,19)]), permutations = 999, method="bray")

library(pairwiseAdonis)
library(ggrepel)
dune.pairwise.adonis <- pairwise.adonis(x=otu, factors=df$group, 
                                        sim.function = "vegdist",
                                        sim.method = "bray",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)

dune.pairwise.adonis
#dune.pairwise.adonis<-dune.pairwise.adonis[order(dune.pairwise.adonis$p.value),]
dune<-ggplot(dune.pairwise.adonis, aes(x=reorder(pairs,dune.pairwise.adonis$p.value), y=p.value))+
  geom_col()+
  theme(axis.text.x = element_text(angle=60, hjust = 1))+
  coord_flip()+ 
  geom_text(aes(label=p.value), 
            position = position_dodge2(width = 0.9, preserve = 'single'), 
            vjust = -0.2, hjust = 1,color="white")+
  labs(x = "pairs", y = "p.value")+
  theme_bw()

ggsave(dune, file="Pvalue_genus_bar.pdf", width=4, height=3,limitsize = FALSE) 

dune_adonis <- paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)


centroids <- aggregate(cbind(V1,V2)~group,df,mean)
f <- function(z) qt(0.025,df=length(z)-1, lower.tail=F)* sd(z)/sqrt(length(z)) 
se<- aggregate(cbind(se.x=V1,se.y=V2)~group,df,f)
centroids <- merge(centroids,se,by="group")
centroids$AI<-centroids$group
centroids$AI<-ifelse(centroids$AI==c("AI-LM","AI-noLM"),"Adjuvant_imatinib","No_Adjuvant_imatinib")
centroids$LM<-centroids$group
centroids$LM<-ifelse(centroids$group==c("noAI-LM","AI-LM"),"LM","noLM")
centroids

p1<-ggplot(data=df,aes(x=V1,y=V2,color=AI))+
   geom_point(data=centroids,size=10,aes(color=AI))+
   geom_text_repel(data = centroids,aes(label = group),
                   position = "identity",colour="black", size = 3.5)+
   geom_errorbar(data=centroids,aes(ymin=V2-se.y,ymax=V2+se.y),width=0.01)+
   geom_errorbarh(data=centroids,aes(xmin=V1-se.x,xmax=V1+se.x),height=0.01)+
   scale_colour_manual(values = c("#56B4E9","#E69F00"))+
   theme_bw()
  #geom_point(size=1.8)+
  #theme(panel.grid = element_blank())+
  #geom_vline(xintercept = 0,lty="dashed")+
  #geom_hline(yintercept = 0,lty="dashed")+
  #geom_text(aes(label=samples, y=V2+0.02,x=V1+0.02, vjust=0),size=3.5)+
  #guides(color=guide_legend(title="Group"))+
  #labs(x=paste0("PC1 ",pc[1],"%"),
  #     y=paste0("PC2 ",pc[2],"%"),title=dune_adonis)+
  #scale_color_manual(values = color) +
  #scale_fill_manual(values = c("#56B4E9","#E69F00","#993333","#003399"))+
  #theme(axis.title.x=element_text(size=12),
  #      axis.title.y=element_text(size=12,angle=90),
  #      axis.text.y=element_text(size=10),
  #      axis.text.x=element_text(size=10),
  #      panel.grid=element_blank())+
  #      stat_ellipse(type = "norm") +
  #      coord_fixed()
#ggsave(p1, file="Live.metastasis_pcoa.pdf", width=6, height=5,limitsize = FALSE)
ggsave(p1, file="Adjuvant_imatinib_genus_pcoa.pdf", width=6, height=5,limitsize = FALSE)

##################20 feature genus barplot#######################
matrix_union<-read.csv("union108_genus_all_matrix.csv",header=TRUE)
#matrix_union<-read.csv("taxa_corrected2_Union_matrix.csv",header=TRUE)
gene_list<-read.csv("selected_genus9.30.csv",header=TRUE)
feature_genes_matrix<-merge(gene_list,matrix_union,by="X",all=FALSE)
feature_genes_list<-melt(feature_genes_matrix)
colnames(feature_genes_list)<-c("genus","Specimen_number","Relative_abundance")
feature_genes_list<-feature_genes_list[feature_genes_list$Relative_abundance>0.0001,]
feature_genes_with_group<-merge(feature_genes_list,group_union,by="samples",all=FALSE)
feature_genes_with_group$count<-1
feature_genes_with_group$RA<-feature_genes_with_group$count
gender_group<-feature_genes_with_group %>% count(genus, Gender,sort = TRUE)
AI_group<-feature_genes_with_group %>% count(genus, Adjuvant_imatinib,sort = TRUE)
group_group<-feature_genes_with_group %>% count(genus, group,sort = TRUE)
order_genus<-sort(tapply(group_group$n,list(group_group$genus),sum))
order_genus<-as.data.frame(order_genus)

#spearman_correlation<-read.csv("spearman_correlation_genus.csv",header=TRUE)
#spearman_correlation_AI<-spearman_correlation[spearman_correlation$Condition=="LM-AI",]
#spearman_correlation_group<-merge(group_group,spearman_correlation_AI,by="genus")

group_group$genus <- factor(group_group$genus,levels=rownames(order_genus))

bar_feature<-ggplot(data=group_group, aes(x=genus, y=n/108, fill=group))+
  geom_bar(stat = "identity")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  #geom_text(aes(label=P.value), vjust=0,size=1.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

bar_feature

ggsave(bar_feature, file="bar_feature_bar.pdf", width=5, height=6,limitsize = FALSE)


####################lefse#############################
library(microeco)
library(file2meco)
library(ggplot2)

#feature_table <- read.csv('feature_table_risk.csv', row.names = 1)
#sample_table <- read.csv('sample_table.csv', row.names = 1)
#tax_table <- read.csv('tax_table_risk.csv', row.names = 1)

feature_table <- read.csv('otu_tree.csv', row.names = 1)
sample_table <- read.csv('sample_table.csv', row.names = 1)
tax_table <- read.csv('21genera_tree.csv', row.names = 1)
tax_table<-tax_table[,-1]

dataset <- microtable$new(sample_table = sample_table,
                          otu_table = feature_table, 
                          tax_table = tax_table)

lefse <- trans_diff$new(dataset = dataset, 
                        method = 'lefse', 
                        group = "Liver_metastasis", 
                        #p_adjust_method = "fdr",
                        alpha = 1)

lefseplot<-lefse$plot_diff_bar(use_number = 1:50, 
                    width = 0.4, 
                    group_order = c("High","Intermediate","Low")) +
  ggsci::scale_color_jama() +
  ggsci::scale_fill_jama()

ggsave(lefseplot, file="lefseplot_Risk.pdf", width=5, height=6,limitsize = FALSE)

#use_labels <- c("Klebsiella", "Escherichia.Shigella","Empedobacter","Shewanella","Listeria","Sphingobacterium","Pseudomonas","Turicibacter","Stenotrophomonas","Dietzia","Staphylococcus","Chryseobacterium","Dorea","Brevundimonas","Coprococcus","Subdoligranulum","Massilia","Dialister","Subgroup_7","Noviherbaspirillum","TRA3.20")

diff_cladogram<-lefse$plot_diff_cladogram(use_taxa_num = 88, 
                          use_feature_num = 88, 
                          #clade_label_level = 8,
                          select_show_labels = NULL)+
                          ggsci::scale_color_jama() +
                          ggsci::scale_fill_jama()

ggsave(diff_cladogram, file="diff_cladogram_21.pdf", width=10, height=10,limitsize = FALSE)
##### lda score ###########
Risk_matrix<-read.csv("all168samples_matrix_genus.csv",header=TRUE)
Risk_list<-melt(Risk_matrix)
colnames(Risk_list)<-c("genus","samples","RA")
Risk_list<-Risk_list[Risk_list$RA>0.0001,]
Risk_list$RA<-Risk_list$RA*100
Risk_group<-read.csv("sample_table_risk.csv",header=TRUE)
row.names(Risk_group)<-Risk_group$samples
Risk_group<-Risk_group[,-1]
#Risk_a<-merge(Risk_group,Risk_list,by="samples")
Risk_c<-tapply(Risk_list$RA,list(Risk_list$samples,Risk_list$genus),median)
Risk_c[is.na(Risk_c)]<-0
Risk_c<-as.data.frame(Risk_c)

ML_matrix<-read.csv("union108samples_matrix_genus0.01.csv",header=TRUE)
ML_list<-melt(ML_matrix)
colnames(ML_list)<-c("genus","samples","RA")
ML_list<-ML_list[ML_list$RA>0.0001,]
ML_list$RA<-ML_list$RA*100
ML_group<-unique(group_union[c(1,6,14)])
colnames(ML_group)<-c("samples","Risk","Liver_metastasis")
rownames(ML_group)<-ML_group$samples
ML_group<-ML_group[,-1]
ML_group<-as.data.frame(ML_group)
#ML_a<-merge(ML_group,ML_list,by="samples")
ML_c<-tapply(ML_list$RA,list(ML_list$samples,ML_list$genus),median)
ML_c[is.na(ML_c)]<-0
ML_c<-as.data.frame(ML_c)

library(MASS)

Risk.lda <- lda(x =Risk_c, grouping = Risk_group$Risk)
ML.lda <- lda(x =ML_c, grouping = ML_group$Liver_metastasis)

RiskScore<-Risk.lda$means
LMScore<-ML.lda$means


###############heatmap and tree##################
matrix_union<-read.csv("union108_genus_all_matrix.csv",header=TRUE)
#matrix_union<-read.csv("taxa_corrected2_Union_matrix.csv",header=TRUE)
gene_list<-read.csv("ALL_features_true.csv",header=TRUE)
feature_genes_list<-melt(matrix_union)
colnames(feature_genes_list)<-c("Genus","samples","RA")
feature_genes_matrix<-merge(gene_list,feature_genes_list,by="Genus",all=FALSE)


list_group<-merge(feature_genes_matrix,group_union,by="samples",all=FALSE)
list_group<-list_group[list_group$RA>0.0001,]
list_matrix<-tapply(list_group$RA,list(list_group$samples,list_group$Genus),median)
list_matrix[is.na(list_matrix)]<-0
library(psych)
library(stringr)

group<-read.csv("correlation_group.csv",header=TRUE)
rownames(group)<-group$Specimen_number
group<-group[,-1]

otu_genus <-list_matrix
#genus_table<-t(otu_genus)

cortest_psy_sdj <- corr.test(otu_genus, group, method = "spearman", adjust = "fdr")
colnames(cortest_psy_sdj$p.adj)<-c("AI.LM.p.adj","AI.noLM.p.adj","noAI.LM.p.adj","noAI.noLM.p.adj")

##############
OTU_matrix<-read.csv("4group_matrix.csv",header=TRUE)
corr_matrix<-read.csv("cortest_psy_sdj.csv",header=TRUE)
padj_matrix<-read.csv("cortest_psy_sdj.padj.csv",header=TRUE)

OTU_list<-melt(OTU_matrix)
colnames(OTU_list)<-c("group","genus","RA")
corr_list<-melt(corr_matrix)
colnames(corr_list)<-c("genus","group","corr")
padj_list<-melt(padj_matrix)
colnames(padj_list)<-c("genus","group","padj")

plot<-merge(OTU_list,corr_list,by=c("genus","group"))
plot2<-merge(plot,padj_list,by=c("genus","group"))

colnames(gene_list)<-c("genus","feature","Intracellular","Oxygen")

heatmaptoplot<-merge(plot2,gene_list,by="genus")
heatmaptoplot$logRA<-log2(heatmaptoplot$RA*100)
heatmaptoplot$logRA[which(heatmaptoplot$logRA =='-Inf')] <- '0'
write.csv(heatmaptoplot,"heatmaptoplot2.csv")
#################plot heatmap###################
heatmaptoplot<-read.csv("heatmaptoplot.csv",header = TRUE)
group3<-unique(heatmaptoplot[c(1,6,7,8)])
M_p<-tapply(heatmaptoplot$RA,list(heatmaptoplot$group,heatmaptoplot$genus),median)
phr <- hclust(dist(t(M_p))) %>% 
  ggtree(layout="rectangular",branch.length="none")
phc <- hclust(dist(M_p)) %>% 
  ggtree() + layout_dendrogram()

Type_A <- group3$genus %>% as.data.frame() %>%
  mutate(group=group3$feature) %>%
  mutate(group3="")%>%
  ggplot(aes(group3,.,fill=group))+
  scale_fill_brewer(palette="Set3")+
  geom_tile() + 
  scale_x_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank())+
  labs(fill = "feature")

Type_B <- group3$genus %>% as.data.frame() %>%
  mutate(group=group3$Intracellular) %>%
  mutate(group3="")%>%
  ggplot(aes(group3,.,fill=group))+
  scale_fill_brewer(palette="Dark2")+
  geom_tile() + 
  scale_x_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL)+
  labs(fill = "Intracellular")

Type_C <- group3$genus %>% as.data.frame() %>%
  mutate(group=group3$Oxygen) %>%
  mutate(group3="")%>%
  ggplot(aes(group3,.,fill=group))+
  scale_fill_brewer(palette="Set2")+
  geom_tile() + 
  scale_x_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL)+
  labs(fill = "Oxygen")

heatmaptoplot$group <- factor(heatmaptoplot$group,levels=c("AI.LM","noAI.LM","AI.noLM","noAI.noLM"))

p2 <- ggplot(heatmaptoplot,aes(x=group,y=genus))+
  scale_color_gradientn(colours = c('#3A5FCD','white','#FA8072'))+
  theme_bw()+
  geom_point(aes(size=`corr`,
                 color=`logRA`))+
  geom_text(aes(label=padj),col ="black",size=3)+
  theme(panel.grid = element_blank(),axis.text.x =element_text(angle =90,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)+
  geom_vline(xintercept=2.5,size=.8)
  p2<-p2%>%
    insert_left(Type_A,width=.07)%>%
    insert_left(Type_B,width=.07)%>%
    insert_left(Type_C,width=.07)%>%
    insert_left(phr,width=.2)

ggsave(p2, file="heatmap_with_corr.pdf", width=6, height=6,limitsize = FALSE)


###################画热图###########################
matrix_union<-read.csv("union108_genus_all_matrix.csv",header=TRUE)
matrix_all168<-read.csv("all168samples_matrix_genus.csv",header=TRUE)

list_all168<-melt(matrix_all168)
colnames(list_all168)<-c("genus","samples","RA")
list_all168<-list_all168[list_all168$RA>0,]
list_union<-melt(matrix_union)
colnames(list_union)<-c("genus","samples","RA")
list_union<-list_union[list_union$RA>0,]
a1<-merge(list_union,gene_list,by="genus")
a1$count<-1
a2<-merge(list_all168,gene_list,by="genus")
a2$count<-1
b1<-tapply(a1$count,list(a1$genus),sum)
b2<-tapply(a2$count,list(a2$genus),sum)

gene_list<-read.csv("compare_enrichment_socre_2group.csv",header=TRUE)


group_union<-Final_description_patient_info
colnames(group_union)<-c("samples","Gender","Age","Tumor_size","Mitotic_index","Risk", "Adjuvant_imatinib", "Gene_test","PLR","NLR","FIB","PNI","Presence_of.necrosis","Liver_metastasis","Time", "Batch","group")
union_list<-melt(matrix_union)
colnames(union_list)<-c("genus","samples","RA")
union_list<-union_list[union_list$RA>0.0001,]
union_list_9.30<-merge(gene_list,union_list,by="genus",all=FALSE)


plo<-read.csv("compare_effectsize_60genus2group.csv",header=TRUE)

dd<-ggplot(data=plo, aes(x=LM_Size, y=Risk_size))+
  geom_point(aes(color=featureGroup),size=3)+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_vline(xintercept = 0.474,color="grey",lty="dashed")+
  geom_hline(yintercept = 0.474,color="grey",lty="dashed")+
  scale_color_brewer(palette="Set1")+
  #geom_text(aes(label=ifelse(Risk_size>0.474,as.character(genus),'')),hjust=0,vjust=0)+
  #geom_text(aes(label=ifelse(LM_Size>0.474,as.character(genus),'')),hjust=0,vjust=0)+
  geom_smooth(method = 'lm', formula = y ~ x, se = T)+
  stat_cor(data=plo, method = "pearson")+
  theme_bw()+
  theme_minimal()

ggsave(dd, file="dot_with_corr.pdf", width=8, height=7,limitsize = FALSE)

###########################heatmap all Samples############################
############Risk#####
matrix_all168<-read.csv("all168samples_matrix_genus.csv",header=TRUE)
matrix_all168_list<-melt(matrix_all168)
risk_genus<-read.csv("selectGenus_Risk_0.936_.csv",header=TRUE)
group_risk<-read.csv("sample_table_risk.csv",header=TRUE)

risk_genus_select<-merge(risk_genus,matrix_all168_list,by="X")
risk_genus_select<-risk_genus_select[risk_genus_select$value>0.0001,]
colnames(risk_genus_select)<-c("genus","group","Oxygen","samples","RA")
risk_genus_group<-merge(risk_genus_select,group_risk,by="samples")

risk_genus_group$group <- factor(risk_genus_group$group,levels=c("Positive","Negative"))
#risk_genus_group$Risk<-factor(risk_genus_group$Risk,levels=c("High","Intermediate","Low"))

risk_genus_matrix<-tapply(risk_genus_group$RA,list(risk_genus_group$Risk,risk_genus_group$genus),median)
risk_genus_matrix[is.na(risk_genus_matrix)]<-0

g<-unique(risk_genus_group[c(2,3)])
rownames(g)<-g$genus
#g<-g[,-1]


library(pheatmap)
Rh<-pheatmap(risk_genus_matrix,
         scale = "column", 
         
         cluster_row = FALSE,
         clustering_method = "ward.D2")

ggsave(Rh, file="Risk_heatmap.pdf", width=8, height=4,limitsize = FALSE)

#############LM############
matrix_all108<-read.csv("union108_genus_all_matrix.csv",header=TRUE)
matrix_all108_list<-melt(matrix_all108)
LM_genus<-read.csv("108_weight_feature0.930.csv",header=TRUE)
group_LM<-read.csv("Final_patient_info.csv",header=TRUE)
group_LM<-group_LM[c(1,6,14)]
colnames(group_LM)<-c("variable","Risk","Liver_metastasis")
colnames(matrix_all108_list)<-c("genus","variable","value")

lm_genus_select<-merge(LM_genus,matrix_all108_list,by="genus")
lm_genus_select<-lm_genus_select[lm_genus_select$value>0.001,]
#colnames(lm_genus_select)<-c("genus","group","Oxygen","samples","RA")
lm_genus_group<-merge(lm_genus_select,group_LM,by="variable")

lm_genus_group$group <- factor(lm_genus_group$group,levels=c("Positive","Negative"))
#risk_genus_group$Risk<-factor(risk_genus_group$Risk,levels=c("High","Intermediate","Low"))

lm_genus_matrix<-tapply(lm_genus_group$value,list(lm_genus_group$Liver_metastasis,lm_genus_group$genus),mean)
lm_genus_matrix[is.na(lm_genus_matrix)]<-0

Lh<-pheatmap(lm_genus_matrix,
             scale = "column", 
             cluster_row = FALSE,
             #cluster_col = FALSE,
             clustering_method = "ward.D2")

ggsave(Lh, file="LM_heatmap.pdf", width=6, height=3,limitsize = FALSE)


c<-aggregate(group_union$count,list(group_union$Batch,group_union$Adjuvant_imatinib,group_union$Liver_metastasis),sum)
colnames(c)<-c("Batch","Adjuvant_imatinib","Liver_metastasis","count")

##################feature genera tree##############
library(ggtree)
library(treeio)
library(ggsci)
library(ggraph)
library(igraph)
library(tidyverse)
library(ggcor)
library(readr)
library(dplyr)
library(igraph)
library(ggraph)
library(colorspace)

otu <- read.csv('otu_tree.csv',header=TRUE)
#row.names(otu)<-otu$Genus
#otu<-otu[,-1]
otu_melt<-melt(otu)
colnames(otu_melt)<-c("Genus","Specimen_number","RA")
otu_melt<-otu_melt[otu_melt$RA>0.01,]
otu_melt$RA<-log10(otu_melt$RA)
otu_melt<-merge(otu_melt,Final_description_patient_info,by="Specimen_number",all=FALSE)
otu_melt$Genus<- factor(otu_melt$Genus,levels=c("g_Brevundimonas","g_Escherichia.Shigella","g_Klebsiella","g_Massilia","g_Noviherbaspirillum","g_TRA3.20","g_Pseudomonas","g_Shewanella","g_Stenotrophomonas","g_Chryseobacterium","g_Empedobacter","g_Sphingobacterium","g_Coprococcus","g_Dorea","g_Subdoligranulum","g_Dialister","g_Listeria","g_Staphylococcus","g_Turicibacter","g_Dietzia","g_Subgroup_7"))
otu<-tapply(otu_melt$RA,list(otu_melt$Genus,otu_melt$Batch),median)
otu[is.na(otu)]<-0

da <- cor_tbl(otu,cluster = F)
d2 <- data.frame(otu_melt$Liver_metastasis)

b<-quickcor(da,circular = T,cluster = F,grid.colour = 'white',inner = 3,
            outer = 0.3,open = 0.1,width = 0.2) +
  geom_colour(colour = FALSE) +
  scale_fill_continuous_sequential(palette = "Viridis", rev = FALSE) +
  guides(fill = guide_colorbar(title = 'log10 RA'))+
  #anno_col_tree() +
  #anno_row_tree() +
  set_p_yaxis()
  #set_p_xaxis()

 b              
ggsave(b, file="heatmap_D2.pdf", width=12, height=12,limitsize = FALSE)  

group_tree<-read.csv("21genera_tree.csv",header=TRUE)
#group_tree<-group_tree[,-1]


edges_level1_2 <- group_tree %>% select(
  Kingdom, Phylum) %>% unique %>% rename(
    from=Kingdom, to=Phylum)


edges_level2_3 <- group_tree %>% select(
  Phylum, Class) %>% unique %>% rename(
    from=Phylum, to=Class)


edges_level3_4 <- group_tree %>% select(
  Class, Order) %>% unique %>% rename(
    from=Class, to=Order)


edges_level4_5 <- group_tree %>% select(
  Order,Family) %>% unique %>% rename(
    from=Order, to=Family)


edges_level5_6 <- group_tree %>% select(
  Family,Genus) %>% unique %>% rename(
    from=Family, to=Genus)



edge_list=rbind(edges_level1_2, edges_level2_3,edges_level3_4,edges_level4_5,edges_level5_6)
edge_list<-merge(edge_list,lll,by="to",all=FALSE)

edge_list<-cbind(edge_list[2],edge_list[1],edge_list[3])

class_l<-read.csv("21genera_class.csv",header=TRUE)

edge_list<-merge(edge_list,class_l,by="to",all=FALSE)

name <- unique(c(as.character(edge_list$from), as.character(edge_list$to)))
#graph = graph_from_data_frame(edge_list$to, vertices = edge_list$from)
set_graph_style()
#edge_list <- tbl_graph(edge_list$from, edge_list$class)
lay<-create_layout(edge_list, 'circlepack')
#edge_list$RA<-as.numeric(edge_list$RA)
edge_list<-unique(edge_list)
colnames(lll)<-c("name","RA")

edges <- edge_list %>%
  mutate(corr = sample(-1:1, size = n(), replace = TRUE))

nodes <- data.frame(
  name = unique(union(edges$from, edges$to))
)
nodes<-merge(nodes,lll,by="name")
library( tidygraph)
nodes
g <- tbl_graph(nodes = nodes, edges = edges)
flareGraph <- tbl_graph(flare$vertices, flare$edges) %>%
  mutate(
    class = map_bfs_chr(node_is_root(), .f = function(node, dist, path, ...) {
      if (dist <= 1) {
        return(shortName[node])
      }
      path$result[[nrow(path)]]
    })
  )
#extrafont::loadfonts()
a<-ggraph(g, 'partition', circular = TRUE,height=0.5) +
  geom_node_arc_bar(aes(fill = log10(RA)),linewidth =0.1) +
  #geom_node_tile()+
  coord_fixed()+
  theme_void()+
  #geom_node_text(aes(label=name),size=2.5,position = "identity")+
  geom_node_text(aes(label=name), size=3.5) +
  scale_fill_continuous_sequential(palette = "PinkYl", rev = TRUE)
  #theme_graph()
  #theme_void()
a<-a+ggraph(g, 'dendrogram', circular = TRUE) + 
  geom_edge_elbow() +
  layout_tbl_graph_unrooted()+
  coord_fixed() +
  theme_graph()
a

ggsave(a, file="21genera_tree_ref_shape_text.pdf", width=16, height=8,limitsize = FALSE)   

g#########################treemap###############
RA_matrxi<-read.csv("otu_tree.csv")
RA_l<-melt(RA_matrxi)
RA_matrxi<-aggregate(RA_l$value,list(RA_l$Genus),median)
colnames(RA_matrxi)<-c("Genus","RA")
group_list<-read.csv("21genera_tree.csv")
all_sun<-merge(RA_matrxi,group_list,by="Genus")
all_sum_list<-melt(all_sun)
genus_sum<-aggregate(all_sum_list$value,list(all_sum_list$Genus),sum)
kingdom_sum<-aggregate(all_sum_list$value,list(all_sum_list$Kingdom),sum)
phylum_sum<-aggregate(all_sum_list$value,list(all_sum_list$Phylum),sum)
class_sum<-aggregate(all_sum_list$value,list(all_sum_list$Class),sum)
family_sum<-aggregate(all_sum_list$value,list(all_sum_list$Family),sum)
order_sum<-aggregate(all_sum_list$value,list(all_sum_list$Order),sum)
lll<-rbind(genus_sum,family_sum,class_sum,order_sum,phylum_sum,kingdom_sum)
colnames(lll)<-c("to","RA")

