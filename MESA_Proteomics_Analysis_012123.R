rm(list=ls())

library(tidyverse)
library(ggforce)
library(readxl)
library(haven)
library(gtsummary)
library(OlinkAnalyze)
library(bartMachine)
library(bartMachineJARs)
library(ggcorrplot)
library(HDtest)
library(WGCNA)
library(RColorBrewer)
library(scmap)
library(dendextend)


setwd("/file/to/path")

##############################
# 0. settings for ggplot#
##############################
hw <- theme_gray()+ theme(
  plot.title=element_text(hjust=0.5,size=18,face="bold"),
  plot.subtitle=element_text(hjust=0.5,size=12),
  plot.caption=element_text(hjust=-.5,size=10),
  strip.background=element_rect(fill=rgb(.9,.95,1),
                                colour=gray(.5), linewidth=.2),
  panel.border=element_rect(fill=FALSE,colour=gray(.70)),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.spacing.x = unit(0.2,"cm"),
  panel.spacing.y = unit(0.2,"cm"),
  axis.text=element_text(colour="black",size=10),
  axis.text.y=element_text(margin=ggplot2::margin(0,3,0,3)),
  axis.text.x=element_text(margin=ggplot2::margin(-1,0,3,0)),
  axis.title=element_text(size=16,face="bold"),
  legend.text=element_text(size=14),
  legend.title = element_blank(),
  legend.position='none'
)

##############################
# 1. read in data#
##############################

#lge status
lge=na.omit(read_excel("LGE_matched_new_2-2020.xlsx"));
names(lge)=c('LGEneg','LGEpos');
lge_long=pivot_longer(lge,cols=1:2,names_to='LGEstatus',values_to='ID');
lge_long=lge_long[,2:1];
lge_long$LGEstatus = factor(lge_long$LGEstatus,levels=c('LGEneg','LGEpos'));

#ecv status
ecv=read_excel("High_ECV_new_2-2020.xlsx");
names(ecv)=c('ID','ECVtopq');
ecv$ECVtopq = factor(ecv$ECVtopq,levels=c(0,1));

#olink
olink_dta=read_excel("20191239_deFilippi_NPX.xlsx", range = "A8:CO503", col_names = FALSE);
olink_colnames=read_excel("20191239_deFilippi_NPX.xlsx", range = "B4:CO4", col_names = FALSE);
olink_olinkid=read_excel("20191239_deFilippi_NPX.xlsx", range = "B4:CO6", col_names = FALSE);
olink_olinkid = as.data.frame(t(olink_olinkid));
colnames(olink_olinkid) = c('Assay','UniProt','OlinkID')

olink_lod=as.numeric(read_excel("20191239_deFilippi_NPX.xlsx", range = "B505:CO505", col_names = FALSE));
names(olink_dta)=c('ID',as.character(olink_colnames));

#mesa
mesa_dta=read_dta('Liver and cardiac T1 data.dta');
mesa_dta=mesa_dta[,c('idno','age1c','gender1','race1c','bmi5c','sbp5c','dbp5c','cig5c',
                     'chol5','hdl5','dm035c','htnmed5c','lipid5c','income5',
                     'cepgfr5c','ecvsyn','myoscar5')]

#recodings
#combine treated and untreated diabetes
mesa_dta$gcmetstatus=factor(recode(as.numeric(mesa_dta$dm035c),'3'=2),levels=c(0,1,2)); 
mesa_dta$incomerc=factor(recode(as.numeric(mesa_dta$income5),
                                '1'=0,'2'=0,'3'=0,'4'=0,'5'=0,
                                '6'=1,'7'=1,'8'=1,'9'=1,'10'=1,
                                '11'=2,'12'=2,'13'=2,'14'=2,'15'=2),
                         levels=c(0,1,2)); #collapse to <$20k,[$20k,50k),$50k+
names(mesa_dta)[1]='ID';

#class variables approriately
mesa_dta$ID = as.integer(mesa_dta$ID);
mesa_dta$age1c = as.integer(mesa_dta$age1c);
mesa_dta$gender1 = factor(mesa_dta$gender1,levels=c(0,1));
mesa_dta$race1c = factor(mesa_dta$race1c,levels=c(1,2,3,4));
mesa_dta$bmi5c = as.numeric(mesa_dta$bmi5c);
mesa_dta$sbp5c = as.numeric(mesa_dta$sbp5c);
mesa_dta$dbp5c = as.numeric(mesa_dta$dbp5c);
mesa_dta$cig5c = factor(mesa_dta$cig5c,levels=c(0,1,2));
mesa_dta$chol5 = as.numeric(mesa_dta$chol5);
mesa_dta$hdl5 = as.numeric(mesa_dta$hdl5);
mesa_dta$dm035c = factor(mesa_dta$dm035c,levels=c(0,1,2,3));
mesa_dta$htnmed5c = factor(mesa_dta$htnmed5c,levels=c(0,1));
mesa_dta$lipid5c = factor(mesa_dta$lipid5c,levels=c(0,1));
mesa_dta$income5 = factor(mesa_dta$income5,levels=1:15);
mesa_dta$cepgfr5c = as.numeric(mesa_dta$cepgfr5c);
mesa_dta$ecvsyn = as.numeric(mesa_dta$ecvsyn);
mesa_dta$myoscar5 = factor(na_if(mesa_dta$myoscar5,9),
                           levels=c(0,1));

#merge mesa and ecv/lge
ecvmesa_dta=merge(x=ecv,y=mesa_dta,by='ID',all.x=TRUE);
lgemesa_dta=merge(x=lge_long,y=mesa_dta,by='ID',all.x=TRUE);

#merge olink and ecv/lge
ecv_dtawide=merge(x=ecv,y=olink_dta,by='ID',all.x=TRUE);
lge_dtawide=merge(x=lge_long,y=olink_dta,by='ID',all.x=TRUE);

#convert to long format
lge_dtalong=pivot_longer(data=lge_dtawide,cols='TNFRSF14':'CCL16',
                         names_to='Protein',values_to='ExpressionLevel')
ecv_dtalong=pivot_longer(data=ecv_dtawide,cols='TNFRSF14':'CCL16',
                         names_to='Protein',values_to='ExpressionLevel')

##############################
#2. violin plots#
##############################
p1=ggplot(data=ecv_dtalong,mapping=aes(x=ECVtopq, y=ExpressionLevel, fill=ECVtopq))+
  geom_violin()+facet_wrap_paginate('Protein',scales='free',ncol=3,nrow=4,page=1)+hw;

pdf("Output/SuppFigure1_ECVviolinplots.pdf", onefile = TRUE)
for (i in 1:n_pages(p1)){
  print(ggplot(data=ecv_dtalong,mapping=aes(x=ECVtopq, y=ExpressionLevel, fill=ECVtopq))+
          geom_violin()+facet_wrap_paginate('Protein',scales='free',ncol=3,nrow=4,page=i)+hw);
}
dev.off()

p2=ggplot(data=lge_dtalong,mapping=aes(x=LGEstatus, y=ExpressionLevel, fill=LGEstatus))+
  geom_violin()+facet_wrap_paginate('Protein',scales='free',ncol=3,nrow=4,page=1)+hw;

pdf("Output/SuppFigure2_LGEviolinplots.pdf", onefile = TRUE)
for (i in 1:n_pages(p2)){
  print(ggplot(data=lge_dtalong,mapping=aes(x=LGEstatus, y=ExpressionLevel, fill=LGEstatus))+
          geom_violin()+facet_wrap_paginate('Protein',scales='free',ncol=3,nrow=4,page=i)+hw);
}
dev.off()

##############################
#3. summary tables#
##############################
t1 <- tbl_summary(ecvmesa_dta[,c(-1,-15)],
            by=ECVtopq,
            type = all_continuous() ~ "continuous2",
            statistic = list(all_continuous() ~ c("{mean} ({sd})",
                                                  "{median} ({p25}, {p75})", 
                                                  "{min}, {max}"),
                             all_categorical() ~ "{n} / {N} ({p}%)")) %>% 
  add_p(list(all_continuous() ~ "t.test",all_categorical() ~ "fisher.test")) %>% add_n()
t1
gt::gtsave(as_gt(t1), file = "Output/Table1_ECV.html")
gt::gtsave(as_gt(t1), file = "Output/Table1_ECV.pdf")

t2 <- tbl_summary(lgemesa_dta[,c(-1,-15)],by=LGEstatus,
            type = all_continuous() ~ "continuous2",
            statistic = list(all_continuous() ~ c("{mean} ({sd})",
                                                  "{median} ({p25}, {p75})", 
                                                  "{min}, {max}"),
                             all_categorical() ~ "{n} / {N} ({p}%)")) %>% 
  add_p(list(all_continuous() ~ "t.test",all_categorical() ~ "fisher.test")) %>% add_n()
t2
gt::gtsave(as_gt(t2), file = "Output/Table2_LGE.html")
gt::gtsave(as_gt(t2), file = "Output/Table2_LGE.pdf")

t3 <- tbl_summary(ecv_dtawide[,-1],by=ECVtopq,
            type = all_continuous() ~ "continuous2",
            statistic = list(all_continuous() ~ c("{mean} ({sd})",
                                                  "{median} ({p25}, {p75})", 
                                                  "{min}, {max}"),
                             all_categorical() ~ "{n} / {N} ({p}%)")) %>% 
  add_p(all_continuous() ~ "t.test") %>% add_n()
t3
gt::gtsave(as_gt(t3), file = "Output/Table3_ECV_proteins.html")
gt::gtsave(as_gt(t3), file = "Output/Table3_ECV_proteins.pdf")

t4 <- tbl_summary(lge_dtawide[,-1],by=LGEstatus,
            type = all_continuous() ~ "continuous2",
            statistic = list(all_continuous() ~ c("{mean} ({sd})",
                                                  "{median} ({p25}, {p75})", 
                                                  "{min}, {max}"),
                             all_categorical() ~ "{n} / {N} ({p}%)")) %>% 
  add_p(all_continuous() ~ "t.test") %>% add_n()
t4
gt::gtsave(as_gt(t4), file = "Output/Table4_LGE_proteins.html")
gt::gtsave(as_gt(t4), file = "Output/Table4_LGE_proteins.pdf")

##############################
#4. volcano plots#
##############################
#ECV#
#format data for use of olink ttest package
ecv_dtaolinkform=merge(ecv_dtalong,olink_olinkid,by.x='Protein',by.y='Assay')
colnames(ecv_dtaolinkform)[1:2]=c('Assay','SampleID');
colnames(ecv_dtaolinkform)[4] = 'NPX'
ecv_dtaolinkform$Index = rep(seq(1:unique(table(ecv_dtaolinkform$Assay))),
                             length(unique(ecv_dtaolinkform$Assay)))
ecv_dtaolinkform$Panel = 'Olink CARDIOVASCULAR III(v.6113)'
ecv_dtaolinkform$ECVtopq = factor(ecv_dtaolinkform$ECVtopq,levels=c(1,0));

#use olink ttest package to produce ttest results
ttest_results=olink_ttest(df=ecv_dtaolinkform,
                          variable='ECVtopq',
                          alternative = 'two.sided')
ttest_results$estimate = log2(ttest_results$'1'/ttest_results$'0')
ttest_results$Significance = factor("Not Significant",levels=c("Not Significant",
                                                               "Significant (Unadjusted)",
                                                               "Significant (Adjusted)"))
ttest_results$Significance[ttest_results$Adjusted_pval<=0.05] = "Significant (Adjusted)"
ttest_results$Significance[ttest_results$Adjusted_pval>0.05 & 
                             ttest_results$p.value<0.05] = "Significant (Unadjusted)"

olinkid_list = ttest_results$OlinkID[ttest_results$p.value<0.05]

#plot volcano plot
ecv_volcano_plot <- ttest_results %>%
  ggplot2::ggplot(ggplot2::aes(x = estimate, y = -log10(p.value),
                               color = Significance)) +
  ggplot2::geom_point() +
  ggplot2::labs(x = "Log2-fold change", y = "-log10(p-value)") +
  ggrepel::geom_label_repel(data = subset(ttest_results, OlinkID %in% olinkid_list),
                            ggplot2::aes(label = Assay), box.padding = 1, show.legend = FALSE,max.overlaps=Inf) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype="dotted") +
  OlinkAnalyze::set_plot_theme() + ggtitle('ECV Volcano Plot')+ylim(c(0,8.5))+xlim(c(-0.3,0.3))+
  scale_color_manual(values=c("grey70","darkblue","red"))
pdf("Output/Figure1a_ECVVolcanoPlot.pdf", onefile = TRUE)
print(ecv_volcano_plot)
dev.off()


#LGE#
#format data for use of olink ttest package
lge_dtaolinkform=merge(lge_dtalong,olink_olinkid,by.x='Protein',by.y='Assay')
colnames(lge_dtaolinkform)[1:2]=c('Assay','SampleID');
colnames(lge_dtaolinkform)[4] = 'NPX'
lge_dtaolinkform$Index = rep(seq(1:unique(table(lge_dtaolinkform$Assay))),
                             length(unique(lge_dtaolinkform$Assay)))
lge_dtaolinkform$Panel = 'Olink CARDIOVASCULAR III(v.6113)'
lge_dtaolinkform$LGEstatus = factor(lge_dtaolinkform$LGEstatus,levels=c('LGEpos','LGEneg'));

#use olink ttest package to produce ttest results
ttest_results=olink_ttest(df=lge_dtaolinkform,
                          variable='LGEstatus',
                          alternative = 'two.sided')
ttest_results$estimate = log2(ttest_results$LGEpos/ttest_results$LGEneg)
ttest_results$Significance = factor("Not Significant",levels=c("Not Significant",
                                                               "Significant (Unadjusted)",
                                                               "Significant (Adjusted)"))
ttest_results$Significance[ttest_results$Adjusted_pval<=0.05] = "Significant (Adjusted)"
ttest_results$Significance[ttest_results$Adjusted_pval>0.05 & 
                             ttest_results$p.value<0.05] = "Significant (Unadjusted)"

olinkid_list = ttest_results$OlinkID[ttest_results$p.value<0.05]

#plot volcano plot
lge_volcano_plot <- ttest_results %>%
  ggplot2::ggplot(ggplot2::aes(x = estimate, y = -log10(p.value),
                               color = Significance)) +
  ggplot2::geom_point() +
  ggplot2::labs(x = "Log2-fold change", y = "-log10(p-value)") +
  ggrepel::geom_label_repel(data = subset(ttest_results, OlinkID %in% olinkid_list),
                            ggplot2::aes(label = Assay), box.padding = 1, show.legend = FALSE,max.overlaps=Inf) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype="dotted") +
  OlinkAnalyze::set_plot_theme() + ggtitle('LGE Volcano Plot')+ylim(c(0,8.5))+xlim(c(-0.3,0.3))+
  scale_color_manual(values=c("grey70","darkblue","red"))
pdf("Output/Figure1b_LGEVolcanoPlot.pdf", onefile = TRUE)
print(lge_volcano_plot)
dev.off()

##############################
#4. BART model and variable importance plots#
##############################
##ECV

seed = 46

set.seed(seed)
options(java.parameters = "-Xmx4000m")
set_bart_machine_num_cores(1)


#model fit
y <- factor(ecv_dtawide$ECVtopq[-c(115,285)],levels=c(1,0))
X <- cbind(ecv_dtawide[-c(115,285),3:ncol(ecv_dtawide)],
           ecvmesa_dta[-c(115,285),c(3:14,16,18:20)])
bart_machine_ecv <- bartMachine(X, y,
                            num_trees=100,
                            num_burn_in=500,
                            num_iterations_after_burn_in=2000,
                            use_missing_data=TRUE,
                            mem_cache_for_speed = FALSE,
                            seed=seed)
bart_machine_ecv
pdf("Output/BART_ECV_ConvergenceDiagnostics.pdf", onefile = TRUE)
plot_convergence_diagnostics(bart_machine_ecv)
dev.off()

#variable importance plot (proteins only)
bart_machine_ecv_vi <- investigate_var_importance(bart_machine_ecv,
                                              num_replicates_for_avg = 20,
                                              num_trees_bottleneck = 20,num_var_plot = 30)

protidx = (names(bart_machine_ecv_vi$avg_var_props) %in% colnames(olink_dta))
protidx[(which(cumsum(protidx)==20)[1]+1):length(protidx)] = FALSE
moe=1.96 * bart_machine_ecv_vi$sd_var_props[protidx] / sqrt(20)

pdf("Output/Figure2a_ECV_VariableImportance.pdf", onefile = TRUE)
par(mar = c(7, 6, 3, 0))
bars = barplot(bart_machine_ecv_vi$avg_var_props[protidx], 
               names.arg = names(bart_machine_ecv_vi$avg_var_props)[protidx], 
               las = 2, 
               ylab = "Inclusion Proportion", 
               col = "gray",
               ylim = c(0, max(bart_machine_ecv_vi$avg_var_props[protidx] + moe))
)
conf_upper = bart_machine_ecv_vi$avg_var_props[protidx] + 1.96 * bart_machine_ecv_vi$sd_var_props[protidx] / sqrt(20)
conf_lower = bart_machine_ecv_vi$avg_var_props[protidx] - 1.96 * bart_machine_ecv_vi$sd_var_props[protidx] / sqrt(20)
segments(bars, bart_machine_ecv_vi$avg_var_props[protidx], bars, conf_upper, col = rgb(0.59, 0.39, 0.39), lwd = 3) # Draw error bars
segments(bars, bart_machine_ecv_vi$avg_var_props[protidx], bars, conf_lower, col = rgb(0.59, 0.39, 0.39), lwd = 3)
dev.off()

#significance testing for top 3 vars
# WARNING: THIS CODE CAN TAKE A FEW HOURS TO RUN
pdf("Output/BART_ECV_Top3_Significance.pdf", onefile = TRUE)
par(mar = c(5.1, 4.1, 4.1, 2.1))
topvars <- c("PAI","NT-proBNP","IGFBP-1")
testresults <- cov_importance_test(bart_machine_ecv, covariates = topvars,
                                   num_permutation_samples = 1000, plot=TRUE)
dev.off()

##LGE

seed = 389

set.seed(seed)
options(java.parameters = "-Xmx4000m")
set_bart_machine_num_cores(1)

y2 <- factor(lge_dtawide$LGEstatus,levels=c('LGEpos','LGEneg'))
X2 <- cbind(lge_dtawide[,3:ncol(lge_dtawide)],
            lgemesa_dta[,c(3:14,16,19:20)])

#model fit
bart_machine_lge <- bartMachine(X2, y2,
                            num_trees=100,
                            num_burn_in=500,
                            num_iterations_after_burn_in=2000,
                            use_missing_data=TRUE,
                            mem_cache_for_speed = FALSE,
                            seed=seed)
bart_machine_lge

pdf("Output/BART_LGE_ConvergenceDiagnostics.pdf", onefile = TRUE)
plot_convergence_diagnostics(bart_machine_lge)
dev.off()

#variable importance plot (proteins only)
bart_machine_lge_vi <- investigate_var_importance(bart_machine_lge,
                                                  num_replicates_for_avg = 20,
                                                  num_trees_bottleneck = 20,num_var_plot = 30)

protidx = (names(bart_machine_lge_vi$avg_var_props) %in% colnames(olink_dta))
protidx[(which(cumsum(protidx)==20)[1]+1):length(protidx)] = FALSE
moe=1.96 * bart_machine_lge_vi$sd_var_props[protidx] / sqrt(20)

pdf("Output/Figure2b_LGE_VariableImportance.pdf", onefile = TRUE)
par(mar = c(7, 6, 3, 0))
bars = barplot(bart_machine_lge_vi$avg_var_props[protidx], 
               names.arg = names(bart_machine_lge_vi$avg_var_props)[protidx], 
               las = 2, 
               ylab = "Inclusion Proportion", 
               col = "gray",
               ylim = c(0, max(bart_machine_lge_vi$avg_var_props[protidx] + moe))
)
conf_upper = bart_machine_lge_vi$avg_var_props[protidx] + 1.96 * bart_machine_lge_vi$sd_var_props[protidx] / sqrt(20)
conf_lower = bart_machine_lge_vi$avg_var_props[protidx] - 1.96 * bart_machine_lge_vi$sd_var_props[protidx] / sqrt(20)
segments(bars, bart_machine_lge_vi$avg_var_props[protidx], bars, conf_upper, col = rgb(0.59, 0.39, 0.39), lwd = 3) # Draw error bars
segments(bars, bart_machine_lge_vi$avg_var_props[protidx], bars, conf_lower, col = rgb(0.59, 0.39, 0.39), lwd = 3)
dev.off()

#significance testing for top 2 vars
# WARNING: THIS CODE CAN TAKE A FEW HOURS TO RUN
pdf("Output/BART_LGE_Top3_Significance.pdf", onefile = TRUE)
par(mar = c(5.1, 4.1, 4.1, 2.1))
topvars <- c("NT-proBNP","IGFBP-1")
testresults <- cov_importance_test(bart_machine_lge, covariates = topvars,
                                   num_permutation_samples = 1000, plot=TRUE)
dev.off()

##############################
#5. Covariance matrix significance testing#
##############################

#ECV
#visualize correlation matrices
p1 <- ggcorrplot(cor(ecv_dtawide[ecv_dtawide$ECVtopq==0,c(-1,-2,-36,-41)]),
              title="ECVq1",hc.order=TRUE,hc.method="complete")+hw+
  theme(axis.text = element_text(size=5), axis.title=element_blank(),
        axis.text.x=element_text(angle=90))

pdf("Output/ECVq1_correlationplot.pdf", onefile = TRUE)
p1
dev.off()

#get order of proteins
cormat2 = cor(ecv_dtawide[ecv_dtawide$ECVtopq==1,as.character(p1$data[1:90,1])])

p2=ggcorrplot(cormat2,
              title="ECVq4",hc.order=FALSE,hc.method="complete")+hw+
  theme(axis.text = element_text(size=5), axis.title=element_blank(),
        axis.text.x=element_text(angle=90))
pdf("Output/ECVq4_correlationplot.pdf", onefile = TRUE)
p2
dev.off()

#test equality of covariance matrices ECV
testCov(ecv_dtawide[ecv_dtawide$ECVtopq==0,as.character(p1$data[1:90,1])],
        ecv_dtawide[ecv_dtawide$ECVtopq==1,as.character(p1$data[1:90,1])],
        method="HD",J=5000,alpha=0.05,n.core=1)

#LGE
#visualize correlation matrices
p1=ggcorrplot(cor(lge_dtawide[lge_dtawide$LGEstatus=='LGEneg',c(-1,-2)]),
              title="LGE negative",hc.order=TRUE,hc.method="complete")+hw+
  theme(axis.text = element_text(size=5), axis.title=element_blank(),
        axis.text.x=element_text(angle=90))

pdf("Output/LGEneg_correlationplot.pdf", onefile = TRUE)
p1
dev.off()

#get order of proteins
cormat = cor(lge_dtawide[lge_dtawide$LGEstatus=='LGEpos',as.character(p1$data[1:92,1])])

p2=ggcorrplot(cormat,
              title="LGE positive",hc.order=FALSE,hc.method="complete")+hw+
  theme(axis.text = element_text(size=5), axis.title=element_blank(),
        axis.text.x=element_text(angle=90))

pdf("Output/LGEpos_correlationplot.pdf", onefile = TRUE)
p2
dev.off()

#test equality of covariance matrices LGE
testCov(lge_dtawide[lge_dtawide$LGEstatus=='LGEneg',as.character(p1$data[1:92,1])],
        lge_dtawide[lge_dtawide$LGEstatus=='LGEpos',as.character(p1$data[1:92,1])],
        method="HD",J=5000,alpha=0.05,n.core=4)

##############################
#6. Hierarchical clustering dendrograms, Sankey plot, and tanglegram#
##############################

#ECVq1
adj1.ecvq1=cor(ecv_dtawide[ecv_dtawide$ECVtopq==0,c(-1,-2)],
                use="pairwise.complete.obs")

# Turn adjacency into a measure of dissimilarity and cluster
dissADJ=(1-adj1.ecvq1)/2
hierADJ=hclust(as.dist(dissADJ), method="complete" )
ecvq1.dend = as.dendrogram(hierADJ)

# Plot the dendrogram with module colors
height=sort(hierADJ$height,decreasing=TRUE)[5]
clusthgt=cutreeStatic(hierADJ, cutHeight=height, minSize=1)
ECVq1=brewer.pal(5,"Set3")[clusthgt]

pdf("Output/Figure3a_ECVq1_dendrogram.pdf", onefile = TRUE)
plotDendroAndColors(hierADJ, colors = data.frame(ECVq1),
                    dendroLabels = NULL, abHeight = height,
                    main = "Protein dendrogram and clusters (ECVq1)")
dev.off()

#get clusters
ecvq1clust=data.frame(protein=hierADJ$labels,
                       ecvq1cluster=cutreeStatic(hierADJ, cutHeight=height, minSize=1))

#ECVq4
adj1.ecvq4=cor(ecv_dtawide[ecv_dtawide$ECVtopq==1,c(-1,-2)],
                use="pairwise.complete.obs")

# Turn adjacency into a measure of dissimilarity and cluster
dissADJ=(1-adj1.ecvq4)/2
hierADJ=hclust(as.dist(dissADJ), method="complete" )
ecvq4.dend = as.dendrogram(hierADJ)

# Plot the dendrogram with module colors
height=sort(hierADJ$height,decreasing=TRUE)[5]
clusthgt=cutreeStatic(hierADJ, cutHeight=height, minSize=1)
ECVq4=brewer.pal(5,"Set3")[clusthgt]

pdf("Output/Figure3b_ECVq4_dendrogram.pdf", onefile = TRUE)
plotDendroAndColors(hierADJ, colors = data.frame(ECVq4),
                    dendroLabels = NULL, abHeight = height,
                    main = "Protein dendrogram and clusters (ECVq4)")
dev.off()

#get clusters
ecvq4clust=data.frame(protein=hierADJ$labels,
                      ecvq4cluster=cutreeStatic(hierADJ, cutHeight=height, minSize=1))

#Sankey plot for ECVq1 and ECVq4 clusterings
clustcomp=cbind(ecvq1clust,ecvq4clust)
clustcomp=clustcomp[,-3]
plot(getSankey(clustcomp$ecvq1cluster,clustcomp$ecvq4cluster,
               colors=brewer.pal(5,"Set3")))

#dendrogram comparison for ECVq1 and ECVq4
dl.ecv <- dendlist(ecvq1.dend, ecvq4.dend)

#measure entanglement and plot tanglegram for best alignment
entanglement(dl.ecv)
dl.ecv.ut=untangle(dl.ecv,"step2side") #step2side gives the lowest entanglement score so best alignment
entanglement(dl.ecv.ut)

pdf("Output/SuppFigure3_ECV_tanglegram.pdf", onefile = TRUE)
tanglegram(dl.ecv.ut,sort=FALSE,
           color_lines=brewer.pal(9,"Set1"),
           main_left="ECV Q1",
           main_right="ECV Q4",
           highlight_branches_lwd=FALSE,
           highlight_branches_col=FALSE,
           common_subtrees_color_branches=TRUE,
           lab.cex=.5)
dev.off()




