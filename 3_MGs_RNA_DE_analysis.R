#########################################################################################
# Molecular glues - DE analysis applied on loess-normalized gene expression data
#########################################################################################

rm(list=ls())

require(venn)
require(vioplot)

library(extrafont)
loadfonts(device = "pdf")
fonts()

library(heatmap3)
library("affy")
library("edgeR")

library("limma")

require(ggplot2)



########################################################################################
# load data table with gene names (without ERCCs/SIRVs) after nomralization:
########################################################################################


# load matrix including ERCC-loess-normalized gene expression values for 18 samples (only HGNC genes)

count_df=as.data.frame(read.table("MGs_RNAseq_CPM_ERCC_loess_normalized_matrix.HGNC.only_high_conc.txt",sep="\t", header=TRUE, row.names=1, dec=".", fill = TRUE))



#Used Colors

#DMSO: 153, 153, 153, #999999
#THZ: 222, 241, 245, #DEF1F5
#Hit 4: 241, 146, 34, #F19222
#Hit1: 4, 104, 145, #046891
#Hit2: 94, 195, 235, #5EC3EB
#Hit8: 0, 48, 72, #003047

colors_MGs_18=c(rep("#999999",3),rep("#DEF1F5",3),rep("#F19222",3),rep("#046891",3),rep("#5EC3EB",3),rep("#003047",3))



########################################################################################
# Boxplots representing loess-normalized values of gene expression per sample
########################################################################################

# Boxplots without outliers

pdf("boxplots_reads_ERCCnorm.without_outliers.pdf",6,4)
par(cex.axis=1,cex.main=1,family="Arial",mar=c(7,4.5,4,2))

boxplot(count_df,outline=F
        ,main="Genes, Loess-normalization (log2 cpm)",las=2
        ,names=c("DMSO_rep1","DMSO_rep2","DMSO_rep3","THZ531_rep1","THZ531_rep2","THZ531_rep3","Hit4_rep1","Hit4_rep2","Hit4_rep3","Hit1_rep1","Hit1_rep2","Hit1_rep3","Hit2_rep1","Hit2_rep2","Hit2_rep3","Hit8_rep1","Hit8_rep2","Hit8_rep3")
        ,col=colors_MGs_18,ylab="Loess-normalized log2 cpm")
dev.off()


# Boxplots with outliers

pdf("boxplots_reads_ERCCnorm.with_outliers.pdf",6,4)
par(cex.axis=1,cex.main=1,family="Arial",mar=c(7,4.5,4,2))

boxplot(count_df,outline=T
        ,main="Genes, loess-normalization (log2 cpm)",las=2
        ,names=c("DMSO_rep1","DMSO_rep2","DMSO_rep3","THZ531_rep1","THZ531_rep2","THZ531_rep3","Hit4_rep1","Hit4_rep2","Hit4_rep3","Hit1_rep1","Hit1_rep2","Hit1_rep3","Hit2_rep1","Hit2_rep2","Hit2_rep3","Hit8_rep1","Hit8_rep2","Hit8_rep3")
        ,col=colors_MGs_18,ylab="Loess-normalized log2 cpm")
dev.off()






######################################################
# DE analysis, use lmFit to calculate coeffcients
######################################################

coldata=as.data.frame(read.table("coldata_MGs_RNA_18_samples",sep="\t", header=TRUE, row.names=1, dec="."))

f_ON = factor(as.character(coldata$Condition), levels = unique(coldata$Condition))
designAmp <- model.matrix(~0+f_ON)
colnames(designAmp) = unique(coldata$Condition)

fit <- lmFit(count_df, designAmp)

fit<-eBayes(fit)





########################################################################
# Differentially expressed genes between Treatment vs DMSO
########################################################################

colors_names=as.data.frame(read.table("colors_names.txt",sep="\t", header=TRUE, row.names=1, dec="."))


#treatments=c("THZ531","Hit1_2_5uM","Hit2_7uM","Hit8_3_5uM","Hit4_25uM")



for (treatment_name in c("Hit4_25uM"))
{
  condition_name=as.character(colors_names[treatment_name,3])
  color=as.character(colors_names[treatment_name,2])
  print(paste(treatment_name,condition_name))
  
  fit <- lmFit(count_df, designAmp)
  contrast <- makeContrasts(deg = Hit4_25uM-DMSO, levels = designAmp)
  print(contrast)
  lmfit <- contrasts.fit(fit, contrast)
  
  
  #computes the empirical Bayes statistics for differential analysis
  #it's like a t-test but the standard errors are moderated across genes
  fit <- eBayes(lmfit)
  
  
  # Volcano plot
  
  res <- topTable(fit,adjust="BH",number=nrow(count_df), coef = 'deg', sort.by="none")
  
  pdf(paste("VolcanoPlot_",treatment_name,"_vs_DMSO.highlighted.Adj_pvalues.pdf",sep=""),useDingbats=FALSE)
  par(family = "Arial")
  plot(res$logFC,-log10(res$adj.P.Val),pch=16,main=paste("", condition_name, "vs DMSO"),cex=0.35,
       xlab="Log2 Fold Change", ylab="-log10(adj. P-value)",ylim=c(0,11),xlim=c(-8,6)
       ,family="Arial",col="#B4B4B4")
  
  points(res[which(abs(res$logFC)>=2 & res$adj.P.Val<=0.05),]$logFC,-log10(res[which(abs(res$logFC)>=2 & res$adj.P.Val<=0.05),]$adj.P.Val)
         ,pch=16,cex=0.35,ylim=c(0,11),xlim=c(-8,6),col=color)
  dev.off()
  
  
  write.table(res, file=paste("Comparison_",treatment_name,"_vs_DMSO.HGNC.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)
  
}








###################################################
# Violin plots (log2FC values)
###################################################

# Load a file including log2FC values (from comparison treatment vs DMSO) per each treatment

LFC_rnk_merged = as.data.frame(read.table("LFC_rnk_merged.txt"),sep="\t", header=T, dec=".", fill = TRUE)



#DMSO: 153, 153, 153, #999999
#THZ: 222, 241, 245, #DEF1F5
#Hit 4: 241, 146, 34, #F19222
#Hit1: 4, 104, 145, #046891
#Hit2: 94, 195, 235, #5EC3EB
#Hit8: 0, 48, 72, #003047

MG_treatments_colors=c("#DEF1F5","#046891","#5EC3EB","#003047","#F19222")



pdf("violin_plots_Log2FC_all_genes.lmfit.pdf",5,8)
par(family = "Arial")
vioplot(LFC_rnk_merged[,2:6]
        ,col=MG_treatments_colors,ylab="Log2FC in gene expression (treatment vs DMSO)"
        ,rectCol="#B4B4B4"
        ,lineCol="#B4B4B4"
        ,border="black",family="Arial"
        ,side="both"
        ,range=1
        ,areaEqual=F
        #,plotCentre="line"
        ,ylim=c(-8,6)
        ,las=2
)
abline(h=0, col="black",lty="dashed")


# Add statistical significance determined by two-sided unpaired t-test applied on levels of gene expression in compound-treated cells (3 replicates) vs DMSO-treated cells (3 replicates)
G<-reads_hg38_ERCC_loess[gene_rows,]

t.test(apply(G[,1:3],1,mean),apply(G[,4:6],1,mean) )$p.value #THZ 1.712608e-293
t.test(apply(G[,1:3],1,mean),apply(G[,7:9],1,mean) )$p.value #Hit4 0.2439421
t.test(apply(G[,1:3],1,mean),apply(G[,10:12],1,mean) )$p.value #Hit1 6.922707e-238
t.test(apply(G[,1:3],1,mean),apply(G[,13:15],1,mean) )$p.value #Hit2 3.323458e-304
t.test(apply(G[,1:3],1,mean),apply(G[,16:18],1,mean) )$p.value #Hit8 2.400329e-32

text(x = c(1, 2,3,4,5), y = c(6,6), labels = c("p<0.0001","p<0.0001","p<0.0001","p<0.0001","p=0.24"),cex=0.8)

dev.off()





########################################################
# Waterfall plots for Log2FC values
########################################################

# Load a file including log2FC values (from comparison treatment vs DMSO) per each treatment, ordered based on THZ treatment

merged_ordered_rnk=as.data.frame(read.table("rankings_merged.THZ_ordered.txt",sep="\t", header=TRUE,  dec=".", fill = TRUE))


# Waterfall plots - separately per each sample

pdf("waterfall_plots.pdf",10,8)

par(family = "Arial")
par(mfrow=c(2,3))
barplot(rnk_THZ531[,2],ylim=c(-7,7),col="#DEF1F5",border=NA,main="THZ531",ylab="Log2FC in gene expression")
barplot(rnk_Hit1_2_5uM[,2],ylim=c(-7,7),col="#046891",border=NA,main="Hit1",ylab="Log2FC in gene expression")
barplot(rnk_Hit2_7uM[,2],ylim=c(-7,7),col="#5EC3EB",border=NA,main="Hit2",ylab="Log2FC in gene expression")
barplot(rnk_Hit8_3_5uM[,2],ylim=c(-7,7),col="#003047",border=NA,main="Hit8",ylab="Log2FC in gene expression")
barplot(rnk_Hit4_25uM[,2],ylim=c(-7,7),col="#F19222",border=NA,main="Hit4",ylab="Log2FC in gene expression")

dev.off()


# Waterfall plots - overlapping

pdf("waterfall_plots.allInOne.full.pdf",10,8)
par(family = "Arial")
barplot(rnk_Hit2_7uM[,2],ylim=c(-7,7),col=adjustcolor("#5EC3EB",alpha.f = 1),border=NA,ylab="Log2FC in gene expression",add=FALSE)
barplot(rnk_THZ531[,2],ylim=c(-7,7),col=adjustcolor("#DEF1F5",alpha.f = 1),border=NA,ylab="Log2FC in gene expression",add=TRUE)

barplot(rnk_Hit1_2_5uM[,2],ylim=c(-7,7),col=adjustcolor("#046891",alpha.f = 1),border=NA,ylab="Log2FC in gene expression",add=TRUE)
barplot(rnk_Hit8_3_5uM[,2],ylim=c(-7,7),col=adjustcolor("#003047",alpha.f = 1),border=NA,ylab="Log2FC in gene expression",add=TRUE)
barplot(rnk_Hit4_25uM[,2],ylim=c(-7,7),col=adjustcolor("#F19222",alpha.f = 1),border=NA,ylab="Log2FC in gene expression",add=TRUE)
legend("bottomright",pch=19,legend=c("Hit2","THZ531","Hit1","Hit8","Hit4"),col=c("#5EC3EB","#DEF1F5","#046891","#003047","#F19222"))
dev.off()


# Waterfall plots - overlapping, transparent

pdf("waterfall_plots.allInOne.transparent2.pdf",10,8)
par(family = "Arial")
barplot(rnk_Hit2_7uM[,2],ylim=c(-7,7),col=adjustcolor("#5EC3EB",alpha.f = 0.6),border=NA,ylab="Log2FC in gene expression",add=FALSE)
barplot(rnk_THZ531[,2],ylim=c(-7,7),col=adjustcolor("#DEF1F5",alpha.f = 0.8),border=NA,ylab="Log2FC in gene expression",add=TRUE)

barplot(rnk_Hit1_2_5uM[,2],ylim=c(-7,7),col=adjustcolor("#046891",alpha.f = 0.5),border=NA,ylab="Log2FC in gene expression",add=TRUE)
barplot(rnk_Hit8_3_5uM[,2],ylim=c(-7,7),col=adjustcolor("#003047",alpha.f = 0.5),border=NA,ylab="Log2FC in gene expression",add=TRUE)
barplot(rnk_Hit4_25uM[,2],ylim=c(-7,7),col=adjustcolor("#F19222",alpha.f = 0.4),border=NA,ylab="Log2FC in gene expression",add=TRUE)
legend("bottomright",pch=19,legend=c("Hit2","THZ531","Hit1","Hit8","Hit4"),col=c("#5EC3EB","#DEF1F5","#046891","#003047","#F19222"))
dev.off()








###########################################################
# Generate MDS plots and heatmap
###########################################################

normalizedMat <- count_df


###########


n_top=10000    # number of the most variable genes used for MDS

mds<-plotMDS(normalizedMat,top=n_top)

dev.off()


pdf("MDS_conditions_labels.pdf")

qplot(mds$x,mds$y) + theme_bw() + 
  theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black",size = 1.5),axis.text.x = element_text(size = 17, family="Arial"),axis.text.y = element_text(size = 17, family="Arial")) + 
  geom_point(shape = 21, colour = "black", fill = colors_MGs_18, size = 7, stroke = 1) +
  geom_text(label=c("","DMSO","","","THZ531","","","Hit4","","Hit1","","","Hit2","","","","","Hit8"),family="Arial",size=7,hjust=0.4, vjust=0,col="black") + labs(x="MDS Dimension 1", y="MDS Dimension 2") + 
  
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00, family="Arial")) +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90, family="Arial")) +
  theme(axis.ticks = element_line(size = 1.0)) +
  theme(axis.ticks.length = unit(+.30, "cm") ) 
dev.off()



pdf("MDS_conditions_NoLabels.pdf")
qplot(mds$x,mds$y) + theme_bw() + 
  theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black",size = 1.5),axis.text.x = element_text(size = 17, family="Arial"),axis.text.y = element_text(size = 17, family="Arial")) + 
  geom_point(shape = 21, colour = "black", fill = colors_MGs_18, size = 7, stroke = 1) +
  geom_text(label=c(" "),family="Arial",size=7,hjust=0.4, vjust=0,col="black") + labs(x="MDS Dimension 1", y="MDS Dimension 2")  +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00, family="Arial")) +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90, family="Arial")) +
  theme(axis.ticks = element_line(size = 1.0)) +
  theme(axis.ticks.length = unit(+.30, "cm") ) 
dev.off()






############################################################ 
# Select the most variable genes
############################################################

# Calculate standard deviation
normalizedMat_sd<-as.matrix(apply(normalizedMat,1,sd))

quantile(normalizedMat_sd,probs=seq(0,1,0.1))

# set standard deviation cut-off
SD_cutoff<-1

# select genes that have sd >= SD_cutoff
normalizedMat_varSD<-normalizedMat[normalizedMat_sd[,1]>=SD_cutoff,]

# save gene expression table with genes that have sd >= SD_cutoff
write.table(file=paste("MG_RNAseq_spiked_18samples_variable_regions_sd_above_",SD_cutoff,".txt",sep=""),normalizedMat_varSD,quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)



# Median centering genes

med.normalizedMat_varSD<-apply(normalizedMat_varSD,1,median)
normalizedMat_varSD_medianCenteredGenes<-sweep(normalizedMat_varSD,1,med.normalizedMat_varSD)

# save gene expression table with genes that have sd >= SD_cutoff that are median-centered
write.table(file=paste("RNAseq_spiked_18samples_variable_regions_sd_above_",SD_cutoff,"_MedianCenteredGenes.txt",sep=""),normalizedMat_varSD_medianCenteredGenes,quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE) 




############################################################ 
# Hierarchical clustering and heatmap
############################################################

MED_CenteredGenes=as.data.frame(read.table("RNAseq_spiked_18samples_variable_regions_sd_above_1_MedianCenteredGenes.txt",sep="\t", header=TRUE, dec=".", fill = TRUE))

set_breaks_color_scale=seq(-4, 4, 0.01)
length(set_breaks_color_scale)

my_palette <- colorRampPalette(c("#5757F9", "white", "#FF8081"))(n = length(set_breaks_color_scale)-1)

y=data.matrix(MED_CenteredGenes)[,]



c <- cor(t(y), method="spearman") 
d <- as.dist(1-c)
hr <- hclust(d, method = "average", members=NULL)

hc <- hclust(as.dist(1-cor(y, method="spearman")), method="average")
col.dd<-as.dendrogram(hc)
labels(col.dd)
colnames(y)
order.dendrogram(col.dd)
plot(reorder(col.dd, wts = order(order.dendrogram(col.dd))))
neworder = c("Hit4_rep1", "Hit4_rep2","Hit4_rep3","DMSO_rep1", "DMSO_rep2","DMSO_rep3","Hit1_rep1", "Hit1_rep2","Hit1_rep3","THZ531_rep1", "THZ531_rep2","THZ531_rep3","Hit8_rep1", "Hit8_rep2","Hit8_rep3")
col.dd.reordered  = reorder(col.dd, wts = order(match(neworder, colnames(y))))
plot(col.dd.reordered)
order.dendrogram(col.dd.reordered)







pdf("Heatmap3_SD1_median_Centered_genes.hierarch_clustering.pdf",14,8)
par(cex=0.1,family="Arial")

heatmap3(y,
         col = my_palette,
         Rowv=as.dendrogram(hr),
         Colv=col.dd.reordered,
         labRow=NA ,margins = c(9,1),cexCol=1.5,family="Arial"
         ,scale = "none"
         ,showRowDendro=FALSE, useRaster = FALSE,balanceColor=F
         #,legendfun=function() showLegend(legend = c("eee"),cex=0.1, col=("white"))
)

dev.off()










##################################################################################################################
# Venn diagrams for significantly differentially expressed genes.
#
# Genes with an absolute log2 fold change > 2 and an FDR-adjusted P-value < 0.05 were considered 
# as statistically significantly differentially expressed.
##################################################################################################################



setwd("comparisons/genesets/")

# Venn diagram for Significantly down-regulated genes

Hit1_DN<-read.table("DN/signig_genes_Hit1_2_5uM_vs_DMSO.lfc2_padj05.DN.gmx",sep="\t", header=T, dec=".", fill = F)
Hit2_DN<-read.table("DN/signig_genes_Hit2_7uM_vs_DMSO.lfc2_padj05.DN.gmx",sep="\t", header=T, dec=".", fill = F)
Hit4_DN<-read.table("DN/signig_genes_Hit4_25uM_vs_DMSO.lfc2_padj05.DN.gmx",sep="\t", header=T, dec=".", fill = F)
Hit8_DN<-read.table("DN/signig_genes_Hit8_3_5uM_vs_DMSO.lfc2_padj05.DN.gmx",sep="\t", header=T, dec=".", fill = F)
THZ531_DN<-read.table("DN/signig_genes_THZ531_vs_DMSO.lfc2_padj05.DN.gmx",sep="\t", header=T, dec=".", fill = F)

# With Hit4
gmt_DN<-c(as.list(Hit1_DN),as.list(Hit2_DN),as.list(Hit4_DN),as.list(Hit8_DN),as.list(THZ531_DN))

pdf("Venn_Down_regulated.ellipse.pdf")
par(family="Arial",cex=1.5,par(mfrow=c(1,1)))
venn(gmt_DN, ilab=TRUE, zcolor = c("#046891" , "#5EC3EB", "#F19222", "#003047", "#DEF1F5"),ellipse=F,opacity=0.3,cexil=0.7)

dev.off()




# Without Hit4
gmt_DN<-c(as.list(Hit1_DN),as.list(Hit2_DN),as.list(Hit8_DN),as.list(THZ531_DN))

pdf("Venn_Down_regulated.without_Hit4.ellipse.pdf")
par(family="Arial",cex=1.5,par(mfrow=c(1,1)))
venn(gmt_DN, ilab=TRUE, zcolor = c("#046891" , "#5EC3EB", "#003047", "#DEF1F5"),ellipse=T,opacity=0.3,cexil=0.7)

dev.off()





# Venn diagram for Significantly up-regulated genes

Hit1<-read.table("UP/signig_genes_Hit1_2_5uM_vs_DMSO.lfc2_padj05.UP.gmx",sep="\t", header=T, dec=".", fill = F)
Hit2<-read.table("UP/signig_genes_Hit2_7uM_vs_DMSO.lfc2_padj05.UP.gmx",sep="\t", header=T, dec=".", fill = F)
Hit4<-read.table("UP/signig_genes_Hit4_25uM_vs_DMSO.lfc2_padj05.UP.gmx",sep="\t", header=T, dec=".", fill = F)
Hit8<-read.table("UP/signig_genes_Hit8_3_5uM_vs_DMSO.lfc2_padj05.UP.gmx",sep="\t", header=T, dec=".", fill = F)
THZ531<-read.table("UP/signig_genes_THZ531_vs_DMSO.lfc2_padj05.UP.gmx",sep="\t", header=T, dec=".", fill = F)

# With Hit4
gmt_UP<-c(as.list(Hit1),as.list(Hit2),as.list(Hit4),as.list(Hit8),as.list(THZ531))

pdf("Venn_UP_regulated.ellipse.pdf")
par(family="Arial",cex=1.5,par(mfrow=c(1,1)))
venn(gmt_UP, ilab=TRUE, zcolor = c("#046891" , "#5EC3EB", "#F19222", "#003047", "#DEF1F5") ,ellipse=F,opacity=0.3,cexil=0.7)

dev.off()





# Without Hit4
gmt_UP<-c(as.list(Hit1),as.list(Hit2),as.list(Hit8),as.list(THZ531))

pdf("Venn_UP_regulated.without_Hit4.ellipse.pdf")
par(family="Arial",cex=1.5,par(mfrow=c(1,1)))
venn(gmt_UP, ilab=TRUE, zcolor = c("#046891" , "#5EC3EB", "#003047", "#DEF1F5"),ellipse=T,opacity=0.3,cexil=0.7)

dev.off()






