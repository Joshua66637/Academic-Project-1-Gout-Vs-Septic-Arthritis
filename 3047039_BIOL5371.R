#Assignment

#Loading Libraries
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(ggfortify)

#PART 1: Parsing Files

#Reading Input Files
anno = read.table("C:/DATA/Glasgow/Academics/Statsistics for Bioinformatics/3 Labs/Data for Report Assessment-20241003/Annotations.csv", header=TRUE, row.names=1, sep="\t")
deGoutHc = read.table("C:/DATA/Glasgow/Academics/Statsistics for Bioinformatics/3 Labs/Data for Report Assessment-20241003/DE_GOUT_vs_HC.csv", header=TRUE, row.names=1, sep="\t")
deSaHc = read.table("C:/DATA/Glasgow/Academics/Statsistics for Bioinformatics/3 Labs/Data for Report Assessment-20241003/DE_SA_vs_HC.csv", header=TRUE, row.names=1, sep="\t")
sampleInfo = read.table("C:/DATA/Glasgow/Academics/Statsistics for Bioinformatics/3 Labs/Data for Report Assessment-20241003/Sample_Information.csv", header=TRUE, row.names=1, sep="\t")
exp =data.frame(t ( read.table("C:/DATA/Glasgow/Academics/Statsistics for Bioinformatics/3 Labs/Data for Report Assessment-20241003/Expression_Table.csv", header=TRUE, row.names=1, sep="\t")))

#merging the two differential gene expression files
deMerged = merge(deGoutHc, deSaHc, by=0)
names(deMerged) = c("GeneID", "log2fold_G","p_G","padj_G", "log2fold_SA","p_SA","padj_SA")

#Adding the gene id column
deGoutHc$geneID= rownames(deGoutHc)
deGoutHc= deGoutHc[,c(length(deGoutHc), 1:length(deGoutHc)-1)]
deSaHc$geneID= rownames(deSaHc)
deSaHc= deSaHc[,c(length(deSaHc), 1:length(deSaHc)-1)]

#Merging sample information with gene expression information and converting sex and sample groups to factors
samExp= merge(sampleInfo, exp, by.x=0, by.y =0)
samExp$SAMPLE_GROUP= as.factor(samExp$SAMPLE_GROUP)
samExp$SEX= as.factor(samExp$SEX)
colnames(samExp)[1]= "Sample_Name"
colnames(samExp)[2]= "Samples"

#sorting by most significant genes
deGoutHc= deGoutHc[order(deGoutHc$p.adj),]
deSaHc= deSaHc[order(deSaHc$p.adj),]

#Get the significant genes based on adjusted p-values and fold change
deGoutHcSig = subset(deGoutHc, deGoutHc$p.adj < 0.05 & abs(deGoutHc$log2Fold)>1)
deSaHcSig = subset(deSaHc, deSaHc$p.adj < 0.05 & abs(deGoutHc$log2Fold)>1)
deMergedSig= subset(deMerged, padj_G < 0.05 & padj_SA < 0.05 & abs(deMerged$log2fold_G )>1 & abs(deMerged$log2fold_SA )>1)

#Part 2: Summary stats for clinical measurements
for(i in seq(1,27,9))
{
  print(summary(samExp[i:(i+8),2:4]))
}

#lm of neutrophils vs sex and disease status
SexModel= glm(samExp$NEUTROPHILS ~ samExp$SEX, family = "inverse.gaussian")
GroupModel= glm(samExp$NEUTROPHILS ~ samExp$Samples,family = "inverse.gaussian")
summary(SexModel)
summary(GroupModel)
anova(SexModel)
anova(GroupModel)
hist(rstandard(SexModel))
hist(rstandard(GroupModel))

#Boxplot of neutrophils vs sex segregated by Sample Groups
plot= ggplot(samExp, aes(y= samExp$NEUTROPHILS, x= samExp$SEX))+

  geom_boxplot(aes(colour= Samples, fill=Samples), linewidth= 0.5, alpha= 0.1)+

  stat_summary(fun=median, geom="point",  size=2)+
  stat_summary(aes(label= paste ("Median\n", after_stat(round(y, 2)),"\n\n\n")), fun=median, geom="text",  size=3)+

  xlab("Sex")+
  ylab("Neutrophil Count")+
  ggtitle("Neutrophil Count vs Sex faceted by Sample Groups")+

  ylim(0,20)+

  facet_wrap (~samExp$Samples, nrow = 1, ncol=3)+
  theme_economist()
plot

#Part three: Significant overlap between SA and Gout

#Hypergeometric Test
gout = nrow(deGoutHcSig)
sa = nrow(deSaHcSig)
overlap = nrow(deMergedSig)

saOnly = sa - overlap
goutOnly = gout - overlap
total = nrow(deMerged)

hyper=phyper(overlap-1, sa, total-sa, gout,lower.tail= FALSE)
print(hyper)

#PCA
PCA= prcomp(exp)
PC=summary(PCA)
pca_coordinates = data.frame(PCA$x)
PlotPCA = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = samExp$Samples)) +
          geom_point()+

          xlab(paste("PC1 (", (PC$importance[2]*100), "%)"))+
          ylab(paste("PC2 (", (PC$importance[5]*100), "%)"))+
          ggtitle("Principal Component Analysis")+

          scale_color_discrete(name="Samples")
PlotPCA

#glm for are genes expression vs sample group, sex and neutrophil data

#for Gout
for (gene in seq(1:5))
{
  fit1 = glm((samExp[,rownames(deGoutHcSig)[gene]]+10^-100) ~ samExp$Samples + samExp$SEX + samExp$NEUTROPHILS, family = "inverse.gaussian")#added 10^-100 to deal with 0s
  print(anova(fit1))
  hist(rstandard(fit1))
}

#For SA
for (gene in seq(1:5))
{
  fit1 = glm(samExp[,rownames(deSaHcSig)[gene]]+10^-100 ~ samExp$Samples + samExp$SEX + samExp$NEUTROPHILS, family = "inverse.gaussian")#added 10^-100 to deal with 0s
  print(anova(fit1))
  hist(rstandard(fit1))
}

#Boxplot of top 5 significant gout related genes vs sample groups
GoutSg= anno[rownames(deGoutHcSig)[1:5],1]

for (gene in seq(1:5))
{
  plot= ggplot(samExp, aes(y= log(samExp[,rownames(deGoutHcSig)[gene]]+1), x= samExp$Samples))+

    geom_boxplot(aes(colour=Samples, fill= Samples), linewidth= 0.5, alpha= 0.2)+

    stat_summary(fun.y=median, geom="point",  size=2)+
    stat_summary(aes(label= paste ("Median\n", after_stat(round(y, 2)),"\n\n\n")), fun=median, geom="text",  size=5)+

    xlab("Groups")+
    ylab(paste("Log( ", anno [rownames(deGoutHcSig)[gene],1], ") Expression"))+
    ggtitle(paste("Log( ", anno [rownames(deGoutHcSig)[gene],1], ") Expression vs Sample Groups"))+

    theme_economist()
  print(plot)
}

#Boxplot of top 5 significant SA related genes vs sample groups
SaSg= anno[rownames(deSaHcSig)[1:5],1]

for (gene in seq(1:5))
{
  plot= ggplot(samExp, aes(y= log(samExp[,rownames(deSaHcSig)[gene]]+1), x= samExp$Samples))+

    geom_boxplot(aes(colour=Samples, fill= Samples), linewidth= 0.5, alpha= 0.2)+

    stat_summary(fun.y=median, geom="point",  size=2)+
    stat_summary(aes(label= paste ("Median\n", after_stat(round(y, 2)),"\n\n\n")), fun=median, geom="text",  size=5)+

    xlab("Groups")+
    ylab(paste("Log( ", anno[rownames(deSaHcSig)[gene],1], ") Expression"))+
    ggtitle(paste("Log( ", anno[rownames(deSaHcSig)[gene],1], ") Expression vs Sample Groups"))+

    theme_economist()
  print(plot)
}

#Box plot of all 8 significant genes common to both diseases vs sample groups
BothSg= anno[deMergedSig[,1],1]

for (gene in seq(1:8))
{
  plot= ggplot(samExp, aes(y= log(samExp[,deMergedSig$GeneID[gene]]+1), x= samExp$Samples))+

    geom_boxplot(aes(colour=Samples, fill= Samples), linewidth= 0.5, alpha= 0.2)+

    stat_summary(fun.y=median, geom="point",  size=2)+
    stat_summary(aes(label= paste ("Median\n", after_stat(round(y, 2)),"\n\n\n")), fun=median, geom="text",  size=3)+

    xlab("Groups")+
    ylab(paste("Log( ", anno[deMergedSig$GeneID[gene],1], ") Expression"))+
    ggtitle(paste("Log( ", anno[deMergedSig$GeneID[gene],1], ") Expression vs Sample Groups"))+

    theme_economist()
  print(plot)
}


