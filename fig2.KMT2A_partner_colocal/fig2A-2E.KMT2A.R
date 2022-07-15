library(pdist) #distances between two matrix
library(fgsea) #GSEA
library(data.table) #read large data.frame
library(plyr)
library(ggplot2)
library("reshape2") #melt

#---KMT2A information---#
kmt2a<-as.character(unique(dfusions[dfusions$left=='KMT2A',c("leftchr","leftloc","Location1")]))

bin20k=read.table("20kb.bin.bed",header=F,stringsAsFactors=FALSE)
bin20k<-bin20k[,c(1:2)]
colnames(bin20k)<-c('rightchr','rightloc')
bin20k$Location2<-1:nrow(bin20k)
#---KMT2A genome-wide distances---#
kmt2adic<-matrix(rep(kmt2a,nrow(bin20k)),nrow = nrow(bin20k),byrow = T)
kmt2adic<-cbind(kmt2adic,bin20k)
colnames(kmt2adic)<-c('leftchr','leftloc','Location1','rightchr','rightloc','Location2')

kmt2adic.eud<-kmt2adic

for (i in 1:length(samples)){
  matc=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  matc<-matc[,7:9]
  kmt2am<-matc[kmt2a[3],]
  #colnames(matc) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  patc=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  patc<-patc[,7:9]
  kmt2ap<-patc[kmt2a[3],c("V7","V8","V9")]
  
  matc1<-matc;matc1[matc1$V7=='.',]<-kmt2am
  matc2<-matc;matc2[matc2$V7=='.',]<-kmt2ap
  patc1<-patc;patc1[patc1$V7=='.',]<-kmt2am
  patc2<-patc;patc2[patc2$V7=='.',]<-kmt2ap
 
  eu.d<-cbind(
    as.matrix(pdist(matc1,kmt2am)),
    as.matrix(pdist(matc2,kmt2ap)),
    as.matrix(pdist(patc1,kmt2am)),
    as.matrix(pdist(patc2,kmt2ap)))
  
  rm(matc,matc1,matc2,patc,patc1,patc2)
  
  eu.d<-apply(eu.d,1,function(x){min(x[x>0])})
  eu.d[is.infinite(eu.d)]<-NA
  kmt2adic.eud<-cbind(kmt2adic.eud,eu.d)
  
  rm(eu.d,matc,patc)
}
rm(i)

colnames(kmt2adic.eud)<-c(colnames(kmt2adic),samples)

#GM12878+PBMC
kmt2adic.eud$ave.eud<-apply(kmt2adic.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{mean(x,na.rm = T)}})
kmt2adic.eud$min.eud<-apply(kmt2adic.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{ min(x,na.rm = T)}})
kmt2adic.eud$min.eud<-as.numeric(kmt2adic.eud$min.eud)
kmt2adic.eud$ave.eud<-as.numeric(kmt2adic.eud$ave.eud)

#-GM12878
kmt2adic.eud$ave.gm12878<-apply(kmt2adic.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{mean(x,na.rm = T)}})
kmt2adic.eud$min.gm12878<-apply(kmt2adic.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{ min(x,na.rm = T)}})
kmt2adic.eud$min.gm12878<-as.numeric(kmt2adic.eud$min.gm12878)
kmt2adic.eud$ave.gm12878<-as.numeric(kmt2adic.eud$ave.gm12878)
#PBMC
kmt2adic.eud$ave.pbmc<-apply(kmt2adic.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{mean(x,na.rm = T)}})
kmt2adic.eud$min.pbmc<-apply(kmt2adic.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{ min(x,na.rm = T)}})
kmt2adic.eud$min.pbmc<-as.numeric(kmt2adic.eud$min.pbmc)
kmt2adic.eud$ave.pbmc<-as.numeric(kmt2adic.eud$ave.pbmc)

#Remove loci in chromosome 11(KMT2A)
kmt2adic.eud<-kmt2adic.eud[kmt2adic.eud$leftchr!=kmt2adic.eud$rightchr,] 

#-----GSEA enrichment significance----#
#---Pathway---#
conPathways<-list()
conPathways$kmt2a<-dfusions$Location2[dfusions$left=='KMT2A']
kmt2aGsea<-c()
colgsea<-c(gm12878,pbmc,"ave.eud","min.eud","ave.gm12878","min.gm12878","ave.pbmc","min.pbmc")
for (i in 1:length(colgsea)){
  kmt2a.sort<-kmt2adic.eud[,colgsea[i]]
  names(kmt2a.sort)<-kmt2adic.eud$Location2
  kmt2a.sort<-sort(kmt2a.sort,decreasing =F)
  #GSEA enrichment
  #set.seed(42)        #0.008994539
  kmt2aRes<-fgsea(pathways = conPathways,stats =kmt2a.sort,minSize =1,maxSize = 100,nperm=10000)
  #kmt2aGsea<-c(kmt2aGsea,kmt2aRes$pval)
  #Plot enrichment plot
  #pdf("KMT2A_enrichment_plot.pdf",7,4)
  plotEnrichment(conPathways[["kmt2a"]],100-kmt2a.sort)+labs(title = paste('KMT2A fusion partners(',colgsea[i],')',sep = ""))
  #dev.off()
}
#names(kmt2aGsea)<-colgsea

#-----GSEA enrichment leading Edge(最后分析)----#
kmt2aEdge<-numeric()
for (i in 1:length(colgsea[1:29])){
  kmt2a.sort<-kmt2adic.eud[,colgsea[i]]
  names(kmt2a.sort)<-kmt2adic.eud$Location2
  kmt2a.sort<-sort(kmt2a.sort,decreasing =F)
  kmt2aRes<-fgsea(pathways = conPathways,stats =kmt2a.sort,minSize =1,maxSize = 100,nperm=1000)
  kmt2aEdge<-c(kmt2aEdge,as.numeric(unlist(kmt2aRes$leadingEdge)))
  rm(kmt2aRes)
}
topEdge<-sort(table(kmt2aEdge),decreasing = T)
topEdge<-topEdge[topEdge>15]
topEdge.eud<-kmt2adic.eud[kmt2adic.eud$Location2 %in% names(topEdge),]


####################################################
#####-Figure 3C-----P.values distributions-----#####
####################################################
kmt2aGsea2<-kmt2aGsea
kmt2aGsea2<-kmt2aGsea2[c(gm12878,pbmc,"ave.eud","min.eud")]
kmt2aGsea2<-as.data.frame(kmt2aGsea2)
kmt2aGsea2$sample<-c(cells[row.names(kmt2aGsea2)[1:29]],'Average','Minimum')
colnames(kmt2aGsea2)<-c('-log10(pvalue)','sample')
kmt2aGsea2$`-log10(pvalue)`<-(-log10(kmt2aGsea2$`-log10(pvalue)`))
kmt2aGsea2$sample<-factor(kmt2aGsea2$sample,levels = kmt2aGsea2$sample)
kmt2aGsea2$Size<-factor(c(rep('single-cell',29),'Average','Minimum'))
library(ggplot2)
fig3c<-ggplot(data =kmt2aGsea2 ) + 
  geom_point(aes(x = sample, y = `-log10(pvalue)`,col = Size))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 6.5))+
  geom_hline(aes(yintercept=-log10(0.05)),linetype="dashed")

####################################################
#####------Figure 3D-3F----Three GSEA plots----#####
####################################################
#Result:
#-01, GSEA is not significant at single cell level(only 2/29 with P.value < 0.05).
#-02, GSEA is significant when taking all cells into consideration(ave.eud and min.eud).
#-3 GSEA figures: GSM3271348, ave.eud and min.eud.

#---GSM3271348(GM12878 Cell 2)
kmt2a.sort<-kmt2adic.eud[,'GSM3271348'];names(kmt2a.sort)<-kmt2adic.eud$Location2;kmt2a.sort<-sort(kmt2a.sort,decreasing =F)
gm2<-plotEnrichment(conPathways[["kmt2a"]],100-kmt2a.sort)+labs(title = paste('KMT2A(GM12878 Cell 2), P=0.820',sep = ""))
rm(kmt2a.sort)

#---Average
kmt2a.sort<-kmt2adic.eud[,"ave.eud"];names(kmt2a.sort)<-kmt2adic.eud$Location2;kmt2a.sort<-sort(kmt2a.sort,decreasing =F)
avecell<-plotEnrichment(conPathways[["kmt2a"]],100-kmt2a.sort)+labs(title = paste('KMT2A(Average), P=0.023',sep = ""))
rm(kmt2a.sort)

#---Minimum
kmt2a.sort<-kmt2adic.eud[,"min.eud"];names(kmt2a.sort)<-kmt2adic.eud$Location2;kmt2a.sort<-sort(kmt2a.sort,decreasing =F)
mincell<-plotEnrichment(conPathways[["kmt2a"]],100-kmt2a.sort)+labs(title = paste('KMT2A(Minimum), P=0.013',sep = ""))
rm(kmt2a.sort)

#pdf("../../figures/figure.3/figure3C-F.pdf",width = 16,height = 4)
library("cowplot")
plot_grid(fig3c,gm2,avecell,mincell,
          labels = c('C','D','E','F'),
          rel_widths = c(1.5,1,1,1),
          ncol = 4, nrow = 1)
#dev.off()

###########################################################
###--Figure 3I--verified by 2019,molecular cell paper--####
###########################################################
#vfusion<-as.data.frame(t(dheatmap[c('KMT2A-MLLT1','KMT2A-AFF1','KMT2A-MLLT4','KMT2A-MLLT3'),]))

#--Figure 3I(1)----take 'KMT2A-MLLT1','KMT2A-AFF1','KMT2A-MLLT4','KMT2A-MLLT3' analysis---#
vfusion<-as.data.frame(t(dheatmap[c('KMT2A-MLLT1','KMT2A-AFF1','KMT2A-MLLT4','KMT2A-MLLT3'),]))
vfusion$Samples<-row.names(vfusion)
vfusion<-vfusion[!is.na(vfusion$`KMT2A-MLLT1`) & !is.na(vfusion$`KMT2A-MLLT4`),]

vfusion<-melt(vfusion,id.vars = "Samples")
colnames(vfusion)<-c('Samples','Fusions','Distances')
#----boxplot---#
vf1<-ggplot(vfusion, aes(x=Fusions, y=Distances)) + 
  geom_boxplot(aes(fill = Fusions), alpha = 0.5,show.legend = F) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  geom_jitter(shape=16,size=1,position=position_jitter(0.2))

#-----方差分析：4种融合基因比较-----#
#---符合正太分布---#
ad.test(vfusion$Distances[vfusion$Fusions=='KMT2A-MLLT1'])
ad.test(vfusion$Distances[vfusion$Fusions=='KMT2A-AFF1'])
#---符合方差齐性---#
bartlett.test(Distances~Fusions,data = vfusion)
#---方差分析
summary(aov(formula = Distances~Fusions,data = vfusion))
#P=0.0107

summary(aov(formula = Distances~Fusions,data = vfusion[!vfusion$Fusions %in% 'KMT2A-MLLT1',]))
#P=0.478

#---Figure 3I(2)-density plot---#
mu <- ddply(vfusion, "Fusions", summarise, grp.mean=mean(Distances)) # Add mean lines
vf2<-ggplot(vfusion, aes(x=Distances, color=Fusions)) +
  geom_density()+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Fusions),linetype="dashed")

#---Figure 3I(3)---Focus on KMT2A-MLLT1（closest） and KMT2A-MLLT3(furthest)
v3fusion<-as.data.frame(t(dheatmap[c('KMT2A-MLLT1', 'KMT2A-MLLT3'),]))
v3fusion$Samples<-row.names(v3fusion)
v3fusion<-v3fusion[!is.na(v3fusion$`KMT2A-MLLT1`),]

v3fusion<-melt(v3fusion,id.vars = "Samples")
colnames(v3fusion)<-c('Samples','Fusions','Distances')
#----boxplot---#
vf3<-ggplot(v3fusion, aes(x=Fusions, y=Distances)) + 
  geom_boxplot(aes(fill = Fusions), alpha = 0.5,show.legend = F) + 
  geom_line(aes(group = Samples),alpha = 0.3)+
  geom_jitter(shape=16, position=position_jitter(0))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#---Figure 3I(4)-density plot---#
mu3 <- ddply(v3fusion, "Fusions", summarise, grp.mean=mean(Distances)) # Add mean lines
vf4<-ggplot(v3fusion, aes(x=Distances, color=Fusions)) +
  geom_density()+
  geom_vline(data=mu3, aes(xintercept=grp.mean, color=Fusions),linetype="dashed")

#pdf("../../figures/sup.figure.2/Supplementary Figure 2. Comparisons of 4 common KMT2A partners.gb.pdf",width = 16,height = 4)
library("cowplot")
plot_grid(vf1,vf2,vf3,
          labels = c('A','B','C'),
          rel_widths = c(1,1.5,1),
          ncol = 3, nrow = 1)
#dev.off()

#---Empty figure---#
void<-ggplot()+geom_point()+theme_bw()+theme_void()

#---Make figure 3.
#pdf("../../figures/figure.3/figure3A-I.pdf",width = 16,height = 12)
library("cowplot")

a_b<-plot_grid(void,fig3c,labels =c('A','B'),rel_widths = c(1,1.5))
c_h<-plot_grid(gm2,avecell,mincell,
               vf1,vf2,vf3,
               labels = c('C','D','E','F','H','I'),
               ncol = 3, nrow = 2)

plot_grid(a_b,c_h,rel_heights = c(1,2),ncol = 1)
#dev.off()

#---Supplementary Figure 3. Distances of 14 MLL fusions---#
mllt<-mllpart.eud[mllpart.eud$patients>4,]
row.names(mllt)<-paste('KMT2A',mllt$gene,sep = '-')
mllt<-mllt[,c(gm12878,pbmc)]
mllt<-mllt[order(apply(mllt,1,median)),]

mllrank<-matrix(nrow = 14,ncol = 14)
for (i in 1:nrow(mllrank)){
  pval<-c()
  for (j in 1:nrow(mllrank)){pval<-c(pval,wilcox.test(as.numeric(mllt[i,]),as.numeric(mllt[j,]),paired = T)$p.value)}
  mllrank[i,]<-pval
}
row.names(mllrank)<-row.names(mllt)
colnames(mllrank)<-row.names(mllt)

library("reshape2") #melt
library(ggplot2)
mlltPlot<-melt(t(mllt))
colnames(mlltPlot)<-c('Sample','Fusion','Distance')
mlltPlot$Fusion<-factor(mlltPlot$Fusion,levels = row.names(mllt))

#pdf("../../figures/sup.figure3.pdf",width = 16,height = 8)
library("cowplot")
ggplot(data = mlltPlot)+
  geom_boxplot(aes(x=Fusion,y=Distance,fill=Fusion),show.legend = F,alpha=0.5)+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
#dev.off()
#geom_hline(aes(yintercept=25),linetype="dashed")+
#  geom_hline(aes(yintercept=30),linetype="dashed")

#--First rank: KMT2A-SEPT5, KMT2A-MLLT1, KMT2A-ELL, KMT2A-MLLT11
#--Secound rank:KMT2A-SEPT9, KMT2A-EPS15 and KMT2A-AFF1. 
#--third KMT2A-MLLT4, KMT2A-MLLT10, KMT2A-AFF3, KMT2A-MLLT6, KMT2A-TET1, KMT2A-MLLT3, KMT2A-SEPT6.
