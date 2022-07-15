#######################################################
####---Figure.1A---Heatmap Af all Fusions-----------###
#######################################################
#fusionhmap<-dfusions.eud[,samples]
library(pheatmap)
#annotation
anno.col<-data.frame(cell.line = rep(c('GM12878','PBMC'),c(11,18)),
                     cell.type = c(rep('B lymphoblastoid',11),
                                   'T lymphocyte','B lymphocyte','T lymphocyte','T lymphocyte','T lymphocyte',
                                   'T lymphocyte','T lymphocyte','T lymphocyte','monocyte/neutrophil',
                                   'T lymphocyte','T lymphocyte','T lymphocyte','T lymphocyte','monocyte/neutrophil',
                                   'T lymphocyte','T lymphocyte','T lymphocyte','monocyte/neutrophil'))
anno.col<-anno.col[,c("cell.type","cell.line")]
row.names(anno.col)<-c(gm12878,pbmc)

anno.col2<-data.frame(cell.line = rep(c('GM12878','PBMC'),c(11,18)))
row.names(anno.col2)<-c(gm12878,pbmc)

anno.row<-data.frame(AML=rep('no',58),B.ALL=rep('no',58),T.ALL=rep('no',58))
anno.row<-lapply(anno.row,function(x) as.character(x))
anno.row$AML[grep('acute_myeloid_leukaemia',dfusions.eud$tumor)]<-'yes'
anno.row$B.ALL[grep('acute_lymphoblastic_B_cell',dfusions.eud$tumor)]<-'yes'
anno.row$T.ALL[grep('acute_lymphoblastic_T_cell',dfusions.eud$tumor)]<-'yes'
anno.row<-as.data.frame(anno.row)
row.names(anno.row)<-paste(dfusions.eud$left,dfusions.eud$right,sep = "-")
# Create dip-c matrix
dheatmap<-dfusions.eud[,c(gm12878,pbmc)]
row.names(dheatmap)<-paste(dfusions.eud$left,dfusions.eud$right,sep = "-")

#row:YlOrRd
library(pheatmap)
require("RColorBrewer")
#pdf("../../figures/figure.1A_2.pdf",width = 7,height = 10)
pheatmap(dheatmap,
         color = brewer.pal(n = 7, name ="YlOrRd"),
         annotation_col = anno.col,
         cluster_rows = T,cluster_cols = T,
         fontsize=6,cexCol=2,
         angle_col = "90",
         show_colnames=F,
         treeheight_col=20,
         treeheight_row = 30,
         cellwidth = 12, 
         cellheight = 6,
         scale = "column",
         )
#dev.off()
#filename = "../../figures/figure.1A.pdf"

#######################################################
####---Figure.1B---Significant difference-----------###
#######################################################
#---01---GM12878------#
group<-c(rep('fusion',length(dfusions.eud$ave.gm12878)),rep('ctrl',length(dfusions.eud$ave.gm12878)))
dis<-c(dfusions.eud$ave.gm12878,sample(ctrlsv.eud$ave.gm12878,nrow(dfusions.eud),replace =F))
ave<-data.frame(group,dis,stat=rep('Average',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dis)
####------Using Minimun Distance instead of average between samples-------##
group<-c(rep('fusion',length(dfusions.eud$min.gm12878)),rep('ctrl',length(dfusions.eud$min.gm12878)))
dis<-c(dfusions.eud$min.gm12878,sample(ctrlsv.eud$min.gm12878,nrow(dfusions.eud),replace =F))
mini<-data.frame(group,dis,stat=rep('Minimum',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dis)
sigfusion.gm12878<-rbind(ave,mini);rm(ave,mini)
colnames(sigfusion)<-c('Group','Distances','Stat')

#---02---PBMC------#
group<-c(rep('fusion',length(dfusions.eud$ave.pbmc)),rep('ctrl',length(dfusions.eud$ave.pbmc)))
dis<-c(dfusions.eud$ave.pbmc,sample(ctrlsv.eud$ave.pbmc,nrow(dfusions.eud),replace =F))
ave<-data.frame(group,dis,stat=rep('Average',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dis)
####------Using Minimun Distance instead of average between samples-------##
group<-c(rep('fusion',length(dfusions.eud$min.pbmc)),rep('ctrl',length(dfusions.eud$min.pbmc)))
dis<-c(dfusions.eud$min.pbmc,sample(ctrlsv.eud$min.pbmc,nrow(dfusions.eud),replace =F))
mini<-data.frame(group,dis,stat=rep('Minimum',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dis)
sigfusion.pbmc<-rbind(ave,mini);rm(ave,mini)
colnames(sigfusion)<-c('Group','Distances','Stat')
sigfusion<-rbind(sigfusion.gm12878,sigfusion.pbmc)
sigfusion$cell.type<-c(rep('GM12878',nrow(sigfusion.gm12878)),rep('PBMC',nrow(sigfusion.pbmc)))
colnames(sigfusion)<-c('Group','Distances','Stat','Cell.type')
#--facet_grid several plots in one figure---#
pdf("../../figures/figure.1B.pdf",height =5,width =7)
ggplot(sigfusion, aes(x=Group, y=Distances)) + 
  geom_boxplot(aes(fill = Group), alpha = 0.5,show.legend = F) + 
  facet_grid(Stat ~ Cell.type,scales = "free")
dev.off()

#----old controls:P-value----#
wilcox.test(ctrlsv.eud$ave.eud,dfusions.eud$ave.eud)$p.value
wilcox.test(ctrlsv.eud$min.eud,dfusions.eud$min.eud)$p.value
wilcox.test(ctrlsv.eud$ave.gm12878,dfusions.eud$ave.gm12878)$p.value
wilcox.test(ctrlsv.eud$min.gm12878,dfusions.eud$min.gm12878)$p.value
wilcox.test(ctrlsv.eud$ave.pbmc,dfusions.eud$ave.pbmc)$p.value
wilcox.test(ctrlsv.eud$min.pbmc,dfusions.eud$min.pbmc)$p.value


#######################################################
#########----Figure 1C.3D plot of BCR-ABL1----#########
#######################################################
library("magrittr")
library("scatterplot3d")
library(tidyverse) #数据框格式（Data.Frame）

# Three samples with close distances between BCR-ABL1---#
samples3d<-c('GSM3271366','GSM3271390','GSM3271406')
bcr<-dfusions[dfusions$left=='BCR' & dfusions$right=='ABL1',]$Location1
abl1<-dfusions[dfusions$left=='BCR' & dfusions$right=='ABL1',]$Location2

make3dph<-function(sample3dname){
  mat.hic=read.table(paste('20kb_bin/',sample3dname,'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  pat.hic=read.table(paste('20kb_bin/',sample3dname,'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  hic<-rbind(mat.hic,pat.hic)
  colnames(hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  hic$allele<-c(rep('maternal',nrow(mat.hic)),rep('paternal',nrow(pat.hic)))
  hic<-hic[,c("chr","start","end","x","y","z","allele")]
  hic$x<-as.numeric(hic$x); hic$y<-as.numeric(hic$y); hic$z<-as.numeric(hic$z)
  #---colour---#
  hic$colors<-hic$chr
  hic$colors[hic$chr=='22' & hic$allele=='maternal']<-"#F46D43"
  hic$colors[hic$chr=='22' & hic$allele=='paternal']<-"#FDAE61"
  hic$colors[hic$chr=='9' & hic$allele=='maternal']<-"#74ADD1"
  hic$colors[hic$chr=='9' & hic$allele=='paternal']<-"#4575B4"
  #---shape---#
  hic$pchs<-rep(20,nrow(hic))
  hic[c(bcr,bcr+nrow(mat.hic)),"pchs"]<-17 #triangle emphasis
  hic[c(abl1,abl1+nrow(mat.hic)),"pchs"]<-18 #菱形
  
  #---Point sizes---#
  hic$cexs<-rep(0.03,nrow(hic))
  #---BCR,ABL1:cex=2
  hic[c(bcr,bcr+nrow(mat.hic)),"cexs"]<-2
  hic[c(abl1,abl1+nrow(mat.hic)),"cexs"]<-2.5
  #---focus on chr9 and chr22---#
  hic<-hic[hic$chr %in% c(22,9),]
  return(hic)
  rm(hic,mat.hic,pat.hic)
}

pdf("../../figures/figure.1C_2.pdf",height =6,width =14*1.5)
layout(matrix(c(1,2,3),byrow = T,nrow = 1))
#--sample 1-GSM3271366,GM12878 Cell 15--#
s1hic<-make3dph(samples3d[1])
scatterplot3d(s1hic[,4:6],
              color = s1hic$colors,
              pch = s1hic$pchs,
              cex.symbols = s1hic$cexs,
              col.axis = "blue",col.grid ="lightblue",
              xlim = c(-60,60),ylim = c(-60,60),zlim = c(-60,60),
              box = T,angle=35,
              main = 'GM12878 Cell 15',xlab="",ylab = "",zlab = ""
)
#--sample 2---GSM3271390,PBMC Cell 10 #
s2hic<-make3dph(samples3d[2])
scatterplot3d(s2hic[,4:6],
              color = s2hic$colors,
              pch = s2hic$pchs,
              cex.symbols = s2hic$cexs,
              col.axis = "blue",col.grid ="lightblue",
              xlim = c(-60,60),ylim = c(-60,60),zlim = c(-60,60),
              box = T,angle=25,
              main='PBMC Cell 10',xlab="",ylab = "",zlab = ""
)
#--sample 3---GSM3271406,PBMC Cell 18#
s3hic<-make3dph(samples3d[3])
scatterplot3d(s3hic[,4:6],
              color = s3hic$colors,
              pch = s3hic$pchs,
              cex.symbols = s3hic$cexs,
              col.axis = "blue",col.grid ="lightblue",
              xlim = c(-60,60),ylim = c(-60,60),zlim = c(-60,60),
              box = T,angle=35,
              main='PBMC Cell 18',xlab="",ylab = "",zlab = ""
)
dev.off()


#----Plot in one figure----#
#---Using AI to put figure 1A,B,C together.



########################################
#-----------sup.figure1B.v3------------#
########################################
#-(1)-PBMC and GM12878 together#
fusRatio<-
  data.frame(
    Ratio=c(apply(dfusions.eud[,c(gm12878,pbmc)],1,function(x){x<-x[!is.na(x)];length(x[x<=15])})/length(c(gm12878,pbmc)),
            apply(ofusions.eud[,c(gm12878,pbmc)],1,function(x){x<-x[!is.na(x)];length(x[x<=15])})/length(c(gm12878,pbmc)),
            apply(ctrlsv.eud[,c(gm12878,pbmc)],1,function(x){x<-x[!is.na(x)];length(x[x<=15])})/length(c(gm12878,pbmc))),
    Group=rep(c('fusions','solid tumor fusions','ctrl'),c(nrow(dfusions.eud),nrow(ofusions.eud),nrow(ctrlsv.eud)))
  )
#fusRatio$Group<-factor(fusRatio$Group,levels = c('fusions','solid tumor fusions','ctrl'))
fusRatio$Group<-factor(fusRatio$Group,levels = c('ctrl','solid tumor fusions','fusions'))


#---Plot fig1B.v3---#
#pdf("../../figures/figure.1/figure.1B.v3.pdf",height =3,width =5)
ggplot(fusRatio, aes(x=Group, y=Ratio)) + 
  geom_boxplot(aes(fill = Group), alpha = 0.5,show.legend = F,outlier.shape = T)+
  coord_cartesian(ylim = c(0,0.35))+
  labs()
#dev.off()

#---Significance results---#
#Fusions~Solid fusions:p-value=0.0001021
#Fusions~Ctrl:p-value=7.659e-15
#Solid fusions~Ctrl:p-value=5.46e-08

#-(2)-PBMC and GM12878 separately#
fusRatio2<-data.frame(
  Ratio=c(apply(dfusions.eud[,gm12878],1,function(x){x<-x[!is.na(x)];length(x[x<=15])})/length(gm12878),
          apply(dfusions.eud[,pbmc],1,function(x){x<-x[!is.na(x)];length(x[x<=15])})/length(pbmc),
          apply(ofusions.eud[,gm12878],1,function(x){x<-x[!is.na(x)];length(x[x<=15])})/length(gm12878),
          apply(ofusions.eud[,pbmc],1,function(x){x<-x[!is.na(x)];length(x[x<=15])})/length(pbmc),
          apply(ctrlsv.eud[,gm12878],1,function(x){x<-x[!is.na(x)];length(x[x<=15])})/length(gm12878),
          apply(ctrlsv.eud[,pbmc],1,function(x){x<-x[!is.na(x)];length(x[x<=15])})/length(pbmc)),
  Group=rep(c('fusions','solid tumor fusions','ctrl'),c(nrow(dfusions.eud),nrow(ofusions.eud),nrow(ctrlsv.eud))*2),
  Cell.type=rep(c('GM12878','PBMC','GM12878','PBMC','GM12878','PBMC'),
                c(nrow(dfusions.eud),nrow(dfusions.eud),nrow(ofusions.eud),nrow(ofusions.eud),nrow(ctrlsv.eud),nrow(ctrlsv.eud)))
)
fusRatio2$Group<-factor(fusRatio2$Group,levels = c('fusions','solid tumor fusions','ctrl'))
