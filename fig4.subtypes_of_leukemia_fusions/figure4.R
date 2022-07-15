
####################################################
#####-Figure 4A----heatmap or PCA analysis-----#####
####################################################
#annotation
anno.col<-data.frame(cell.line = rep(c('GM12878','PBMC'),c(11,18)),
                     cell.type = c(rep('B lymphoblastoid',11),
                                   'T lymphocyte','B lymphocyte','T lymphocyte','T lymphocyte','T lymphocyte',
                                   'T lymphocyte','T lymphocyte','T lymphocyte','monocyte/neutrophil',
                                   'T lymphocyte','T lymphocyte','T lymphocyte','T lymphocyte','monocyte/neutrophil',
                                   'T lymphocyte','T lymphocyte','T lymphocyte','monocyte/neutrophil'))
anno.col<-anno.col[,c("cell.type","cell.line")]
row.names(anno.col)<-c(gm12878,pbmc)
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

pheatmap(dheatmap,cluster_cols = T,fontsize=7,
         annotation_col = anno.col,
         annotation_row=anno.row,
         fontsize=8,cexCol=2,
         cellwidth = 40, 
         cellheight = 8,
         angle_col = "45",
#         filename = "02_heatmap.pdf")

#---PCA using tSNE---#
library(Rtsne)
library(DMwR)
tsnemap<-t(dheatmap)
tsnemap<-centralImputation(tsnemap)   #replace NA values

set.seed(42)
tsneout <- Rtsne(tsnemap,dims=2,perplexity=5,theta=0.0,max_iter=1000)

tsneplot <- as.data.frame(tsneout$Y)
colnames(tsneplot) <- c("tSNE1","tSNE2")
tsneplot$labels<-cells[row.names(anno.col)]
tsneplot$cell.type<-anno.col$cell.type
tsneplot$cell.line<-anno.col$cell.line

fig4a<-
  ggplot(tsneplot, aes(x=tSNE1, y=tSNE2, col=cell.line,shape=cell.type))+
  geom_point(size=2.5)+
  theme(legend.position='bottom',
        legend.direction = "vertical")+
  stat_ellipse()

####################################################
#####-Figure 4B-fusion diff between AML and ALL-####
####################################################
#B-cell:GM12878,PBMC(2)
#T-cell:pbmc(1,4:8,10,12,13,15,16)
#monocyte:PBMC(9,14,18)

#-Figure.4B---AML vs ALL---#
#-AML:pbmc[c(9,14,18)] include GSM3271388,GSM3271398 and GSM3271406
#-ALL:B-ALL and T-ALL
aml.sig<-c()
for (i in c(1:58)){
  p.value='NA'
  tryCatch({
    p.value<-t.test(as.numeric(dfusions.eud[i,pbmc[c(9,14,18)]]),as.numeric(dfusions.eud[i,c(pbmc[c(1:8,10:13,15:17)],gm12878)]))$p.value
  },error=function(e){cat("error","\n")})
  aml.sig<-c(aml.sig,p.value)
}
row.names(aml.sig)<-paste(dfusions.eud$left,dfusions.eud$right,sep='-')
aml.sig<-data.frame(aml.sig,dfusions.eud$ave.eud)
aml.sig$fusion<-row.names(aml.sig)
colnames(aml.sig)<-c('-log10(p.value)','Average.distance')
aml.sig$`-log10(p.value)`<--log10(aml.sig$`-log10(p.value)`)
####
library(ggrepel)
fig4b1<-ggplot(aml.sig, aes(y =`-log10(p.value)`, x = Average.distance))+
  geom_point() +
  scale_color_manual(values = c("red", "grey","black")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  geom_text_repel(
    data = subset(aml.sig, fusion %in% c('KMT2A-ABI1','KMT2A-MLLT4','KMT2A-SORBS2')),
    aes(label =fusion),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

#-ETV6-ABL1   KMT2A-ABI1  KMT2A-MLLT4 KMT2A-SORBS2 
#-0.006827806  0.039932204  0.008373126  0.017206130

Cell.type<-rep(rep(c('monocyte','lymphocyte'),c(3,26)),3)
Fusion<-rep(c('KMT2A-ABI1','KMT2A-MLLT4','KMT2A-SORBS2'),c(29,29,29))
Distances<-dfusions.eud[c(9,38,50),c(pbmc[c(9,14,18)],pbmc[c(1:8,10:13,15:17)],gm12878)]
Distances<-as.numeric(t(as.matrix(Distances)))
amlfusion<-data.frame(Fusion,Cell.type,Distances)
rm(Cell.type,Fusion,Distances)

#--facet_grid several plots in one figure---#
fig4b2<-ggplot(amlfusion, aes(x=Cell.type, y=Distances)) + 
  geom_boxplot(aes(fill = Cell.type), alpha = 0.5,show.legend = F) + 
  facet_grid( ~ Fusion,scales = "free")


#-----0418---new fig.4B---merge ALL, AML and common fusion genes by boxplot---#
#-AML:pbmc[c(9,14,18)] include GSM3271388,GSM3271398 and GSM3271406
#-ALL:B-ALL and T-ALL


amlgene<-c(22,35,36,39,47,48) #6
allgene<-c(12,13,34) #3
comgene<-c(37,38,24,51,46) #5

# myelocyte:3
#lymphocyte:26


amlgm<-as.numeric(as.matrix(dfusions.eud[amlgene,pbmc[c(9,14,18)]]))
amlgl<-as.numeric(as.matrix(dfusions.eud[amlgene,c(pbmc[c(1:8,10:13,15:17)],gm12878)]))
allgm<-as.numeric(as.matrix(dfusions.eud[allgene,pbmc[c(9,14,18)]]))
allgl<-as.numeric(as.matrix(dfusions.eud[allgene,c(pbmc[c(1:8,10:13,15:17)],gm12878)]))
comgm<-as.numeric(as.matrix(dfusions.eud[comgene,pbmc[c(9,14,18)]]))
comgl<-as.numeric(as.matrix(dfusions.eud[comgene,c(pbmc[c(1:8,10:13,15:17)],gm12878)]))

amll.dif<-data.frame(
  Fusions<-c(rep(c('AML fusions','ALL fusions','Common fusions'),c(6*29,3*29,5*29))),
  Cell.types<-c(rep(rep(c('myelocyte','lymphocyte'),3),c(6*3,6*26,3*3,3*26,5*3,5*26))),
  Distances<-c(amlgm,amlgl,allgm,allgl,comgm,comgl)
)

new4b1<-ggplot(amll.dif, aes(x=Fusions, y=Distances)) + 
  geom_boxplot(aes(fill = Cell.types), alpha = 0.5,show.legend = T)+
  theme(legend.position = "top")
new4b1

wilcox.test(amlgm,amlgl) #ns, p-value = 0.8919
wilcox.test(allgl,allgm) #ns, p-value = 0.8727
wilcox.test(comgl,comgm) #ns, p-value = 0.6103

#---generate new fig.4b1
pdf("../../figures/figure.4/figure4b1.v2.pdf",width = 4.52,height = 3)
new4b1
dev.off()

###################################################################
#---KMT2A-AFF1，KMT2A-MLLT3（AF9） and KMT2A-MLLT10---#
dfusions.eud[dfusions.eud$right %in% c('AFF1','MLLT3','MLLT10'),c('left','right')]
#Row 12,37,35
Cell.type<-rep(rep(c('monocyte','lymphocyte'),c(3,26)),3)
Fusion<-rep(c('KMT2A-AFF1','KMT2A-MLLT3','KMT2A-MLLT10'),c(29,29,29))
Distances<-dfusions.eud[c(12,37,35),c(pbmc[c(9,14,18)],pbmc[c(1:8,10:13,15:17)],gm12878)]
Distances<-as.numeric(t(as.matrix(Distances)))
alfusion<-data.frame(Fusion,Cell.type,Distances)
rm(Cell.type,Fusion,Distances)

#-As sup figure-#-facet_grid several plots in one figure---#
sup4a<-
  ggplot(alfusion, aes(x=Cell.type, y=Distances)) + 
  geom_boxplot(aes(fill = Cell.type), alpha = 0.5,show.legend = F) + 
  facet_grid( ~ Fusion,scales = "free")

#---B-ALL vs T-ALL---#
all.sig<-c()
for (i in c(1:58)){
  p.value='NA'
  tryCatch({
    p.value<-t.test(as.numeric(dfusions.eud[i,gm12878]),as.numeric(dfusions.eud[i,c(pbmc[c(1:8,10:13,15:17)])]))$p.value
  },error=function(e){cat("error","\n")})
  all.sig<-c(all.sig,p.value)
}
names(all.sig)<-paste(dfusions.eud$left,dfusions.eud$right,sep=',')

all.sig.raw<-all.sig
all.sig<-data.frame(all.sig)
row.names(all.sig)<-paste(dfusions.eud$left,dfusions.eud$right,sep='-')
all.sig<-data.frame(all.sig,dfusions.eud$ave.eud)
all.sig$fusion<-row.names(all.sig)
colnames(all.sig)<-c('-log10(p.value)','Average.distance')
all.sig$`-log10(p.value)`<--log10(all.sig$`-log10(p.value)`)
colnames(all.sig)[3]<-'fusion'
all.sig<-all.sig[grepl('KMT2A',all.sig$fusion),]

library(ggrepel)
sup.fig4b<-ggplot(all.sig, aes(y =`-log10(p.value)`, x = Average.distance))+
  geom_point() +
  scale_color_manual(values = c("red", "grey","black")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  geom_text_repel(
    data = subset(all.sig, fusion %in% c('KMT2A-NCKIPSD','KMT2A-AFF1')),
    aes(label =fusion),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

#---Figure.4C and 4D---associations between the number of patients and distances----#
#-fig.4C
library("biomaRt") #retrive information for ensembl

mll=read.table("fig5.distances_and_occurences/all_MLL_fusions.txt",sep = "\t",header = T,stringsAsFactors = F)

mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
geneinfo<-
  getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "strand"),
        filters = "external_gene_name",
        values = mll$gene,
        mart=mart)
geneinfo<-geneinfo[!grepl('_',geneinfo$chromosome_name),]
setdiff(mll$gene,geneinfo$external_gene_name) #unmatched ID

#Unmatched ID: [1] "KNL1"         "FOXO3A"       "KIAS1524"     "AP2S2"        "CT45S2"       "LOC100131626"
#                   CASC5           FOXO3          NO             NO             CT45A2         NO
# remove KIAS1524, AP2S2 and LOC100131626
colnames(geneinfo)<-c('gene','chr','loc','strand')
mllpart<-merge(mll,geneinfo);rm(geneinfo)

dis=read.table("20kb.bin.bed",sep = "\t",stringsAsFactors = F)
colnames(dis)<-c('chr','start','end')
mllpart$Location1<-rep(104479,nrow(mllpart))

#right parter
mllpart$Location2<-
  apply(mllpart[,c("chr","loc")],1,function(x){
    any.chr<-as.character(x[1]); any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>dis$start & any.pos<dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )

#---Distances---#
mllpart.eud<-mllpart
for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(mllpart[,c("Location1","Location2")],1,function(x){
      ## 7:9 refer to x,y,z coordinates
      mat1<-as.numeric(mat.hic[x[1],7:9]); pat1<-as.numeric(pat.hic[x[1],7:9])   #two alleles of location 1
      mat2<-as.numeric(mat.hic[x[2],7:9]); pat2<-as.numeric(pat.hic[x[2],7:9])   #two alleles of location 2
      tmpd=min(c(
        sqrt((mat1[1]-mat2[1])^2+(mat1[2]-mat2[2])^2+(mat1[3]-mat2[3])^2),
        sqrt((mat1[1]-pat2[1])^2+(mat1[2]-pat2[2])^2+(mat1[3]-pat2[3])^2),
        sqrt((pat1[1]-mat2[1])^2+(pat1[2]-mat2[2])^2+(pat1[3]-mat2[3])^2),
        sqrt((pat1[1]-pat2[1])^2+(pat1[2]-pat2[2])^2+(pat1[3]-pat2[3])^2)),
        na.rm = T)
      tmpd[is.infinite(tmpd)]<-NA
      return(tmpd);rm(tmpd)
    }
    )
  mllpart.eud<-cbind(mllpart.eud,eu.d)
  rm(eu.d)
}
rm(i,mat.hic,pat.hic)

colnames(mllpart.eud)<-c(colnames(mllpart),samples)

#GM12878+PBMC
mllpart.eud$ave.eud<-apply(mllpart.eud[,c(gm12878,pbmc)],1,function(x){mean(x,na.rm = T)})
mllpart.eud$min.eud<-apply(mllpart.eud[,c(gm12878,pbmc)],1,function(x){ min(x,na.rm = T)})  #generate inf value
#mllpart.eud$min.eud[is.infinite(mllpart.eud$min.eud)]<-'NA'  #remove inf value
mllpart.eud$min.eud<-as.numeric(mllpart.eud$min.eud)
mllpart.eud$ave.eud<-as.numeric(mllpart.eud$ave.eud)

#-GM12878
mllpart.eud$ave.gm12878<-apply(mllpart.eud[,gm12878],1,function(x){mean(x,na.rm = T)})
mllpart.eud$min.gm12878<-apply(mllpart.eud[,gm12878],1,function(x){ min(x,na.rm = T)})  #generate inf value
#mllpart.eud$min.gm12878[is.infinite(mllpart.eud$min.gm12878)]<-'NA'  #remove inf value
mllpart.eud$min.gm12878<-as.numeric(mllpart.eud$min.gm12878)
mllpart.eud$ave.gm12878<-as.numeric(mllpart.eud$ave.gm12878)
#PBMC
mllpart.eud$ave.pbmc<-apply(mllpart.eud[,pbmc],1,function(x){mean(x,na.rm = T)})
mllpart.eud$min.pbmc<-apply(mllpart.eud[,pbmc],1,function(x){ min(x,na.rm = T)})  #generate inf value
#mllpart.eud$min.pbmc[is.infinite(mllpart.eud$min.pbmc)]<-'NA'  #remove inf value
mllpart.eud$min.pbmc<-as.numeric(mllpart.eud$min.pbmc)
mllpart.eud$ave.pbmc<-as.numeric(mllpart.eud$ave.pbmc)

#---Remove fusions on the same chromosome with MLL(chr11)---#
#mllpartall.eud<-mllpart.eud
mllpart.eud<-mllpart.eud[mllpart.eud$chr!=11,]
#52 in all
#32 in common with COSMIC


####################################################
#####-Figure 4C and 4D---cor of dis and occurence-----#####
####################################################
library(ggpubr)
fig4c<-ggscatter(mllpart.eud, x = "ave.eud", y = "patients",add = "reg.line",size = 1,
                 conf.int = TRUE,add.params = list(fill = "lightgray"),
                 xlab='Average Distance',ylab = 'Occurences of MLL Fusion',
                 title='    MLL fusions')+
  scale_y_log10(limits=c(1,1000))+
  stat_cor(method = "pearson")

####################################################
#####-Figure 4D-dis and occurence in SEC complex-###
####################################################
#fig.4D
sec<-data.frame(c('AFF1','MLLT3','MLLT1','MLLT10','ELL','MLLT6','AFF3','AFF4'),
                c(839,449,302,197,97,14,8,2),stringsAsFactors = F)
colnames(sec)<-c('Fusion','Frequency')
sec<-merge(dfusions.eud,sec,by.x ="right",by.y ="Fusion",all.x = F)

fig4d<-ggscatter(sec, x = "ave.eud", y = "Frequency",add = "reg.line",size = 1,
                 conf.int = TRUE,add.params = list(fill = "lightgray"),
                 xlab='Average Distance',ylab = 'Occurences',
                 title='    MLL and SEC complexes')+
  scale_y_log10(limits=c(1,1000))+
  stat_cor(method = "pearson")

library("cowplot")
fig4cd<-plot_grid(fig5a,fig5b,
          labels = c('C','D'),
          rel_widths = c(1,1),
          ncol = 2, nrow = 1)

fig4b<-plot_grid(fig4b1,fig4b2,ncol = 1,rel_heights = c(1.2,1))
fig4ab<-plot_grid(fig4a,fig4b,
          labels = c('A','B'),
          rel_widths = c(1,1.3),
          ncol = 2, nrow = 1)


pdf("../../figures/figure.4/figure4.pdf",width = 12,height = 10)
plot_grid(fig4ab,fig4cd,
          rel_heights = c(2,1),
          ncol = 1, nrow = 2)
dev.off()

pdf("../../figures/sup.figure.4/sup4.pdf",width = 8,height = 6)
plot_grid(sup4a,sup.fig4b,
          rel_heights = c(1,1),
          ncol = 1, nrow = 2)
dev.off()