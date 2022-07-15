library("magrittr")
library("scatterplot3d")
library(qqman)

#load distance matrix from SPRITE_less/
#load("06_allSVs_distances_on_spritecons.Rdata")

#load fusions info from('dfusions' dataframe)
#load("../02_diphic_dis/02_fusion_distances.Rdata")

#hg19 1M bin
bin1m<-read.table("hg19_1M_bin/hg19_1M_bin.bed",header = F,stringsAsFactors = F)
colnames(bin1m)<-c('chr','start','end')
bin1m$chr[bin1m$chr=='X']<-23 # chrX->23, match with dfusions

#SPRITE inter-chromosomal contacts
intercon <- read.table("GSE114242_human_inter_1Mb/human_inter_1Mb_nover2_1000000_iced.txt",sep="\t",
                       header=F,stringsAsFactors=F) 
intercon <- as.matrix(intercon)
bintag<-paste('chr',paste(bin1m$chr,bin1m$start/10^6,sep = '.'),'m',sep = '')
row.names(intercon)<-bintag
colnames(intercon)<-bintag

#-------fusion gene locations in dis dataframe--------# load from 02_diphic
#location each parters of SVs in spritecon(dataframe)
#left parter
dfusions$Location1<-
  apply(dfusions[,c("leftchr","leftloc")],1,function(x){
    any.chr<-x[1]; any.pos<-as.numeric(x[2])
    tmpnum<-which(bin1m$chr==any.chr & any.pos>bin1m$start & any.pos<bin1m$end)[1]
    return(tmpnum);rm(tmpnum) }
  )
#right parter
dfusions$Location2<-
  apply(dfusions[,c("rightchr","rightloc")],1,function(x){
    any.chr<-x[1]; any.pos<-as.numeric(x[2])
    tmpnum<-which(bin1m$chr==any.chr & any.pos>bin1m$start & any.pos<bin1m$end)[1]
    return(tmpnum);rm(tmpnum) }
  )

##--retrive the contacts between two parters--##
dfusions$spritecon<-apply(dfusions[,c("Location1","Location2")],1,function(x){intercon[x[1],x[2]]})

###------------Remove duplicated fusions in dfusions, such as KMT2A-AFF1 and AFF1-KMT2A------------###
cosmic=read.table("../02_diphic_dis/cosmic_all_fusions.txt",header = T,stringsAsFactors = F,sep = "\t")
cosmic<-cosmic[,c(1,5,9,10,11)]
dfusions<-merge(dfusions,cosmic)

###################################################
#######----pairs(fusion gene and random gene) as control---######
###################################################
nright<-sample(1:nrow(bin1m),58*2,replace=F)
ctrlfus<-as.data.frame(cbind(c(dfusions$Location1,dfusions$Location2),nright))
colnames(ctrlfus)<-c("Location1","Location2")
which(bin1m$chr[ctrlfus$Location1]==bin1m$chr[ctrlfus$Location2]) #18,95行相同染色体
# bin1m$chr[ctrlfus$Location1[c(18,95)]]
ctrlfus$Location2[c(18,95)]<-c(278,983) #random替代
#rm(nright)
ctrlfus$spritecon<-apply(ctrlfus,1,function(x){intercon[x[1],x[2]]})

###------------Significant difference--------------###
wilcox.test(dfusions$spritecon,ctrlfus$spritecon)
#p-value = 1.807e-05

#####----Figure 3------#######
#--Fig.3A.correlations between di-c and SPRITE
dfusions.spr<-dfusions
dfusions.hic<-dfusions.eud
dfusions.spr$right==dfusions.hic$right
dfusions.all<-cbind(dfusions.hic,dfusions.spr)[,c(1:26,36:38)]
dcor1<-dfusions.all[,c("spritecon","ave.eud")];colnames(dcor1)<-c('spritecon','eud')
dcor2<-dfusions.all[,c("spritecon","min.eud")];colnames(dcor2)<-c('spritecon','eud')
dcor<-rbind(dcor1,dcor2);dcor$type<-c(rep('average',nrow(dcor1)),rep('minimum',nrow(dcor2)))
rm(dcor1,dcor2,dfu)
dcor$fusion<-paste(dfusions.all$left,dfusions.all$right,sep="-")
dcor$type[dcor$type=='average']<-'Average'
dcor$type[dcor$type=='minimum']<-'Minimum'

p1<-ggplot(data = dcor, mapping = aes(x=spritecon, y = eud,col = type)) +
  geom_point(size=1) + 
  stat_smooth(method = 'lm') +
  theme(legend.position="top") +
  labs(x='SPRITE contacts',y='Distances')
#      geom_text(aes(label=fusion), size=2)

#--fig3b,dot plot--#
group<-c(rep('fusion',length(dfusions$spritecon)),rep('ctrl',length(dfusions$spritecon)))
dis<-c(dfusions$spritecon,sample(ctrlfus$spritecon,nrow(dfusions),replace =F))
cons<-data.frame(group,dis)  # data.fram[1-group + 2-distance]
cons<-cons[!(is.na(cons$dis)),]
rm(group,dis)
ylim1 = boxplot.stats(cons$dis)$stats[c(1, 5)]
p2<-ggplot(cons, aes(x=group, y=dis)) + 
  scale_color_manual(values=c("#999999", "#E69F00")) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  coord_cartesian(ylim = c(12,45))+
  labs(x='group',y='SPRITE contacts')

library(cowplot)
pdf("03_SPRITE.newctrl.pdf",10,6)
#--plot several figures in one picture--#
plot_grid(p1,p2,labels = LETTERS[1:2])
dev.off()
#save(dfusions.all,dcor,file = "dip-c_sprite_cor.rda")

######----focus on hotspot fusion gene (KMT2A)-----########
#KMT2A:1942,chr11.118m
#AFF1:781,chr4.87m
#MLLT1:2679,chr19.6m

#All intercon with 1942
exc11<-setdiff(1:ncol(intercon),which(bin1m$chr==11)) #exclude chr11
kmt2a<-intercon['chr11.118m',exc11]
kmt2a.sort<-sort(kmt2a,decreasing =T)
conPathways<-list()
conPathways$kmt2a<-bintag[dfusions$Location2[dfusions$left=='KMT2A']]

#GSEA-similar enrichment
library(fgsea)
library(data.table)
library(ggplot2)

#GSEA enrichment
set.seed(42)        #P.value=0.008994539
kmt2aRes<-fgsea(pathways = conPathways,stats =kmt2a.sort,minSize =1,maxSize = 100,nperm=10000)

#Plot enrichment plot
pdf("KMT2A_enrichment_plot.pdf",7,6)
plotEnrichment(conPathways[["kmt2a"]],kmt2a.sort)+labs(title = "KMT2A fusion partners")
dev.off()
