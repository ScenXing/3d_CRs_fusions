#-Hi-C matrix
ko=read.table("MCF7_shRUNX1_hg19_1Mb.matrix.txt",sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
gfp=read.table("MCF7_shGFP_hg19_1Mb.matrix.txt",sep = "\t",stringsAsFactors = F,header = T,row.names = 1)

bin1m=read.csv("bin1m.txt",header = F,stringsAsFactors = F)
bin1m$num=seq(1:nrow(bin1m))
colnames(bin1m)[1:3]=c('chr','start','end')
bin1m$name<-paste(bin1m$chr,ceiling(bin1m$end/10^6),sep = '_')

#Change row and column names
row.names(ko)=bin1m$name;colnames(ko)=bin1m$name;
row.names(gfp)=bin1m$name;colnames(gfp)=bin1m$name;
ko[ko=='NaN']=NA
gfp[gfp=='NaN']=NA

#Load 02_fusion_distances.rda
#--------Annotations of crgenes------------#
crgenesko=crgenes

#Locations of genes
crgenesko$loc.cr<-
  apply(crgenesko[,c("chr","start")],1,function(x){
    any.chr<-as.numeric(x[1]); any.pos<-as.numeric(x[2])
    tmpnum<-which(bin1m$chr==any.chr & any.pos>=bin1m$start & any.pos<=bin1m$end)[1]
    return(tmpnum);rm(tmpnum) }
  )
#add fusion partners
crgenesko$loc.partner<-
  apply(crgenesko[,c("Chromosome","chromStart")],1,function(x){
    any.chr<-x[1]; any.pos<-as.numeric(x[2])
    tmpnum<-which(bin1m$chr==any.chr & any.pos>=bin1m$start & any.pos<=bin1m$end)[1]
    return(tmpnum);rm(tmpnum) }
  )
crgenesko$loc.kmt2a<-rep(1942,nrow(crgenesko)) #Add KMT2A

#CRGs相关的分析（对应Fig5B）表明KO连接数更多
#crcombineko<-crcombine
library(hash)
#location
kotmp<-data.frame(gene=c(crgenesko$fusion,crgenesko$gene,'KMT2A'),loc=c(crgenesko$loc.partner,crgenesko$loc.cr,1942))
kotmp<-unique(kotmp);kotmp<-kotmp[!is.na(kotmp$loc),]

genelocko<-hash(kotmp$gene,kotmp$loc)
crcombineko$Location1<-values(genelocko,keys=crcombineko$left)
crcombineko$Location2<-values(genelocko,keys=crcombineko$right)
crcombineko$contacts.gfp<-
  apply(crcombineko[,c("Location1","Location2")],1,function(x){return(gfp[x[1],x[2]])})
crcombineko$contacts.ko<-
  apply(crcombineko[,c("Location1","Location2")],1,function(x){return(ko[x[1],x[2]])})
crcombineko2<-crcombineko[crcombineko$leftchr!=crcombineko$rightchr,]
crcombineko2$diff<-crcombineko2$contacts.ko-crcombineko2$contacts.gfp

crkoplot<-data.frame(
  contacts=c(crcombineko2$contacts.gfp,crcombineko2$contacts.ko),
  group=c(rep('Ctrl',nrow(crcombineko2)),rep('KO',nrow(crcombineko2))))
  
fig7g<-ggplot(crkoplot, aes(x=group, y=log(contacts))) + 
  geom_boxplot(aes(fill = group), alpha = 0.5,show.legend = F) +
  labs(x='',
       y='log2(normalized contacts between CRGs)',
       title = "proximity of CRGs affected by RUNX1 knockout")
pdf("../../../figures/figure.7.CRGs.coexp/figure7G.RUNX1_KO.pdf",width = 7,height = 4)
fig7g
dev.off()

#-------focus on the gene pairs of CRGs that are closer to each other-----
#---Retrive the 
crcombinetop<-crcombine.eud[crcombine.eud$ave.eud<23,1:2] #23(~1.5μm in the nucleus)
crcombinetop<-merge(crcombineko2,crcombinetop)
#Paired wilcox.test: p-value = 
library(plyr)
nfig7g<-ggplot(crcombinetop, aes(x = diff)) +
  geom_density(color = "black",fill = "red",alpha=0.2) +
  geom_vline(aes(xintercept =mean(crcombinetop$diff)),linetype="dashed") +
  scale_y_continuous(limits = c(0,0.8)) +
  scale_x_continuous(limits = c(-2,2)) +
  theme_bw() +
  theme(panel.grid=element_blank())
pdf("../../../figures/figure.7.CRGs.coexp/figure7G.RUNX1_KO.gb.pdf",width = 7,height = 4)
nfig7g
dev.off()
library(lsr)
cohensD(crcombinetop$diff,mu = 0)
#--------Annotations of MLL partners（结果不显著,不再分析）------------#
mllpartko=mllpart
#Locations of genes
mllpartko$Location1<-rep(1942,nrow(mllpartko))
mllpartko$Location2<-
  apply(mllpartko[,c("chr","loc")],1,function(x){
    any.chr<-x[1]; any.pos<-as.numeric(x[2])
    tmpnum<-which(bin1m$chr==any.chr & any.pos>=bin1m$start & any.pos<=bin1m$end)[1]
    return(tmpnum);rm(tmpnum)}
  )
mllpartko$contacts.ko<-
  apply(mllpartko[,c("Location1","Location2")],1,function(x){
    sum(ko[x[1],x[2]])})
mllpartko$contacts.gfp<-
  apply(mllpartko[,c("Location1","Location2")],1,function(x){
    sum(gfp[x[1],x[2]])})
mllkoplot<-mllpartko[mllpartko$chr!=11,]

save.image("../02_fusion_distances.rda")
