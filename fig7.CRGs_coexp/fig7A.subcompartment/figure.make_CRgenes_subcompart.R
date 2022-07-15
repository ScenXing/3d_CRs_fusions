library("biomaRt") #retrive information for ensembl
library(hash) #perl hash
library("reshape2") #melt

crgenes=read.table("figure7.complex_rearrangements/MLL_multi_genes.txt",sep = "\t",header = T,stringsAsFactors = F)
colnames(crgenes)<-c('gene','fusion','old')

#crch37(hg19):get info from ensembl
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
crtmp<-
  getBM(attributes = c("external_gene_name", "chromosome_name","start_position","end_position","strand"),
        filters = "external_gene_name",
        values = crgenes$gene,
        mart=mart)
crtmp<-crtmp[!grepl('_',crtmp$chromosome_name),]
setdiff(crgenes[[1]],crtmp$external_gene_name) #unmatched ID

colnames(crtmp)<-c('gene','chr','start','end','strand')
#--transciption start site
crtmp$tss<-crtmp$start
crtmp$tss[crtmp$strand==(-1)]<-crtmp$end[crtmp$strand==(-1)]

crgenes<-merge(crgenes,crtmp)
rm(crtmp)

#right parter
crgenes$Location<-
  apply(crgenes[,c("chr","start")],1,function(x){
    any.chr<-as.character(x[1]); any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>dis$start & any.pos<dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )

#----Sub-Compartment
write.table(crgenes[,c(4:6,2,1,10)],"subcompartment/crgenes.bed",sep = "\t",row.names = F,col.names = F,quote = F)

#---Get MLL partners genes---#
load("../02_fusion_distances.rda")
mlltmp<-
  getBM(attributes = c( "chromosome_name","start_position","end_position","strand","external_gene_name"),
        filters = "external_gene_name",
        values = mllpart$gene,
        mart=mart)
mlltmp<-mlltmp[!grepl('_',mlltmp$chromosome_name),]
setdiff(mllpart[[1]],mlltmp$external_gene_name) #unmatched ID
colnames(mlltmp)<-c('chr','start','end','strand','gene')
mlltmp<-mlltmp[order(mlltmp$chr),]
write.table(mlltmp,"subcompartment/mllgenes.bed",sep = "\t",row.names = F,col.names = F,quote = F)

#----Enrichment of CR-genes and MLL partners  in sub-Compartments
crgenesCompart=read.table("subcompartment/crgenes_in_compartment.bed",sep = "\t",stringsAsFactors = F)
crgenesSimu=read.table("subcompartment/crgenes_simu_in_compartment.bed",sep = "\t",stringsAsFactors = F)
mllgenesCompart=read.table("subcompartment/mllgenes_in_compartment.bed",sep = "\t",stringsAsFactors = F)

compartStat=as.data.frame(cbind(c(table(crgenesCompart$V10),0),table(mllgenesCompart$V9),table(crgenesSimu$V10)))
row.names(compartStat)[6]<-'B4'
compartStat$type=row.names(compartStat)
colnames(compartStat)<-c('crgenes','KMT2A.parters','random','type')
compartStat$crgenes<-compartStat$crgenes/148
compartStat$KMT2A.parters<-compartStat$KMT2A.parters/68
compartStat$random<-compartStat$random/137
library(reshape2)
compartStat<-melt(compartStat)
colnames(compartStat)<-c('type','group','value')

ggplot(compartStat, aes(x = group, y = value,fill = type))+
  geom_bar(stat ="identity",width = 0.6,position ="stack")+
  labs(x = "",y = "the Percentage", title = "")+    
  theme_bw() +
  theme(panel.grid=element_blank())+
  coord_cartesian(ylim = c(0,1.1))+
  scale_y_continuous(breaks=seq(0, 1, 0.2))
ggsave("compartmentPattern.pdf",width = 5.5, height = 4)

#---DSBs per genes---#
dsbGene=read.table("crgenes_dsb/GSM3444989_K562_ETO_gene_based.bed",sep = "\t",stringsAsFactors = F)
dsbGene$V5=dsbGene$V7
dsbGene$V7=dsbGene$V5/((dsbGene$V3-dsbGene$V2)*4.876289/10^3)  #RPKM, sum(dsbGene$V5)/10^6=4.876289
colnames(dsbGene)<-c('chr','start','end','gene','dsbcount','strand','dsbfpkm')
dsbGene<-dsbGene[!duplicated(dsbGene$gene),]

##----KMT2A partners       ----#
dsbMll<-dsbGene[dsbGene$gene %in% mllpart$gene,]

#----Complex rearrangements----#
#dsbCr<-dsbGene[dsbGene$gene %in% c(crgenes$fusion,crgenes$gene,'KMT2A'),]
dsbCr<-dsbGene[dsbGene$gene %in% c(crgenes$gene),]



#--Fig8B-DSBs per gene---#
dsbPlot<-data.frame(DSBperGene=c(dsbCr$dsbfpkm,dsbMll$dsbfpkm,dsbGene$dsbfpkm),
                    Group=c(rep('CR genes',nrow(dsbCr)),
                            rep('MLL partners',nrow(dsbMll)),
                            rep('All genes',nrow(dsbGene))))
library(ggpubr)
fig8dsb<-ggboxplot(dsbPlot, "Group", "DSBperGene", fill = "Group",
                 palette = "npg",
                 outlier.shape = NA,
                 add.params = list(fill = "white"))
fig8dsb<-ggpar(fig8dsb,yscale = "none",legend = "none",ylab = "DSBs per gene",ylim=c(0,5))
fig8dsb<-fig8dsb+stat_compare_means(comparisons=list(c('All genes','MLL partners'),c('CR genes','MLL partners'),c('All genes','CR genes')),
                                label = "p.format",label.y = c(580,580,680))




