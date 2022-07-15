##---Figure 1---Gene Annotations----##
library(biomaRt)
library(regioneR)
gene.symbols <- unique(c(dfusions$left,dfusions$right))
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
fusiongenes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
seqlevelsStyle(fusiongenes) <- "UCSC"

#-head(fusiongenes)

pdf("gene_on_the_chromosome.pdf",9,9)
library(karyoploteR)
kp <- plotKaryotype(genome="hg19")
kpPlotMarkers(kp, data=fusiongenes, labels=fusiongenes$hgnc_symbol, text.orientation = "horizontal",
              r1=0.5, cex=0.5, adjust.label.position = F,label.color="red")
dev.off()

# modified by AI

##---Figure 2-------##
leftgene<-dfusions[,c('left','leftchr')]
colnames(leftgene)<-c('gene','chr')

rightgene<-dfusions[,c('right','rightchr')]
colnames(rightgene)<-c('gene','chr')

genechr<-rbind(leftgene,rightgene)
genechr<-unique(genechr)
genechr<-table(genechr$chr)
genechr<-data.frame(genechr)
colnames(genechr)<-c('chr','num')

chrsize=read.table("chrom.sizes.hg19.txt")
numsize<-merge(genechr,chrsize,by.x = 'chr',by.y = 'V1')
numsize<-numsize[,2:3]
colnames(numsize)<-c('num','size')
numsize$size<-numsize$size/10^6
numsize$num<-as.numeric(numsize$num)
pdf("num_size.pdf",width = 4,height = 3)
  ggplot(numsize, aes(x=size,y=num)) +
  geom_point(size=1,shape=21) + 
  scale_x_continuous(limits=c(0,300)) +
  scale_y_continuous(limits=c(0,11))
dev.off()


#------gene density plot------#
chrden<-read.table("PCG_num_in_chr.txt")
colnames(chrden)<-c('numPCG','chr')
chrden<-chrden[,c('chr','numPCG')]
chrden$chr[chrden$chr=='X']<-23

numden<-merge(genechr,chrden)
numden<-numsize[,2:3]

library(ggplot2)
library(ggpmisc)
pdf("num_genedensity.v2.pdf",width = 4,height = 3)
ggplot(numden, aes(x=numPCG,y=num)) +
  geom_point(size=1,shape=20) + 
  scale_x_continuous(limits=c(0,2500)) +
  scale_y_continuous(limits=c(0,11)) +
  theme_classic(base_size = 16)+
  geom_smooth(method = 'lm', formula = y ~0+x)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')), formula = y ~ 0+x, parse = T)
#geom_text(aes(x=numPCG,y=num+0.5,label=chr))+
dev.off()

fit=lm(numden$num~numden$numPCG+0)
summary(fit)

save.image("Supplementary_figure_1.rda")
