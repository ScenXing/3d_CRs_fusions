#source("http://bioconductor.org/biocLite.R");biocLite("biomaRt")
library("biomaRt") #retrive information for ensembl
library(hash) #perl hash

#-----Leukemia fusion genes-------#
fusions<-read.table("Cosmic_Fusions_table.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(fusions)<-c('left','right','tumor')
#---fusions in Leukaemia---#
fusions<-fusions[grep('leukaemia',fusions$tumor,value =F),]

#crch37(hg19):get info from ensembl
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
geneinfo<-
  getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "strand"),
        filters = "external_gene_name",
        values = unique(c(fusions$left,fusions$right)),
        mart=mart)
geneinfo<-geneinfo[!grepl('_',geneinfo$chromosome_name),]
setdiff(unique(c(fusions$left,fusions$right)),geneinfo$external_gene_name) #unmatched ID

chr2cen<-hash(keys=centromeres[[1]],values=centromeres[[2]])
gene2chr<-hash(keys=geneinfo[[1]],values=geneinfo[[2]])
gene2loc<-hash(keys=geneinfo[[1]],values=geneinfo[[3]])
gene2strand<-hash(keys=geneinfo[[1]],values=geneinfo[[4]])

fusions$leftchr<-values(gene2chr,keys=fusions$left)
fusions$leftloc<-values(gene2loc,keys=fusions$left)
fusions$leftstrand<-values(gene2strand,keys=fusions$left)
fusions$leftcentromere<-values(chr2cen,keys=fusions$leftchr)

fusions$rightchr<-values(gene2chr,keys=fusions$right)
fusions$rightloc<-values(gene2loc,keys=fusions$right)
fusions$rightstrand<-values(gene2strand,keys=fusions$right)
fusions$rightcentromere<-values(chr2cen,keys=fusions$rightchr)

#Fusions on different chromsomes
dfusions<-fusions[fusions$leftchr!=fusions$rightchr,]
sfusions<-fusions[fusions$leftchr==fusions$rightchr,]

