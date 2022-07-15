library("biomaRt") #retrive information for ensembl
library(hash) #perl hash
library("reshape2") #melt


#crch37(hg19):get info from ensembl
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
coexptmp<-
  getBM(attributes = c("external_gene_name","ensembl_gene_id","chromosome_name","start_position","end_position","strand"),
        filters = "external_gene_name",
        values = crgenes$gene,
        mart=mart)
colnames(coexptmp)<-c('ensembl','gene','chr','start','end','strand')
coexptmp<-coexptmp[!grepl('_',coexptmp$chr),]
setdiff(crgenes$gene,coexptmp$gene) #unmatched ID

coexpdat=read.table("co_exp/GSE107011_TPM.txt",sep = "\t",header = T,stringsAsFactors = F)
#--Get match ENSEMBL ID--#
conum<-unlist(apply(coexptmp,1,function(x){
  as.numeric(grep(x[1],coexpdat$X))}
  ))
coexpdat<-cbind(coexptmp[names(conum),],coexpdat[conum,])
coexp<-coexpdat[,9:ncol(coexpdat)]
row.names(coexp)=coexpdat$gene


coexp<-coexp[rowMeans(coexp)>1,]
  
library(Hmisc)
coexpStat <- rcorr(t(as.matrix(coexp)))
coexpStat$r[is.na(coexpStat$r)]<-0
coexpStat$P[is.na(coexpStat$P)]<-0

corrplot(coexpStat$r, type="upper", order="hclust", p.mat = coexpStat$P, sig.level = 0.05, insig = "blank")

pdf("co_exp/coexp_of_crgenes.pdf",width = 8.5,height = 7)
pheatmap(coexpStat$r,
         cluster_rows = T,
         cluster_cols = T,
         show_colnames=F,
         show_rownames=T,
         treeheight_col=20,
         treeheight_row = 30
)
dev.off()

#---corelation figure within several genes---#
# Insignificant correlations are leaved blank
#corrplot(res2$r, type="upper", order="hclust", p.mat = res2$P, sig.level = 0.01, insig = "blank")


#--calculate the percentage distribution--#
percentilerank<-function(x){
  rx<-rle(sort(x))
  smaller<-cumsum(c(0, rx$lengths))[seq(length(rx$lengths))]
  larger<-rev(cumsum(c(0, rev(rx$lengths))))[-1]
  rxpr<-smaller/(smaller+larger)
  rxpr[match(x, rx$values)]
}

percendat=read.table("co_exp/GSE107011_TPM.txt",sep = "\t",header = T,stringsAsFactors = F)
percendat<-percendat[,2:ncol(percendat)]
percendat=apply(coexpdat,2,percentilerank)

pbmc<-cbind(coexptmp[names(conum),],percendat[conum,])
row.names(pbmc)<-pbmc$gene
pbmc<-pbmc[,grep('PBMC',colnames(pbmc))]

pdf("co_exp/crgene_pbmc_percentile_nobordercol.pdf",width = 8,height = 7)
pheatmap(pbmc,
         cluster_rows = T,
         cluster_cols = F,
         fontsize_row = 4, 
         fontsize_col = 7,
         angle_col = "0",
         treeheight_row = 20,
         border_color=NA)
dev.off()
#save.image(file="co_exp/coexp.rda")
