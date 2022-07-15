library(GOplot)
library(ggplot2)
library(enrichplot)
install.packages("ggvenn")
library(ggvenn)

chea=read.table("126_genes_Enrichr_ChEA_2016_table.txt",sep = "\t",header = T,stringsAsFactors = F)
rownames(chea)<-chea$ID
chea$Description<-chea$ID
archs4=read.table("126_genes_Enrichr_ARCHS4_TFs_Coexp_table.txt",sep = "\t",header = T,stringsAsFactors = F)
rownames(archs4)<-archs4$ID
archs4$Description<-archs4$ID
perturb=read.table("126_genes_Enrichr_TF_Perturbations_Followed_by_Expression_table.txt",sep = "\t",header = T,stringsAsFactors = F)
rownames(perturb)<-perturb$ID
perturb$Description<-perturb$ID
#------Venn diagram of three databases------#
vennDat <- list(
  ChEA_TF = chea$ID, 
  ARCHS4_TF_Coexp = archs4$ID, 
  TF_Exp_Perturbations = perturb$ID
)
pdf("venn.pdf",width = 6,height = 6)
ggvenn(vennDat,stroke_size = 0.1)
dev.off()

#------Cluster enrichment-----#
#生成agrigo结果的enrichResult类
runx1=read.table("RUNX1_Enrichr_TFs.txt",sep = "\t",header = T,stringsAsFactors = F)
rownames(runx1)<-runx1$ID
geneID_all <- unlist(apply(as.matrix(runx1$geneID),1,function(x) unlist(strsplit(x,'/'))))
runx1Data <- new("enrichResult",
                 result         = runx1,
                 pvalueCutoff   = 0.05,
                 pAdjustMethod  = "BH",
                 qvalueCutoff   = 0.05,
                 gene           = geneID_all,
                 universe       = 'Unknown',
                 geneSets       = list(),
                 organism       = "Unknown",
                 keytype        = "Unknown",
                 readable       = FALSE
)
pdf("RUNX1_regulation.pdf",width = 7,height = 4)
dotplot(runx1Data,font.size = 9) + ggtitle("dotplot for RUNX1 regulation enrichment")
dev.off()

#------------------TF-------------#
geneID_all <- unlist(apply(as.matrix(chea$geneID),1,function(x) unlist(strsplit(x,'/'))))
cheaData <- new("enrichResult",
                result         = chea,
                pvalueCutoff   = 0.05,
                pAdjustMethod  = "BH",
                qvalueCutoff   = 0.05,
                gene           = geneID_all,
                universe       = 'Unknown',
                geneSets       = list(),
                organism       = "Unknown",
                keytype        = "Unknown",
                readable       = FALSE
)
cp1<-dotplot(cheaData,font.size = 9) + ggtitle("dotplot for RUNX1 regulation enrichment")

geneID_all <- unlist(apply(as.matrix(archs4$geneID),1,function(x) unlist(strsplit(x,'/'))))
archs4Data <- new("enrichResult",
                  result         = archs4,
                  pvalueCutoff   = 0.05,
                  pAdjustMethod  = "BH",
                  qvalueCutoff   = 0.05,
                  gene           = geneID_all,
                  universe       = 'Unknown',
                  geneSets       = list(),
                  organism       = "Unknown",
                  keytype        = "Unknown",
                  readable       = FALSE
)
cp2<-dotplot(archs4Data,font.size = 9) + ggtitle("dotplot for ARCHS4 regulation enrichment")

geneID_all <- unlist(apply(as.matrix(perturb$geneID),1,function(x) unlist(strsplit(x,'/'))))
perturbData <- new("enrichResult",
                   result         = perturb,
                   pvalueCutoff   = 0.05,
                   pAdjustMethod  = "BH",
                   qvalueCutoff   = 0.05,
                   gene           = geneID_all,
                   universe       = 'Unknown',
                   geneSets       = list(),
                   organism       = "Unknown",
                   keytype        = "Unknown",
                   readable       = FALSE
)
cp3<-dotplot(perturbData,font.size = 9) + ggtitle("dotplot for perturb enrichment")

cp<-cowplot::plot_grid(cp1,cp2,cp3,ncol=1, labels=LETTERS[1:3])
pdf("three_enrichments.pdf",width = 6, height = 12)
cp
dev.off()

pdf("RUNX1_network.pdf",width = 6,height = 4)
cnetplot(runx1Data,
         cex_label_category=0.45,
         cex_label_gene=0.35,
         cex_gene=0.45
)
dev.off()
