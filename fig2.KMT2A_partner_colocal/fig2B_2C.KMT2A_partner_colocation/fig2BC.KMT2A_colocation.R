kmt2aclus<-dfusions.eud[dfusions.eud$left %in% 'KMT2A',]
row.names(kmt2aclus)<-paste(paste('chr',kmt2aclus$rightchr,sep = ""),kmt2aclus$right,sep = ".")
kmt2aclus<-kmt2aclus[order(as.numeric(kmt2aclus$rightchr)),]
kmt2aclus<-kmt2aclus[,c(gm12878,pbmc)]
colnames(kmt2aclus)<-cells[colnames(kmt2aclus)]

kmt2aclus[kmt2aclus<15]<-1 # spatial proximity or not
kmt2aclus[kmt2aclus>15]<-0
kmt2aclus[is.na(kmt2aclus)]<-0
kmt2aclus<-kmt2aclus[rowSums(kmt2aclus)>2,]

library(pheatmap)
require("RColorBrewer")

kmt2acor<-cor(t(kmt2aclus))
pdf("../../figures/figure.3/KMT2A_cluster/3I_spatial_cor.pdf",width = 7,height = 7)
pheatmap(kmt2acor,
             color = brewer.pal(n = 3, name ="YlOrRd"),
             cluster_rows = T,cluster_cols = T,
             fontsize=6,cexCol=2,
             show_colnames=T,
             treeheight_col=20,
             treeheight_row = 30,
             cellwidth = 8, 
             cellheight = 8
         )
dev.off()

pdf("../../figures/figure.3/KMT2A_cluster/3I_co_locations.pdf",width = 8,height = 7)
pheatmap(kmt2aclus,
         color = brewer.pal(n = 7, name ="YlOrRd")[c(1,2,6)],
         cluster_rows = F,cluster_cols = F,
         fontsize=6,cexCol=2,
         angle_col = "270",
         show_colnames=T,
         treeheight_col=20,
         treeheight_row = 30,
         cellwidth = 12, 
         cellheight = 8
)
dev.off()
#write.table(kmt2aclus,file = "co_locations.txt",sep = "\t",quote = F)