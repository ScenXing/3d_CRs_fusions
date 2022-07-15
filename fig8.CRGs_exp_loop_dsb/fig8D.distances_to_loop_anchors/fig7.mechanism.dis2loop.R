#---Gene Distance to Loop Anchors---#
#geneAnchor=read.table("mechanism/gene_relative_to_anchor.txt",header = F,stringsAsFactors = F,sep = "\t")
geneAnchor=read.table("mechanism/gene_relative_to_anchor.txt",sep = "\t",header = F,stringsAsFactors = F)
colnames(geneAnchor)<-c('Gene','Distance','Start','End','CTCF','whichLoop')
geneAnchor$CTCF=read.table("mechanism/ctcf_loc.txt",sep = "\t",header = F,stringsAsFactors = F)
colnames(geneAnchor)<-c('Gene','Distance','Start','End','CTCF','whichLoop')
#geneAnchor$loopAnchor<-geneAnchor$D1+geneAnchor$coverage/2
#-Left Anchors
leftAnchor<- geneAnchor[geneAnchor$whichLoop==16,]
leftAnchor$left <- leftAnchor$Start - leftAnchor$CTCF
leftAnchor$right <-  leftAnchor$End - leftAnchor$CTCF
#-Right Anchors
rightAnchor<-geneAnchor[geneAnchor$whichLoop==0,]
rightAnchor$left <- rightAnchor$CTCF - rightAnchor$End
rightAnchor$right <-  rightAnchor$CTCF - rightAnchor$Start

geneAnchor<-rbind(leftAnchor,rightAnchor)
geneAnchor$left <- trunc(geneAnchor$left/1000)
geneAnchor$right<- trunc(geneAnchor$right/1000)
geneAnchor$Group<-rep('All genes',nrow(geneAnchor))
geneAnchor<-geneAnchor[abs(geneAnchor$left)<1000 & abs(geneAnchor$right)<1000,]

allRange<-unlist(apply(geneAnchor[,c("left","right")],1,function(x){seq(x[1],x[2],by=1)}))
#allRange<-allRange[allRange > (-201) &allRange < 201]

##----KMT2A partners       ----#
mllAnchor<-geneAnchor[geneAnchor$Gene %in% mllpart$gene,]
mllAnchor$Group<-rep('KMT2A partners',nrow(mllAnchor))
mllAnchor<-mllAnchor[abs(mllAnchor$left)<1000 & abs(mllAnchor$right)<1000,]
mllRange<-unlist(apply(mllAnchor[,c("left","right")],1,function(x){seq(x[1],x[2],by=1)}))

#----Complex rearrangements----#
crAnchor<-geneAnchor[geneAnchor$Gene %in% c(crgenes$fusion,crgenes$gene,'KMT2A'),]
crAnchor$Group<-rep('CR genes',nrow(crAnchor))
crAnchor<-crAnchor[abs(crAnchor$left)<1000 & abs(crAnchor$right)<1000,]

crRange<-unlist(apply(crAnchor[,c("left","right")],1,function(x){seq(x[1],x[2],by=1)}))
#crRange<-crRange[crRange > (-201) &crRange < 201]

#---PLot gene distributions---#
anchorPlot<-data.frame(Cov=c(allRange,crRange,mllRange),Group=c(rep('All genes',length(allRange)),
                                                                rep('CR genes',length(crRange)),
                                                                rep('MLL partners',length(mllRange))))
#anchorPlot<-data.frame(Cov=c(allRange,crRange),Group=c(rep('All genes',length(allRange)),
#                                                      rep('CR genes',length(crRange))
#                                                                ))

#colnames(anchorPlot)<-c('Gene','Distance','Group')
p<-ggplot(anchorPlot, aes(x=Cov, fill=Group)) +
  geom_density(alpha=0.4)

#--Fig8A-Distance to Loop Anchor---#
disPlot<-data.frame(DistanceToLoopAnchor=c(crAnchor$Distance,mllAnchor$Distance,geneAnchor$Distance),
                    Group=c(rep('CR genes',nrow(crAnchor)),
                            rep('MLL partners',nrow(mllAnchor)),
                            rep('All genes',nrow(geneAnchor))))
library(ggpubr)
fig8a<-ggboxplot(disPlot, "Group", "DistanceToLoopAnchor", fill = "Group",
                palette = "npg",
                outlier.shape = NA,
                add.params = list(fill = "white"))
fig8a<-ggpar(fig8a,yscale = "none",legend = "none",ylab = "Distances to Loop Anchors",ylim=c(0,800))
fig8a<-fig8a+stat_compare_means(comparisons=list(c('All genes','MLL partners'),c('CR genes','MLL partners'),c('All genes','CR genes')),
                                label = "p.format",label.y = c(580,580,680))
#ggsave("mechanism/8A.distances_to_loop_anchors.pdf")



