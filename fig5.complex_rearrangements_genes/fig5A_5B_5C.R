
####################################################
#####-Figure 7A--CTRL of complex rearrangement--####
####################################################
library("reshape2") #melt

afgenes<-rbind(af4boxp,af9boxp,af10boxp)
afgenes<-afgenes[,c('ave.eud','min.eud','group')]
afgenes$fusion<-c(rep('KMT2A-AF4',nrow(af4boxp)),rep('KMT2A-AF9',nrow(af9boxp)),rep('KMT2A-AF10',nrow(af10boxp)))
afgenes<-melt(afgenes,id.vars = c("fusion","group"))
colnames(afgenes)<-c('Fusion','Group','Stat','Distances')
afgenes$Stat<-as.character(afgenes$Stat)
afgenes$Stat[afgenes$Stat=='ave.eud']<-'Average'
afgenes$Stat[afgenes$Stat=='min.eud']<-'Minimum'
afgenes$Group[afgenes$Group=='complex_rearragements']<-'CR'


fig7b<-ggplot(afgenes, aes(x=Fusion, y=Distances)) + 
  geom_boxplot(aes(fill = Group), alpha = 0.5,show.legend = T) + 
  facet_grid(Stat~.,scales = "free")+
  labs(x='Complex Rearrangements')
  #theme(legend.position = c(0.9,0.7))
pdf("../../figures/figure.5/figure5B.pdf",width = 8,height = 5)
fig7b
dev.off()
####################################################
#####-Figure 7B-AF4 complex rearrangement in genome--####
####################################################
af4tmp<-af4boxp[af4boxp$chr %in% unique(af4part.eud$chr),c(1:31,34,35)]
af4tmp<-melt(af4tmp,id.vars = c("chr","start","loci","group"))
af4tmp<-arrange(af4tmp,af4tmp$chr,af4tmp$start)
af4tmp<-af4tmp[,3:6]

colnames(af4tmp)<-c('loci','group','sample','value')
af4tmp$loci<-as.character(af4tmp$loci)
af4tmp$loci<-factor(af4tmp$loci,level=unique(af4tmp$loci))
af4tmp$chr<-af4tmp$loci
af4tmp$chr<-as.character(af4tmp$chr)
af4tmp$chr<-as.data.frame(matrix(unlist(strsplit(af4tmp$chr,"\\.")),ncol = 2,byrow = T))[[1]]

fig7c<-
  ggplot(af4tmp, aes(x=loci, y=value)) + 
  geom_boxplot(aes(fill = group), alpha = 0.6,show.legend = T,outlier.shape = NA)+
  theme(axis.text.x = element_text(face = "bold",size = 5, vjust = 0.5, hjust = 0.5, angle = 90))+
  geom_hline(aes(yintercept=median(af4ctrl.eud$ave.eud)))+
  theme(panel.grid.major=element_line(colour=NA))+
  labs(x='genome',y='Distances',title = '                                          KMT2A-AF4 complex rearrangements')+
  theme(legend.position = c(0.9,0.9))
  
#---用chromoplexy或cytoscape图替换
fig7a <- grid::circleGrob()

library("cowplot")
pdf("../../figures/figure.5/figure5.pdf",width = 12,height = 10)
fig7ab<-plot_grid(fig7a,fig7b,labels=c('A','B'),rel_widths = c(1,1.5),nrow = 1)
plot_grid(fig7ab,fig7c,
          labels = c('','C'),
          rel_heights =c(1,1.2),
          ncol = 1, nrow = 2)
dev.off()

#----Figure 5. the ratio of co-localizations of CRs
afratio<-rbind(af4boxp,af9boxp,af10boxp)
afratio$ratio<-apply(afratio[,c(gm12878,pbmc)],1,function(x){x<-x[!is.na(x)];length(x[x<=45])})
afratio<-afratio[,c('ratio','group')]
afratio$fusion<-c(rep('KMT2A-AF4',nrow(af4boxp)),rep('KMT2A-AF9',nrow(af9boxp)),rep('KMT2A-AF10',nrow(af10boxp)))
afratio$group[afratio$group=='complex_rearragements']<-'CR'
colnames(afratio)<-c('Ratio','Group','Fusion')

#--Using barplot similar to qPCR--#
#--EuD<45 as criteria--#
afratioPlot<-data.frame(
  Fusion=c('KMT2A-MLLT10','KMT2A-MLLT10','KMT2A-AF9','KMT2A-AF9','KMT2A-AF4','KMT2A-AF4'),
  Group=c('CR','ctrl','CR','ctrl','CR','ctrl'),
  Ratio=c(4/16,29/229,6/19,34/236,5/21,13/221)
)
afratio$Fusion<-factor(afratio$Fusion,level=c('KMT2A-MLLT10','KMT2A-AF9','KMT2A-AF4'))

af.hist = ggplot(afratioPlot, aes(x = Fusion, y = Ratio,fill = Group))+
  geom_bar(stat ="identity",width = 0.6,position ="dodge")+     
  scale_fill_manual(values = c("orange","darkblue"))+              
  labs(x = "",y = "Percentage of genes", title = "")+    
  theme_bw() +
  theme(panel.grid=element_blank())+
  theme(legend.position = 'top')+
  coord_cartesian(ylim = c(0,0.4))
pdf("../../figures/figure.5/figure5D.pdf",width = 4,height = 4)
af.hist
dev.off()

#matrix(c(15, 46, 76, 686), nrow = 2)
#X-squared = 10.828, df = 1, p-value = 0.0009998
