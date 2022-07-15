bulkdfus<-dfusions.eud[,c(4,5,7,8)]
bulkdfus$group<-rep('fusion',nrow(bulkdfus))
bulkctrl<-cbind(dis[ctrlfus.eud$Location1,1:2],dis[ctrlfus.eud$Location2,1:2])
bulkctrl$group<-rep('ctrl',nrow(bulkctrl))
colnames(bulkctrl)<-colnames(bulkdfus)
bulkdfus<-rbind(bulkdfus,bulkctrl)

bulkdfus$leftchr<-as.numeric(bulkdfus$leftchr)
bulkdfus$rightchr<-as.numeric(bulkdfus$rightchr)

for (i in 1:nrow(bulkdfus)){if(bulkdfus[i,1]>bulkdfus[i,3]){
  bulkdfus[i,]<-bulkdfus[i,c(3,4,1,2,5)]}}

bulkdfus$Location1<-round(bulkdfus$leftloc/250000,0)
bulkdfus$Location2<-round(bulkdfus$rightloc/250000,0)
bulkdfus$rightchr[bulkdfus$rightchr==23]<-'X'
bulkdfus$chrs<-paste(bulkdfus$leftchr,bulkdfus$rightchr,sep = "_")


path<-"/work/programs/leukemia_distances/positive_control/gm12878/GM12878_combined_inter/250kb_interchr/inter_contacts/"
setwd(path)
filenames <- dir(path)


chrnames<-c()
for (i in 1:length(filenames)){
  chrs<-unlist(strsplit(gsub('chr','',filenames[i]),'_'))[1:2]
  chrnames<-c(chrnames,paste(chrs,collapse='_'))
}

library(data.table)
datatablelist = list()
for(i in 1:length(filenames)){
datatablelist[[chrnames[i]]] = fread(filenames[i])
chrs<-unlist(strsplit(gsub('chr','',filenames[i]),'_'))[1:2]
datatablelist[[chrnames[i]]]$V1<- datatablelist[[chrnames[i]]]$V1/250000
datatablelist[[chrnames[i]]]$V2<- datatablelist[[chrnames[i]]]$V2/250000
}

bulkdfus$contacts<-apply(bulkdfus,1,function(x){
  x<-as.character(x)
  x[8]<-as.character(x[8])
  setkey(datatablelist[[x[8]]],V1,V2)
  leftloc<-as.numeric(x[6])
  rightloc<-as.numeric(x[7])
  num<-datatablelist[[x[8]]][list(leftloc,rightloc)]$V3
}
)

bulkdfus$rightchr[bulkdfus$rightchr %in% 'X']<-23
#---comparisions between single-cell Hi-C and bulk Hi-C---#
dfuscor<-dfusions.eud[,c(1:2,4,5,7,8,46,47)]
dfuscor$leftchr<-as.numeric(dfuscor$leftchr)
dfuscor$rightchr<-as.numeric(dfuscor$rightchr)
for (i in 1:nrow(dfuscor)){if(dfuscor[i,3]>dfuscor[i,5]){dfuscor[i,]<-dfuscor[i,c(2,1,5,6,3,4,7,8)]}}
dfuscor<-merge(dfuscor,bulkdfus)

#--correlation--#
library(ggplot2)
p1<-ggplot(data = dfuscor, mapping = aes(x=ave.eud, y = contacts)) +
  geom_point(size=1) + 
  stat_smooth(method = 'lm') +
  theme(legend.position="top") +
  scale_y_continuous(limits = c(0,150))+
  labs(x='Distances by single-cell Hi-C',y='Bulk Hi-C contacts')
#      geom_text(aes(label=fusion), size=2)
#p-value = 6.004e-07, cor=-0.6012728.

#--significance--#
p2<-ggplot(bulkdfus,aes(x=group, y=contacts)) + 
  geom_boxplot(aes(fill = group), alpha = 0.5,show.legend = F) +
  scale_y_continuous(limits = c(0,150))+
  labs(x='Group',y='Bulk Hi-C contacts')

library("cowplot")
pdf("../../../figures/figure.bulk_Hi-C.pdf",width = 10,height = 4.5)
plot_grid(p1,p2,labels=c('A','B'),rel_widths = c(1,1),nrow = 1)
dev.off()

library(lsr) #cohensD
#--p-value = 4.664e-12,cohenD=1.3774
