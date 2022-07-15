library(ggrepel)
library("cowplot")
#save(crgenes,mllpart,file="mechanism.rda")
load("mechanism.rda")
k562=read.table("table.K562.100kb.txt",sep = "\t",header = T,stringsAsFactors = F)
k562$DSB=k562$DSB*100
k562$EXPRESSION=log10(as.numeric(k562$EXPRESSION))
k562cr<-k562[k562$GENE %in% c(crgenes$fusion,crgenes$gene,'KMT2A'),]
k562cr<-k562cr[!duplicated(k562cr$GENE),]
k562mll<-k562[k562$GENE %in% mllpart$gene,]
k562mll<-k562mll[!duplicated(k562mll$GENE),]
k562plot<-plot3cor(k562,k562cr)

cd34=read.table("table.CD34.100kb.txt",sep = "\t",header = T,stringsAsFactors = F)
cd34$DSB=cd34$DSB*100
cd34$EXPRESSION=log10(as.numeric(cd34$EXPRESSION))
cd34cr<-cd34[cd34$GENE %in% c(crgenes$fusion,crgenes$gene,'KMT2A'),]
cd34cr<-cd34cr[!duplicated(cd34cr$GENE),]
cd34mll<-cd34[cd34$GENE %in% mllpart$gene,]
cd34plot<-plot3cor(cd34,cd34cr)

tk6=read.table("table.TK6.100kb.txt",sep = "\t",header = T,stringsAsFactors = F)
tk6$DSB=tk6$DSB*100
tk6$EXPRESSION=log10(as.numeric(tk6$EXPRESSION))
tk6cr<-tk6[tk6$GENE %in% c(crgenes$fusion,crgenes$gene,'KMT2A'),]
tk6cr<-tk6cr[!duplicated(tk6cr$GENE),]
tk6mll<-tk6[tk6$GENE %in% mllpart$gene,]
tk6plot<-plot3cor(tk6,tk6cr)

pdf("dis_dsb_exp.test.pdf",width = 18, height = 5)
plot_grid(k562plot,cd34plot,tk6plot,
          labels = LETTERS[1:3],rel_widths = c(1,1,1),
          ncol =3, nrow = 1)
dev.off()

plot3cor<-function(dat,datcr){
pplot = ggplot()
pplot = pplot + geom_point(data=dat, aes(x=DISTANCE, y=DSB, color=EXPRESSION), shape=20, size=1, alpha = 0.25)
pplot = pplot + scale_colour_gradient2(na.value = "lightblue", low = "lightblue", mid = "lightblue", high = "red", midpoint = 0, limits = c(-2.5, 4.5), breaks = c(-2, 4), labels = c(-2, 4), name = "Expression\nlog10(RPKM)", guide = "colourbar")
pplot = pplot + xlab("Distance from the closest loop anchor (kb)") 
pplot = pplot + ylab("Normalized DSBs per gene")
#pplot = pplot + ggtitle(paste(dat,' cell',sep = ""))
pplot = pplot + theme_classic()
pplot = pplot + scale_x_continuous(breaks = seq(0,100, 25), labels = seq(0,100, 25), limits = c(0,100))
pplot = pplot + ylim(0, 35)
pplot = pplot + theme(legend.title=element_text(size=8), plot.title = element_text(hjust = 0.5))
pplot = pplot + guides(colour = guide_colourbar(order = 1))
#Add CR genes
pplot = pplot + geom_point(data=datcr, aes(x=DISTANCE, y=DSB), shape=21, size=2, alpha = 1,colour="darkgreen")
pplot = pplot + geom_text_repel(data=datcr[datcr$DSB>0.6 & datcr$DISTANCE<90,],aes(x=DISTANCE, y=DSB,label = GENE),size=2)
return(pplot)
}


all3cells<-rbind(k562cr,tk6cr,cd34cr)
#·½²î·ÖÎö
summary(aov(DSB~DISTANCE*EXPRESSION,data = all3cells))
#Df    Sum Sq   Mean Sq F value   Pr(>F)    
#  DISTANCE              1  14025556  14025556   7.707  0.00608 ** 
#  EXPRESSION            1 362806684 362806684 199.372  < 2e-16 ***
#  DISTANCE:EXPRESSION   1  71865884  71865884  39.492 2.41e-09 ***
#  Residuals           180 327554682   1819748                     


#save.image("mechanism.rda")
