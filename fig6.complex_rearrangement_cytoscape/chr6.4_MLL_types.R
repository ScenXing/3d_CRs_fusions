
library(RCircos)
crg<-crgenes[crgenes$chr!=crgenes$Chromosome,]
crg<-crg[!is.na(crg$fusion),]
#---limited to 4 genes on chr6 ('CMAHP','BTN3A1','TNXB','DHX16')---#
chr6gene<-c('CMAHP','BTN3A1','TNXB','DHX16')

chr6.link<-circos.link[circos.link$left %in% chr6gene | circos.link$right %in% chr6gene,]
#---CIRCOS Basic chromosome plot---#
BasicPlot <- function(...){ 
  data(UCSC.HG19.Human.CytoBandIdeogram);
  head(UCSC.HG19.Human.CytoBandIdeogram);
  #chr.exclude <- c("chrY");
  chr.exclude <- NULL;
  cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
  tracks.inside <-3;  
  tracks.outside<-0;
  RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
}

OtherGenePlot <- function(...){
  #---Add Other gene name---#
  #-plot more genes-#
  params <- RCircos.Get.Plot.Parameters() #$char.width
  params$char.width <- 100 #default 500
  params$base.per.unit<-50000 #default 30000
  params$text.size <- 0.2 #default 0.4
  RCircos.Reset.Plot.Parameters(params)
  
  circos.other<-unique(crg[,c(4:6,2)])
  colnames(circos.other)<-c('Chromosome','chromStart','chromEnd','Gene')
  
  partg<-unique(crg[,c(10:12,1)])
  partg<-partg[!is.na(partg$fusion),]
  colnames(partg)<-c('Chromosome','chromStart','chromEnd','Gene')
  circos.other<-rbind(circos.other,partg)
  
  circos.other$Chromosome<-paste('chr',circos.other$Chromosome,sep = "")
  side<-"in"
  track.num<-1
  RCircos.Gene.Connector.Plot(circos.other, track.num, side)
  name.col <- 4
  track.num <- 2
  RCircos.Gene.Name.Plot(circos.other, name.col,track.num, side)
  rm(circos.other)
}

fusGenePlot <- function(...){
  #---Add fusion genes---#
  params$text.color<-'red'
  params$char.width <- 220 #default 500
  params$text.size <- 0.2 #default 0.4
  RCircos.Reset.Plot.Parameters(params)
  
  circos.fus<-unique(crg[,c(10:12,1)])
  colnames(circos.fus)<-c('Chromosome','chromStart','chromEnd','Gene')
  circos.fus<-circos.fus[!is.na(circos.fus$Chromosome),]
  circos.fus$Chromosome<-paste('chr',circos.fus$Chromosome,sep = "")
  circos.fus[nrow(circos.fus)+1,]<-c('chr11','118307205','118397539','KMT2A')
  circos.fus$chromStart<-as.numeric(circos.fus$chromStart)
  circos.fus$chromEnd<-as.numeric(circos.fus$chromEnd)
  
  name.col <- 4
  track.num <- 2
  RCircos.Gene.Name.Plot(circos.fus, name.col,track.num, side)
}

fusLinkPlot <- function(...){
  #-Real fusions-#
  params$line.color<-'red'
  RCircos.Reset.Plot.Parameters(params)
  track.num<-3
  circos.realfus=crg[crg$gene %in% chr6gene,c(4:6,10:12)]
  circos.realfus$chr<-paste('chr',circos.realfus$chr,sep = "")
  circos.realfus$Chromosome<-paste('chr',circos.realfus$Chromosome,sep = "")
  colnames(circos.realfus)<-c("ChromA","chromStartA","chromEndA",
                              "ChromB","chromStartB","chromEndB")
  circos.realfus<-circos.realfus[!is.na(circos.realfus$chromStartB),]
  circos.realfus<-circos.realfus[circos.realfus$ChromA!=circos.realfus$ChromB,]
  RCircos.Link.Plot(circos.realfus,track.num,FALSE,lineWidth = 0.1)
}

closeLinkPlot <- function(...){ 
  #---Add lines ---#
  #-genes with Eud<10-#
  params$line.color<-'green'
  RCircos.Reset.Plot.Parameters(params)
  track.num <- 3
  RCircos.Link.Plot(chr6.link[,3:8],track.num,TRUE,lineWidth = 0.1)
}

#---figure.6/adjacent.genes---#
#pdf(file="fig7a.complex_rearrangement_cytoscape/cr.circos.adjacent.genes.pdf", height=14, width=14, compress=TRUE)
BasicPlot()
OtherGenePlot()
fusGenePlot()
closeLinkPlot()
#dev.off()

#---figure.6/fusion.link---#
#pdf(file="fig7a.complex_rearrangement_cytoscape/cr.circos.fusion.link.pdf", height=14, width=14, compress=TRUE)
pdf("tmp.pdf")
BasicPlot()
OtherGenePlot()
fusGenePlot()
fusLinkPlot()
dev.off()
