#install.packages('RCircos')
library(RCircos)
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
circos.gene<-
  getBM(attributes = c("chromosome_name", "start_position","end_position","external_gene_name"),
        filters = "external_gene_name",
        values = c(crpart,crgenes$gene,'KMT2A'),
        mart=mart)
circos.gene<-circos.gene[!grepl('_',circos.gene$chromosome_name),]
setdiff(c(crpart,crgenes$gene),circos.gene$external_gene_name) #unmatched ID
colnames(circos.gene)<-c('Chromosome','chromStart','chromEnd','Gene')
circos.gene$Chromosome<-paste('chr',circos.gene$Chromosome,sep = "")

#---Make interactions between genes---#
circos.link<-crsigcombine.eud[,c("left","right")]
circos.link<-merge(circos.link,circos.gene,by.x = 'left', by.y = 'Gene')
circos.link<-merge(circos.link,circos.gene,by.x = 'right',by.y = 'Gene')
colnames(circos.link)<-c("right","left",
                         "ChromA","chromStartA","chromEndA",
                         "ChromB","chromStartB","chromEndB")


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
  
  circos.other<-unique(crgenes[,c(4:6,2)])
  colnames(circos.other)<-c('Chromosome','chromStart','chromEnd','Gene')
  circos.other$Chromosome<-paste('chr',circos.other$Chromosome,sep = "")
  side<-"in"
  track.num<-1
  RCircos.Gene.Connector.Plot(circos.other, track.num, side)
  name.col <- 4
  track.num <- 2
  RCircos.Gene.Name.Plot(circos.other, name.col,track.num, side)
}

fusGenePlot <- function(...){
  #---Add fusion genes---#
  params$text.color<-'red'
  params$char.width <- 220 #default 500
  params$text.size <- 0.4 #default 0.4
  RCircos.Reset.Plot.Parameters(params)
  
  circos.fus<-unique(crgenes[,c(10:12,1)])
  colnames(circos.fus)<-c('Chromosome','chromStart','chromEnd','Gene')
  circos.fus<-circos.fus[!is.na(circos.fus$Chromosome),]
  circos.fus$Chromosome<-paste('chr',circos.fus$Chromosome,sep = "")
  circos.fus[nrow(circos.fus)+1,]<-c('chr11','118307205','118397539','KMT2A')
  circos.fus$chromStart<-as.numeric(circos.fus$chromStart)
  circos.fus$chromEnd<-as.numeric(circos.fus$chromEnd)
  
  name.col <- 4
  track.num <- 3
  RCircos.Gene.Name.Plot(circos.fus, name.col,track.num, side)
}
 
fusLinkPlot <- function(...){
  #-Real fusions-#
  params$line.color<-'green'
  RCircos.Reset.Plot.Parameters(params)
  track.num<-4
  circos.realfus=crgenes[,c(4:6,10:12)]
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
  track.num <- 4
  RCircos.Link.Plot(circos.link[,3:8],track.num,TRUE,lineWidth = 0.1)
}

#---figure.6/adjacent.genes---#
pdf(file="fig7a.complex_rearrangement_cytoscape/cr.circos.adjacent.genes.pdf", height=14, width=14, compress=TRUE)
BasicPlot()
OtherGenePlot()
fusGenePlot()
closeLinkPlot()
dev.off()

#---figure.6/fusion.link---#
pdf(file="fig7a.complex_rearrangement_cytoscape/cr.circos.fusion.link.pdf", height=14, width=14, compress=TRUE)
BasicPlot()
OtherGenePlot()
fusGenePlot()
fusLinkPlot()
dev.off()










