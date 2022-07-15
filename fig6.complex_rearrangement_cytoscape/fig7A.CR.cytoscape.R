library("biomaRt") #retrive information for ensembl
library(hash) #perl hash
library("reshape2") #melt
####################################################
#---Make complex rearrangements associated genes---#
####################################################
crgenes=read.table("MLL_multi_genes.txt",sep = "\t",header = T,stringsAsFactors = F)
colnames(crgenes)<-c('gene','fusion','old')

#crch37(hg19):get info from ensembl
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
crtmp<-
  getBM(attributes = c("external_gene_name", "chromosome_name","start_position","end_position","strand"),
        filters = "external_gene_name",
        values = crgenes$gene,
        mart=mart)
crtmp<-crtmp[!grepl('_',crtmp$chromosome_name),]
setdiff(crgenes[[1]],crtmp$external_gene_name) #unmatched ID

colnames(crtmp)<-c('gene','chr','start','end','strand')
#--transciption start site
crtmp$tss<-crtmp$start
crtmp$tss[crtmp$strand==(-1)]<-crtmp$end[crtmp$strand==(-1)]

crgenes<-merge(crgenes,crtmp)
rm(crtmp)

#right parter
crgenes$Location<-
  apply(crgenes[,c("chr","start")],1,function(x){
    any.chr<-as.character(x[1]); any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>dis$start & any.pos<dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )

#-mapping gene coresponding location and chr using HASH---#
#Location
crtmp<-crgenes[,c("gene","chr","Location")];colnames(crtmp)<-c('gene','chr','location')
#crpart<-c('MLLT1','ELL','MLLT10','AFF1','MLLT3','SEPT6','MLLT11','AFF3',
#  'MLLT4','TNRC18','ABI1','MYO1F','PICALM','EPS15','SEPT9','VAV1','MLLT6')
#oldname<-c('ENL','ELL','AF10','AF4','AF9','SEPT6','AF1Q','LAF4','AF6',
#           'TNRC18','ABI1','MYO1F','PICALM','EPS15','SEPT9','VAV1','AF17')
crpart<-data.frame(newname=crpart,oldname=oldname)
colnames(crpart)<-c('newname','fusion')

crgenes<-merge(crgenes,circos.gene,by.x = 'fusion',by.y = 'Gene',all.x = T)
crgenes$Chromosome<-unlist(gsub('chr','',crgenes$Chromosome))

crtmp2<-mllpart[mllpart$gene %in% crpart,c("gene","chr","Location2")];colnames(crtmp2)<-c('gene','chr','location')
crtmp<-unique(rbind(crtmp,crtmp2));rm(crtmp2)

crtmp[nrow(crtmp)+1,]<-c('KMT2A',11,104479)

#crcombine<-as.data.frame(t(combn(crtmp$gene,2)))   #All combinations
#colnames(crcombine)<-c('left','right')
#location
geneloc<-hash(crtmp$gene,crtmp$location)
crcombine$Location1<-values(geneloc,keys=crcombine$left)
crcombine$Location2<-values(geneloc,keys=crcombine$right)
 #chromosome
genechr<-hash(crtmp$gene,crtmp$chr)
crcombine$leftchr<-values(genechr,keys=crcombine$left)
crcombine$rightchr<-values(genechr,keys=crcombine$right)

crcombine<-crcombine[crcombine$leftchr!=crcombine$rightchr,]
rm(crtmp2,crtmp)

crcombine.eud<-crcombine
for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  eu.d<-
    apply(crcombine[,c("Location1","Location2")],1,function(x){
      ## 7:9 refer to x,y,z coordinates
      x<-as.numeric(as.character(x))
      dmat<-rbind(mat.hic[x[1],7:9],pat.hic[x[1],7:9],
                  mat.hic[x[2],7:9],pat.hic[x[2],7:9])
      tmpd=min(dist(dmat)[2:5],na.rm = T)
      tmpd[is.infinite(tmpd)]<-NA
      return(tmpd);rm(tmpd)
    }
    )
  crcombine.eud<-cbind(crcombine.eud,eu.d)
  rm(eu.d,mat.hic,pat.hic)
}
rm(i)

colnames(crcombine.eud)<-c(colnames(crcombine),samples)

#GM12878+PBMC
crcombine.eud$ave.eud<-apply(crcombine.eud[,c(gm12878,pbmc)],1,function(x){mean(x,na.rm = T)})
crcombine.eud$min.eud<-apply(crcombine.eud[,c(gm12878,pbmc)],1,function(x){ min(x,na.rm = T)})  #generate inf value
crcombine.eud$min.eud<-as.numeric(crcombine.eud$min.eud)
crcombine.eud$ave.eud<-as.numeric(crcombine.eud$ave.eud)

#-GM12878
crcombine.eud$ave.gm12878<-apply(crcombine.eud[,gm12878],1,function(x){mean(x,na.rm = T)})
crcombine.eud$min.gm12878<-apply(crcombine.eud[,gm12878],1,function(x){ min(x,na.rm = T)})  #generate inf value
crcombine.eud$min.gm12878<-as.numeric(crcombine.eud$min.gm12878)
crcombine.eud$ave.gm12878<-as.numeric(crcombine.eud$ave.gm12878)
#PBMC
crcombine.eud$ave.pbmc<-apply(crcombine.eud[,pbmc],1,function(x){mean(x,na.rm = T)})
crcombine.eud$min.pbmc<-apply(crcombine.eud[,pbmc],1,function(x){ min(x,na.rm = T)})  #generate inf value
crcombine.eud$min.pbmc<-as.numeric(crcombine.eud$min.pbmc)
crcombine.eud$ave.pbmc<-as.numeric(crcombine.eud$ave.pbmc)

#alleud<-as.data.frame(cbind(eud=as.numeric(c(mllpart.eud$ave.eud,crcombine.eud$ave.eud)),
#     type=c(rep('fusion',nrow(mllpart.eud)),rep('combine',nrow(crcombine.eud)))))
#alleud$eud<-as.numeric(as.character(alleud$eud))

#mu <- ddply(alleud, "type", summarise, grp.mean=median(eud,na.rm = T)) # Add mean lines
#pdf("fig3a.cytoscape/sup.fig3.KMT2A.combine.pdf",width = 8,height = 5)
#ggplot(alleud, aes(x=eud, color=type)) +
#  geom_density()+
#  geom_vline(data=mu, aes(xintercept=grp.mean, color=type),linetype="dashed")
#dev.off()
#Result:
#-Eud distances are similar between fusion and crcombine.
#-There are two hills in fusions, while one hill in all combine.

#-Result:cutoff<15, get sigcombine.genes.v2.tsv
crcombine.eud$num<-apply(crcombine.eud[,c(gm12878,pbmc)],1,function(x){x<-x[!is.na(x)];length(x[x<=10])})
crsigcombine.eud<-crcombine.eud[crcombine.eud$num>2,] # at least in 2 samples
write.table(crsigcombine.eud[,c("left","right")],"fig7a.complex_rearrangement_cytoscape/fig7a.cr.sigcombine.genes.tsv",quote = F,row.names = F,sep = "\t")

#-Followed by cytoscope toolkit--#




