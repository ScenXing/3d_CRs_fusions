#source("http://bioconductor.org/biocLite.R");biocLite("biomaRt")
library("biomaRt") #retrive information for ensembl
library(hash) #perl hash

#########################################################
####--------Add fusion in solid tumors to fig1B------####
#########################################################
#----Get fusions in solid tumors---#
fusions<-read.table("Cosmic_Fusions_table.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(fusions)<-c('left','right','tumor')
#---fusions not in Leukaemia---#
ofusions<-fusions[!grepl('leukaemia',fusions$tumor),]

#---Fusions reported in new COSMIC release---#
cosfusions<-read.table("cosmic_all_fusions_annotated.txt",header = T,sep = "\t")
ofusions<-merge(ofusions,cosfusions)  

#-138/214 Fusions are on different chromsomes
ofusions<-ofusions[ofusions$leftchr!=ofusions$rightchr,]

write.csv(ofusions,file="solid_tumors_fusions_on_different_chromosomes.csv",quote = F,row.names = F)

#---Calculate the distances between two fusion partners---#
#-left parter
ofusions$Location1<-
  apply(ofusions[,c("leftchr","leftloc")],1,function(x){
    any.chr<-x[1]; 
    any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>dis$start & any.pos<dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )
#right parter
ofusions$Location2<-
  apply(ofusions[,c("rightchr","rightloc")],1,function(x){
    any.chr<-x[1]; any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>dis$start & any.pos<dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )

##--calculate the short euclidian distances between two parters--##
ofusions.eud<-ofusions
for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(ofusions[,c("Location1","Location2")],1,function(x){
      ## 7:9 refer to x,y,z coordinates
      mat1<-as.numeric(mat.hic[x[1],7:9]); pat1<-as.numeric(pat.hic[x[1],7:9])   #two alleles of location 1
      mat2<-as.numeric(mat.hic[x[2],7:9]); pat2<-as.numeric(pat.hic[x[2],7:9])   #two alleles of location 2
      tmpd=min(c(
        sqrt((mat1[1]-mat2[1])^2+(mat1[2]-mat2[2])^2+(mat1[3]-mat2[3])^2),
        sqrt((mat1[1]-pat2[1])^2+(mat1[2]-pat2[2])^2+(mat1[3]-pat2[3])^2),
        sqrt((pat1[1]-mat2[1])^2+(pat1[2]-mat2[2])^2+(pat1[3]-mat2[3])^2),
        sqrt((pat1[1]-pat2[1])^2+(pat1[2]-pat2[2])^2+(pat1[3]-pat2[3])^2)),
        na.rm = T)
      tmpd[is.infinite(tmpd)]<-NA
      return(tmpd);rm(tmpd)
    }
    )
  ofusions.eud<-cbind(ofusions.eud,eu.d)
  rm(eu.d)
}
rm(i)

colnames(ofusions.eud)<-c(colnames(ofusions),samples)
#ofusions.eud$min.eud<-as.numeric(ofusions.eud$min.eud)
#ofusions.eud$ave.eud<-as.numeric(ofusions.eud$ave.eud)

#GM12878+PBMC
ofusions.eud$ave.eud<-apply(ofusions.eud[,c(gm12878,pbmc)],1,function(x){mean(x,na.rm = T)})
ofusions.eud$min.eud<-apply(ofusions.eud[,c(gm12878,pbmc)],1,function(x){ min(x,na.rm = T)})  #generate inf value
#ofusions.eud$min.eud[is.infinite(ofusions.eud$min.eud)]<-'NA'  #remove inf value
ofusions.eud$min.eud<-as.numeric(ofusions.eud$min.eud)
ofusions.eud$ave.eud<-as.numeric(ofusions.eud$ave.eud)

#-GM12878
ofusions.eud$ave.gm12878<-apply(ofusions.eud[,gm12878],1,function(x){mean(x,na.rm = T)})
ofusions.eud$min.gm12878<-apply(ofusions.eud[,gm12878],1,function(x){ min(x,na.rm = T)})  #generate inf value
#ofusions.eud$min.gm12878[is.infinite(ofusions.eud$min.gm12878)]<-'NA'  #remove inf value
ofusions.eud$min.gm12878<-as.numeric(ofusions.eud$min.gm12878)
ofusions.eud$ave.gm12878<-as.numeric(ofusions.eud$ave.gm12878)
#PBMC
ofusions.eud$ave.pbmc<-apply(ofusions.eud[,pbmc],1,function(x){mean(x,na.rm = T)})
ofusions.eud$min.pbmc<-apply(ofusions.eud[,pbmc],1,function(x){ min(x,na.rm = T)})  #generate inf value
#ofusions.eud$min.pbmc[is.infinite(ofusions.eud$min.pbmc)]<-'NA'  #remove inf value
ofusions.eud$min.pbmc<-as.numeric(ofusions.eud$min.pbmc)
ofusions.eud$ave.pbmc<-as.numeric(ofusions.eud$ave.pbmc)


#--Make dataframe of ofusions for fig 1B
othernum<-nrow(ofusions.eud) #138
ofusionStat<-data.frame(
  Group=rep('Other fusion',138*4),
  Distances=c(ofusions.eud$ave.gm12878,ofusions.eud$min.gm12878,ofusions.eud$ave.pbmc,ofusions.eud$min.pbmc),
  Stat=rep(c('Average','Minimum','Average','Minimum'),c(138,138,138,138)),
  Cell.type=rep(c('GM12878','PBMC'),c(138*2,138*2))
)

ofusionStat<-rbind(ofusionStat,sigfusion)
ofusionStat$Group<-factor(ofusionStat$Group,levels = c('ctrl','Other fusion','fusion'))

#---Plot fig 1B(New)
#pdf("../../figures/figure.1B_v2.pdf",height =5,width =7)
ggplot(ofusionStat, aes(x=Group, y=Distances)) + 
  geom_boxplot(aes(fill = Group), alpha = 0.5,show.legend = F) + 
  facet_grid(Stat ~ Cell.type,scales = "free")
#dev.off()