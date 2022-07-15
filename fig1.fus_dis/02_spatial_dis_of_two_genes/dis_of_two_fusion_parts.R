library("magrittr")
library("scatterplot3d")
library("reshape2") #melt
library("ggplot2")
library(hash)
library(plyr)

# remove chr in 'chr1', change X to 23
  dfusions$leftchr[dfusions$leftchr=='X']<-23  #change X to 23
dfusions$rightchr[dfusions$rightchr=='X']<-23  #change X to 23
# remove 'chr' in 'chr1'
#allsv$Chromosome1<-gsub('chr','',allsv$Chromosome1);allsv$Chromosome2<-gsub('chr','',allsv$Chromosome2);
#allsv$Chromosome1[allsv$Chromosome1=='X']<-23;allsv$Chromosome2[allsv$Chromosome2=='X']<-23

#----------------------this start----------------------#
#---Location each parters of SVs in dis(dataframe)
 #-left parter
dfusions$Location1<-
  apply(dfusions[,c("leftchr","leftloc")],1,function(x){
    any.chr<-as.numeric(x[1]); any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>dis$start & any.pos<dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )
 #right parter
dfusions$Location2<-
  apply(dfusions[,c("rightchr","rightloc")],1,function(x){
    any.chr<-as.numeric(x[1]); any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>dis$start & any.pos<dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )

##--calculate the short euclidian distances between two parters--##
dttt<-dfusions.eud #back up of dfusions
#dfusions<-dfusions.eud[,1:11]
dfusions.eud<-dfusions
samples=read.table("samples.txt",header = F,sep = "\t",stringsAsFactors = F)[[1]]
gm12878<-c("GSM3271348","GSM3271349","GSM3271351","GSM3271352","GSM3271353","GSM3271356",
                   "GSM3271362","GSM3271364","GSM3271366","GSM3271368","GSM3271370")
pbmc<-samples[17:34]

#---Change Sample Names to Cell Names---#
# 1-14:gm12878,15:32:pbmc
library(readxl)
cell.names <- read_excel("cell_names.xlsx",col_types = c("text", "text"))
cells<-cell.names$cells
names(cells)<-cell.names$samples;rm(cell.names)

#---Calculate the spatial distances between two genes or loci in single cells---#
for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(dfusions[,c("Location1","Location2")],1,function(x){
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
  dfusions.eud<-cbind(dfusions.eud,eu.d)
  rm(eu.d)
}
rm(i)

colnames(dfusions.eud)<-c(colnames(dfusions),samples)
#dfusions.eud$min.eud<-as.numeric(dfusions.eud$min.eud)
#dfusions.eud$ave.eud<-as.numeric(dfusions.eud$ave.eud)

#GM12878+PBMC
dfusions.eud$ave.eud<-apply(dfusions.eud[,c(gm12878,pbmc)],1,function(x){mean(x,na.rm = T)})
dfusions.eud$min.eud<-apply(dfusions.eud[,c(gm12878,pbmc)],1,function(x){ min(x,na.rm = T)})  #generate inf value
#dfusions.eud$min.eud[is.infinite(dfusions.eud$min.eud)]<-'NA'  #remove inf value
dfusions.eud$min.eud<-as.numeric(dfusions.eud$min.eud)
dfusions.eud$ave.eud<-as.numeric(dfusions.eud$ave.eud)

#-GM12878
dfusions.eud$ave.gm12878<-apply(dfusions.eud[,gm12878],1,function(x){mean(x,na.rm = T)})
dfusions.eud$min.gm12878<-apply(dfusions.eud[,gm12878],1,function(x){ min(x,na.rm = T)})  #generate inf value
#dfusions.eud$min.gm12878[is.infinite(dfusions.eud$min.gm12878)]<-'NA'  #remove inf value
dfusions.eud$min.gm12878<-as.numeric(dfusions.eud$min.gm12878)
dfusions.eud$ave.gm12878<-as.numeric(dfusions.eud$ave.gm12878)
#PBMC
dfusions.eud$ave.pbmc<-apply(dfusions.eud[,pbmc],1,function(x){mean(x,na.rm = T)})
dfusions.eud$min.pbmc<-apply(dfusions.eud[,pbmc],1,function(x){ min(x,na.rm = T)})  #generate inf value
#dfusions.eud$min.pbmc[is.infinite(dfusions.eud$min.pbmc)]<-'NA'  #remove inf value
dfusions.eud$min.pbmc<-as.numeric(dfusions.eud$min.pbmc)
dfusions.eud$ave.pbmc<-as.numeric(dfusions.eud$ave.pbmc)

###################################################
#######----Random positions as control---######
###################################################
nleft<-sample(1:nrow(dis),1000,replace=F)
nright<-sample(1:nrow(dis),1000,replace=F)

ctrlsv<-as.data.frame(cbind(nleft,nright))
ctrlsv<-ctrlsv[which(dis$chr[nleft]!=dis$chr[nright]),]
colnames(ctrlsv)<-c("Location1","Location2")
rm(nleft,nright)

for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(ctrlsv,1,function(x){
      ## 7:9 refer to x,y,z coordinates
      mat1<-as.numeric(mat.hic[x[1],7:9]); pat1<-as.numeric(pat.hic[x[1],7:9])   #two alleles of location 1
      mat2<-as.numeric(mat.hic[x[2],7:9]); pat2<-as.numeric(pat.hic[x[2],7:9])   #two alleles of location 2
      tmpd=min(c(
        sqrt((mat1[1]-mat2[1])^2+(mat1[2]-mat2[2])^2+(mat1[3]-mat2[3])^2),
        sqrt((mat1[1]-pat2[1])^2+(mat1[2]-pat2[2])^2+(mat1[3]-pat2[3])^2),
        sqrt((pat1[1]-mat2[1])^2+(pat1[2]-mat2[2])^2+(pat1[3]-mat2[3])^2),
        sqrt((pat1[1]-pat2[1])^2+(pat1[2]-pat2[2])^2+(pat1[3]-pat2[3])^2)
      ))
      return(tmpd);rm(tmpd)
    }
    )
  ctrlsv<-cbind(ctrlsv,eu.d)
}
rm(i,eu.d)

colnames(ctrlsv)<-c(colnames(ctrlsv)[1:2],samples)
ctrlsv$ave.eud<-apply(ctrlsv[,samples],1,function(x){mean(x,na.rm = T)})
ctrlsv$min.eud<-apply(ctrlsv[,samples],1,function(x){ min(x,na.rm = T)})  #generate inf value
ctrlsv$min.eud[is.infinite(ctrlsv$min.eud)]<-'NA'  #remove inf value
ctrlsv$min.eud<-as.numeric(ctrlsv$min.eud)
ctrlsv<-ctrlsv[ctrlsv$ave.eud!='NaN',]

#----ctrlsv.eud:remove rows with much NA----#
#--ctrlsv.eud:33 coloum, location1, location2, 11 GM12878 and 18 PBMC.
na5<-apply(ctrlsv,1,function(x){if(length(which(is.na(x)))<5){return(TRUE)}else{return(FALSE)}}) #romove rows with >=5 NAs
ctrlsv.eud<-ctrlsv[na5,c("Location1","Location2",gm12878,pbmc,'ave.eud','min.eud')]
rm(na5)
#GM12878+PBMC
ctrlsv.eud$ave.eud<-apply(ctrlsv.eud[,c(gm12878,pbmc)],1,function(x){mean(x,na.rm = T)})
ctrlsv.eud$min.eud<-apply(ctrlsv.eud[,c(gm12878,pbmc)],1,function(x){ min(x,na.rm = T)})  #generate inf value
ctrlsv.eud$min.eud[is.infinite(ctrlsv.eud$min.eud)]<-'NA'  #remove inf value
ctrlsv.eud$min.eud<-as.numeric(ctrlsv.eud$min.eud)
#-GM12878
ctrlsv.eud$ave.gm12878<-apply(ctrlsv.eud[,gm12878],1,function(x){mean(x,na.rm = T)})
ctrlsv.eud$min.gm12878<-apply(ctrlsv.eud[,gm12878],1,function(x){ min(x,na.rm = T)})  #generate inf value
ctrlsv.eud$min.gm12878[is.infinite(ctrlsv.eud$min.gm12878)]<-'NA'  #remove inf value
ctrlsv.eud$min.gm12878<-as.numeric(ctrlsv.eud$min.gm12878)
#PBMC
ctrlsv.eud$ave.pbmc<-apply(ctrlsv.eud[,pbmc],1,function(x){mean(x,na.rm = T)})
ctrlsv.eud$min.pbmc<-apply(ctrlsv.eud[,pbmc],1,function(x){ min(x,na.rm = T)})  #generate inf value
ctrlsv.eud$min.pbmc[is.infinite(ctrlsv.eud$min.pbmc)]<-'NA'  #remove inf value
ctrlsv.eud$min.pbmc<-as.numeric(ctrlsv.eud$min.pbmc)

###---Remove duplicated fusions in dfusions, such as KMT2A-AFF1 and AFF1-KMT2A---###
#cosmic=read.table("cosmic_all_fusions.txt",header = T,stringsAsFactors = F,sep = "\t")
#cosmic<-cosmic[,c(1,5,9,10,11)]
#dfusions.eud<-merge(dfusions.eud,cosmic)