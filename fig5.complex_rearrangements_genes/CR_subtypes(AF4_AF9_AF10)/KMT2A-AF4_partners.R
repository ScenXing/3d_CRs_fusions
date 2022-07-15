#source("http://bioconductor.org/biocLite.R");biocLite("biomaRt")
library("biomaRt") #retrive information for ensembl
library(hash) #perl hash
library("reshape2") #melt

#load("fusions_in_leukemia.rda")

#-----Leukemia fusion genes-------#
af4gene<-read.table("complex_rearrangements/KMT2A-AF4_genes.txt",header = F,sep = "\t",stringsAsFactors = F)

#crch37(hg19):get info from ensembl
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
af4tmp<-
  getBM(attributes = c("external_gene_name", "chromosome_name", "start_position"),
        filters = "external_gene_name",
        values = af4gene,
        mart=mart)
af4tmp<-af4tmp[!grepl('_',af4tmp$chromosome_name),]
setdiff(af4gene[[1]],af4tmp$external_gene_name) #unmatched ID
af4gene<-af4tmp
rm(af4tmp)

af4loci=read.table("complex_rearrangements/KMT2A-AF4_loci.txt",sep = "\t",skip = 1)
colnames(af4gene)<-c('gene','chr','start')
colnames(af4loci)<-c('gene','chr','start')
af4part=rbind(af4gene,af4loci)

#-Becase of AF4 is on chromosome4, we remove genes on chr4(AF4) and chr11(KMT2A).
#-We got 21 loci.
af4part<-af4part[af4part$chr!=4,]
af4part<-af4part[af4part$chr!=11,]

#Location each parters of KMT2A-AF4 partners.
#Add = in > and <
af4part$Location<-
  apply(af4part[,c("chr","start")],1,function(x){
    any.chr<-as.numeric(x[1]); any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>=dis$start & any.pos<=dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )

af4part$Location1<-rep(104479,nrow(af4part)) #Add KMT2A
af4part$Location2<-rep(38918,nrow(af4part))  #Add AFF1(AF4)

#-----sum of sides of triangle----#
#sum(dist(x, method = "euclidean",upper = FALSE))
MMM,MMF
MFM,MFF
FMM,FMF
FFM,FFF


af4part.eud<-af4part
for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(af4part[,c("Location","Location1","Location2")],1,function(x){
      ## 7:9 refer to x,y,z coordinates
      mat1<-as.numeric(mat.hic[x[1],7:9]); pat1<-as.numeric(pat.hic[x[1],7:9])   #two alleles of Location
      mat2<-as.numeric(mat.hic[x[2],7:9]); pat2<-as.numeric(pat.hic[x[2],7:9])   #two alleles of location1
      mat3<-as.numeric(mat.hic[x[3],7:9]); pat3<-as.numeric(pat.hic[x[3],7:9])   #two alleles of location2
      
      tmpd=min(c(
        sum(dist(rbind(mat1,mat2,mat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(mat1,mat2,pat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(mat1,pat2,mat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(mat1,pat2,pat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(pat1,mat2,mat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(pat1,mat2,pat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(pat1,pat2,mat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(pat1,pat2,pat3), method = "euclidean",upper = FALSE))),
        na.rm = T)
      tmpd[is.infinite(tmpd)]<-NA
      return(tmpd);rm(tmpd)
    }
    )
  af4part.eud<-cbind(af4part.eud,eu.d)
  rm(eu.d)
}

rm(i)
colnames(af4part.eud)<-c(colnames(af4part),samples)
#GM12878+PBMC
af4part.eud$ave.eud<-apply(af4part.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{mean(x,na.rm = T)}})
af4part.eud$min.eud<-apply(af4part.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{ min(x,na.rm = T)}})
af4part.eud$min.eud<-as.numeric(af4part.eud$min.eud)
af4part.eud$ave.eud<-as.numeric(af4part.eud$ave.eud)

#-GM12878
af4part.eud$ave.gm12878<-apply(af4part.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{mean(x,na.rm = T)}})
af4part.eud$min.gm12878<-apply(af4part.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{ min(x,na.rm = T)}})
af4part.eud$min.gm12878<-as.numeric(af4part.eud$min.gm12878)
af4part.eud$ave.gm12878<-as.numeric(af4part.eud$ave.gm12878)
#PBMC
af4part.eud$ave.pbmc<-apply(af4part.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{mean(x,na.rm = T)}})
af4part.eud$min.pbmc<-apply(af4part.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{ min(x,na.rm = T)}})
af4part.eud$min.pbmc<-as.numeric(af4part.eud$min.pbmc)
af4part.eud$ave.pbmc<-as.numeric(af4part.eud$ave.pbmc)

#----------Controls of KMT2A-AFF1(AF4)-------------#
#---Make controls of AF4 fusions---#
cnum<-sort(sample(nrow(dis),300,replace = F))
af4ctrl<-dis[cnum,]
af4ctrl$Location<-cnum
rm(cnum)
af4ctrl$Location1<-rep(104479,nrow(af4ctrl)) #Add KMT2A
af4ctrl$Location2<-rep(38918,nrow(af4ctrl))  #Add AFF1(AF4)
af4ctrl<-af4ctrl[!af4ctrl$chr %in% c(4,11,'X'),] #remove chr4,chr11 and chrX

#-Distances
af4ctrl.eud<-af4ctrl
for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(af4ctrl[,c("Location","Location1","Location2")],1,function(x){
      ## 7:9 refer to x,y,z coordinates
      mat1<-as.numeric(mat.hic[x[1],7:9]); pat1<-as.numeric(pat.hic[x[1],7:9])   #two alleles of Location
      mat2<-as.numeric(mat.hic[x[2],7:9]); pat2<-as.numeric(pat.hic[x[2],7:9])   #two alleles of location1
      mat3<-as.numeric(mat.hic[x[3],7:9]); pat3<-as.numeric(pat.hic[x[3],7:9])   #two alleles of location2
      
      tmpd=min(c(
        sum(dist(rbind(mat1,mat2,mat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(mat1,mat2,pat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(mat1,pat2,mat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(mat1,pat2,pat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(pat1,mat2,mat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(pat1,mat2,pat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(pat1,pat2,mat3), method = "euclidean",upper = FALSE)),
        sum(dist(rbind(pat1,pat2,pat3), method = "euclidean",upper = FALSE))),
        na.rm = T)
      tmpd[is.infinite(tmpd)]<-NA
      return(tmpd);rm(tmpd)
    }
    )
  af4ctrl.eud<-cbind(af4ctrl.eud,eu.d)
  rm(eu.d)
}

rm(i)
colnames(af4ctrl.eud)<-c(colnames(af4ctrl),samples)
#GM12878+PBMC
af4ctrl.eud$ave.eud<-apply(af4ctrl.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{mean(x,na.rm = T)}})
af4ctrl.eud$min.eud<-apply(af4ctrl.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{ min(x,na.rm = T)}})
af4ctrl.eud$min.eud<-as.numeric(af4ctrl.eud$min.eud)
af4ctrl.eud$ave.eud<-as.numeric(af4ctrl.eud$ave.eud)

#-GM12878
af4ctrl.eud$ave.gm12878<-apply(af4ctrl.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{mean(x,na.rm = T)}})
af4ctrl.eud$min.gm12878<-apply(af4ctrl.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{ min(x,na.rm = T)}})
af4ctrl.eud$min.gm12878<-as.numeric(af4ctrl.eud$min.gm12878)
af4ctrl.eud$ave.gm12878<-as.numeric(af4ctrl.eud$ave.gm12878)
#PBMC
af4ctrl.eud$ave.pbmc<-apply(af4ctrl.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{mean(x,na.rm = T)}})
af4ctrl.eud$min.pbmc<-apply(af4ctrl.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{ min(x,na.rm = T)}})
af4ctrl.eud$min.pbmc<-as.numeric(af4ctrl.eud$min.pbmc)
af4ctrl.eud$ave.pbmc<-as.numeric(af4ctrl.eud$ave.pbmc)

af4ctrl.eud<-af4ctrl.eud[!is.na(af4ctrl.eud$ave.eud),]


#-----Comparisons of ctrl and AF4-----#
wilcox.test(as.numeric(as.matrix(af4part.eud[,gm12878])),as.numeric(as.matrix(af4ctrl.eud[,gm12878])))
wilcox.test(as.numeric(as.matrix(af4part.eud[,pbmc])),as.numeric(as.matrix(af4ctrl.eud[,pbmc])))

boxplot(as.numeric(as.matrix(af4part.eud[,c(gm12878,pbmc)])),
        as.numeric(as.matrix(af4ctrl.eud[sample(row.names(af4ctrl.eud),21),c(gm12878,pbmc)])))

t.test(af4part.eud$min.eud,af4ctrl.eud$min.eud)
t.test(af4part.eud$ave.eud,af4ctrl.eud$ave.eud)


af4boxp<-(rbind(af4ctrl.eud[,c("chr","start",gm12878,pbmc,"ave.eud","min.eud")],
                af4part.eud[,c("chr","start",gm12878,pbmc,"ave.eud","min.eud")]))

af4boxp$group<-rep(c('ctrl','complex_rearragements'),c(nrow(af4ctrl.eud),nrow(af4part.eud)))
af4boxp$loci<-paste(af4boxp$chr,round(af4boxp$start/10^6,0),sep = ".")
af4boxp$loci<-paste('chr',af4boxp$loci,'M',sep = "")
#remove duplicated loci
#af4boxp$loci[duplicated(af4boxp$loci)]<-paste(af4boxp$loci[duplicated(af4boxp$loci)],2,sep ='.')
af4boxp<-af4boxp[!duplicated(af4boxp$loci),]

library("dplyr")
af4boxp$chr<-as.numeric(af4boxp$chr)
af4boxp$start<-as.numeric(af4boxp$start)
af4boxp<-arrange(af4boxp,af4boxp$chr,af4boxp$start)

