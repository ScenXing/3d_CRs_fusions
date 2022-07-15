library("biomaRt") #retrive information for ensembl
library(hash) #perl hash
library("reshape2") #melt

#-----Leukemia fusion genes-------#
af9gene<-read.table("../figure7.complex_rearrangements_genes//KMT2A-AF9_genes.txt",header = F,sep = "\t",stringsAsFactors = F)

#crch37(hg19):get info from ensembl
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
af9tmp<-
  getBM(attributes = c("external_gene_name", "chromosome_name", "start_position"),
        filters = "external_gene_name",
        values = af9gene,
        mart=mart)
af9tmp<-af9tmp[!grepl('_',af9tmp$chromosome_name),]
setdiff(af9gene[[1]],af9tmp$external_gene_name) #unmatched ID
af9gene<-af9tmp
rm(af9tmp)

af9loci=read.table("../figure7.complex_rearrangements_genes/KMT2A-AF9_loci.txt",sep = "\t",skip = 1)
colnames(af9gene)<-c('gene','chr','start')
colnames(af9loci)<-c('gene','chr','start')
af9part=rbind(af9gene,af9loci)

#-Becase of af9 is on chromosome 9, we remove genes on chr9(AF9) and chr11(KMT2A).
#-We got 21 loci.
af9part<-af9part[af9part$chr!=9,]
af9part<-af9part[af9part$chr!=11,]

#Location each parters of KMT2A-af9 partners.
#Add = in > and <
af9part$Location<-
  apply(af9part[,c("chr","start")],1,function(x){
    any.chr<-as.numeric(x[1]); any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>=dis$start & any.pos<=dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )

af9part$Location1<-rep(104479,nrow(af9part)) #Add KMT2A
af9part$Location2<-rep(85743,nrow(af9part))  #Add MLLT3(af9)

#-----sum of sides of triangle----#
af9part.eud<-af9part
for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(af9part[,c("Location","Location1","Location2")],1,function(x){
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
  af9part.eud<-cbind(af9part.eud,eu.d)
  rm(eu.d)
}

rm(i)
colnames(af9part.eud)<-c(colnames(af9part),samples)
#GM12878+PBMC
af9part.eud$ave.eud<-apply(af9part.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{mean(x,na.rm = T)}})
af9part.eud$min.eud<-apply(af9part.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{ min(x,na.rm = T)}})
af9part.eud$min.eud<-as.numeric(af9part.eud$min.eud)
af9part.eud$ave.eud<-as.numeric(af9part.eud$ave.eud)

#-GM12878
af9part.eud$ave.gm12878<-apply(af9part.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{mean(x,na.rm = T)}})
af9part.eud$min.gm12878<-apply(af9part.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{ min(x,na.rm = T)}})
af9part.eud$min.gm12878<-as.numeric(af9part.eud$min.gm12878)
af9part.eud$ave.gm12878<-as.numeric(af9part.eud$ave.gm12878)
#PBMC
af9part.eud$ave.pbmc<-apply(af9part.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{mean(x,na.rm = T)}})
af9part.eud$min.pbmc<-apply(af9part.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{ min(x,na.rm = T)}})
af9part.eud$min.pbmc<-as.numeric(af9part.eud$min.pbmc)
af9part.eud$ave.pbmc<-as.numeric(af9part.eud$ave.pbmc)

#----------Controls of KMT2A-MLLT3(af9)-------------#
#---Make controls of af9 fusions---#
cnum<-sort(sample(nrow(dis),300,replace = F))
af9ctrl<-dis[cnum,]
af9ctrl$Location<-cnum
rm(cnum)
af9ctrl$Location1<-rep(104479,nrow(af9ctrl)) #Add KMT2A
af9ctrl$Location2<-rep(85743,nrow(af9ctrl))  #Add MLLT3(af9)
af9ctrl<-af9ctrl[!af9ctrl$chr %in% c(9,11,'X'),] #remove chr9,chr11 and chrX

#-Distances
af9ctrl.eud<-af9ctrl
for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(af9ctrl[,c("Location","Location1","Location2")],1,function(x){
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
  af9ctrl.eud<-cbind(af9ctrl.eud,eu.d)
  rm(eu.d)
}

rm(i)
colnames(af9ctrl.eud)<-c(colnames(af9ctrl),samples)
#GM12878+PBMC
af9ctrl.eud$ave.eud<-apply(af9ctrl.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{mean(x,na.rm = T)}})
af9ctrl.eud$min.eud<-apply(af9ctrl.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{ min(x,na.rm = T)}})
af9ctrl.eud$min.eud<-as.numeric(af9ctrl.eud$min.eud)
af9ctrl.eud$ave.eud<-as.numeric(af9ctrl.eud$ave.eud)

#-GM12878
af9ctrl.eud$ave.gm12878<-apply(af9ctrl.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{mean(x,na.rm = T)}})
af9ctrl.eud$min.gm12878<-apply(af9ctrl.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{ min(x,na.rm = T)}})
af9ctrl.eud$min.gm12878<-as.numeric(af9ctrl.eud$min.gm12878)
af9ctrl.eud$ave.gm12878<-as.numeric(af9ctrl.eud$ave.gm12878)
#PBMC
af9ctrl.eud$ave.pbmc<-apply(af9ctrl.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{mean(x,na.rm = T)}})
af9ctrl.eud$min.pbmc<-apply(af9ctrl.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{ min(x,na.rm = T)}})
af9ctrl.eud$min.pbmc<-as.numeric(af9ctrl.eud$min.pbmc)
af9ctrl.eud$ave.pbmc<-as.numeric(af9ctrl.eud$ave.pbmc)

af9ctrl.eud<-af9ctrl.eud[!is.na(af9ctrl.eud$ave.eud),]

#-----Comparisons of ctrl and af9-----#
t.test(af9part.eud$min.eud,af9ctrl.eud$min.eud)
t.test(af9part.eud$ave.eud,af9ctrl.eud$ave.eud)


af9boxp<-(rbind(af9ctrl.eud[,c("chr","start",gm12878,pbmc,"ave.eud","min.eud")],
                 af9part.eud[,c("chr","start",gm12878,pbmc,"ave.eud","min.eud")]))

af9boxp$group<-rep(c('ctrl','complex_rearragements'),c(nrow(af9ctrl.eud),nrow(af9part.eud)))
af9boxp$loci<-paste(af9boxp$chr,round(af9boxp$start/10^6,0),sep = ".")
af9boxp$loci<-paste('chr',af9boxp$loci,'M',sep = "")
#remove duplicated loci
#af9boxp$loci[duplicated(af9boxp$loci)]<-paste(af9boxp$loci[duplicated(af9boxp$loci)],2,sep ='.')
af9boxp<-af9boxp[!duplicated(af9boxp$loci),]

library("dplyr")
af9boxp$chr<-as.numeric(af9boxp$chr)
af9boxp$start<-as.numeric(af9boxp$start)
af9boxp<-arrange(af9boxp,af9boxp$chr,af9boxp$start)











