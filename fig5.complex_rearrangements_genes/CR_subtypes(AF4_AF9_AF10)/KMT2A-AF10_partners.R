library("biomaRt") #retrive information for ensembl
library(hash) #perl hash
library("reshape2") #melt

#-----Leukemia fusion genes-------#
af10gene<-read.table("complex_rearrangements/KMT2A-AF10_genes.txt",header = F,sep = "\t",stringsAsFactors = F)

#crch37(hg19):get info from ensembl
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
af10tmp<-
  getBM(attributes = c("external_gene_name", "chromosome_name", "start_position"),
        filters = "external_gene_name",
        values = af10gene,
        mart=mart)
af10tmp<-af10tmp[!grepl('_',af10tmp$chromosome_name),]
setdiff(af10gene[[1]],af10tmp$external_gene_name) #unmatched ID
af10gene<-af10tmp
rm(af10tmp)

af10loci=read.table("complex_rearrangements/KMT2A-AF10_loci.txt",sep = "\t",skip = 1)
colnames(af10gene)<-c('gene','chr','start')
colnames(af10loci)<-c('gene','chr','start')
af10part=rbind(af10gene,af10loci)

#-Becase of AF10 is on chromosome4, we remove genes on chr10(AF10) and chr11(KMT2A).
#-We got 21 loci.
af10part<-af10part[af10part$chr!=10,]
af10part<-af10part[af10part$chr!=11,]

#Location each parters of KMT2A-AF10 partners.
#Add = in > and <
af10part$Location<-
  apply(af10part[,c("chr","start")],1,function(x){
    any.chr<-as.numeric(x[1]); any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>=dis$start & any.pos<=dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )

af10part$Location1<-rep(104479,nrow(af10part)) #Add KMT2A
af10part$Location2<-rep(92878,nrow(af10part))  #Add MLLT10(AF10)

#-----sum of sides of triangle----#
af10part.eud<-af10part
for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(af10part[,c("Location","Location1","Location2")],1,function(x){
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
  af10part.eud<-cbind(af10part.eud,eu.d)
  rm(eu.d)
}

rm(i)
colnames(af10part.eud)<-c(colnames(af10part),samples)
#GM12878+PBMC
af10part.eud$ave.eud<-apply(af10part.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{mean(x,na.rm = T)}})
af10part.eud$min.eud<-apply(af10part.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{ min(x,na.rm = T)}})
af10part.eud$min.eud<-as.numeric(af10part.eud$min.eud)
af10part.eud$ave.eud<-as.numeric(af10part.eud$ave.eud)

#-GM12878
af10part.eud$ave.gm12878<-apply(af10part.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{mean(x,na.rm = T)}})
af10part.eud$min.gm12878<-apply(af10part.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{ min(x,na.rm = T)}})
af10part.eud$min.gm12878<-as.numeric(af10part.eud$min.gm12878)
af10part.eud$ave.gm12878<-as.numeric(af10part.eud$ave.gm12878)
#PBMC
af10part.eud$ave.pbmc<-apply(af10part.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{mean(x,na.rm = T)}})
af10part.eud$min.pbmc<-apply(af10part.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{ min(x,na.rm = T)}})
af10part.eud$min.pbmc<-as.numeric(af10part.eud$min.pbmc)
af10part.eud$ave.pbmc<-as.numeric(af10part.eud$ave.pbmc)

#----------Controls of KMT2A-MLLT10(AF10)-------------#
#---Make controls of AF10 fusions---#
cnum<-sort(sample(nrow(dis),300,replace = F))
af10ctrl<-dis[cnum,]
af10ctrl$Location<-cnum
rm(cnum)
af10ctrl$Location1<-rep(104479,nrow(af10ctrl)) #Add KMT2A
af10ctrl$Location2<-rep(92878,nrow(af10ctrl))  #Add MLLT10(AF10)
af10ctrl<-af10ctrl[!af10ctrl$chr %in% c(10,11,'X'),] #remove chr10,chr11 and chrX

#-Distances
af10ctrl.eud<-af10ctrl
for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(af10ctrl[,c("Location","Location1","Location2")],1,function(x){
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
  af10ctrl.eud<-cbind(af10ctrl.eud,eu.d)
  rm(eu.d)
}

rm(i)
colnames(af10ctrl.eud)<-c(colnames(af10ctrl),samples)
#GM12878+PBMC
af10ctrl.eud$ave.eud<-apply(af10ctrl.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{mean(x,na.rm = T)}})
af10ctrl.eud$min.eud<-apply(af10ctrl.eud[,c(gm12878,pbmc)],1,function(x){if(length(which(is.na(x)))>4){NA}else{ min(x,na.rm = T)}})
af10ctrl.eud$min.eud<-as.numeric(af10ctrl.eud$min.eud)
af10ctrl.eud$ave.eud<-as.numeric(af10ctrl.eud$ave.eud)

#-GM12878
af10ctrl.eud$ave.gm12878<-apply(af10ctrl.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{mean(x,na.rm = T)}})
af10ctrl.eud$min.gm12878<-apply(af10ctrl.eud[,gm12878],1,function(x){if(length(which(is.na(x)))>2){NA}else{ min(x,na.rm = T)}})
af10ctrl.eud$min.gm12878<-as.numeric(af10ctrl.eud$min.gm12878)
af10ctrl.eud$ave.gm12878<-as.numeric(af10ctrl.eud$ave.gm12878)
#PBMC
af10ctrl.eud$ave.pbmc<-apply(af10ctrl.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{mean(x,na.rm = T)}})
af10ctrl.eud$min.pbmc<-apply(af10ctrl.eud[,pbmc],1,function(x){if(length(which(is.na(x)))>3){NA}else{ min(x,na.rm = T)}})
af10ctrl.eud$min.pbmc<-as.numeric(af10ctrl.eud$min.pbmc)
af10ctrl.eud$ave.pbmc<-as.numeric(af10ctrl.eud$ave.pbmc)

af10ctrl.eud<-af10ctrl.eud[!is.na(af10ctrl.eud$ave.eud),]

#-----Comparisons of ctrl and af10-----#
t.test(af10part.eud$min.eud,af10ctrl.eud$min.eud)
t.test(af10part.eud$ave.eud,af10ctrl.eud$ave.eud)


af10boxp<-(rbind(af10ctrl.eud[,c("chr","start",gm12878,pbmc,"ave.eud","min.eud")],
                 af10part.eud[,c("chr","start",gm12878,pbmc,"ave.eud","min.eud")]))
af10boxp$group<-rep(c('ctrl','complex_rearragements'),c(nrow(af10ctrl.eud),nrow(af10part.eud)))
af10boxp$loci<-paste(af10boxp$chr,round(af10boxp$start/10^6,0),sep = ".")
af10boxp$loci<-paste('chr',af10boxp$loci,'M',sep = "")
#remove duplicated loci
#af10boxp$loci[duplicated(af10boxp$loci)]<-paste(af10boxp$loci[duplicated(af10boxp$loci)],2,sep ='.')
af10boxp<-af10boxp[!duplicated(af10boxp$loci),]

library("dplyr")
af10boxp$chr<-as.numeric(af10boxp$chr)
af10boxp$start<-as.numeric(af10boxp$start)
af10boxp<-arrange(af10boxp,af10boxp$chr,af10boxp$start)

af10tmp<-af10boxp[af10boxp$chr %in% unique(af10part.eud$chr),c(1:31,34,35)]
af10tmp2<-melt(af10tmp,id.vars = c("chr","start","loci","group"))
af10tmp2<-arrange(af10tmp2,af10tmp2$chr,af10tmp2$start)
af10tmp2<-af10tmp2[,3:6]

colnames(af10tmp2)<-c('loci','group','sample','value')
af10tmp2$loci<-as.character(af10tmp2$loci)
af10tmp2$loci<-factor(af10tmp2$loci,level=unique(af10tmp2$loci))
#-fig7A
ggplot(af10tmp2, aes(x=loci, y=value)) + 
  geom_boxplot(aes(fill = group), alpha = 0.5,show.legend = F)+
  theme(axis.text.x = element_text(face = "bold",size = 6, vjust = 0.5, hjust = 0.5, angle = 90))
#-fig7B
af10tmp3<-af10boxp[,c("ave.eud","min.eud","group")]
af10tmp3<-melt(af10tmp3,id.vars = c("group"))
colnames(af10tmp3)<-c('Group','Stat','Distances')

ggplot(af10tmp3, aes(x=Group, y=Distances)) + 
  geom_boxplot(aes(fill = Group), alpha = 0.5,show.legend = F) + 
  facet_grid(~Stat,scales = "free")







