library("magrittr")
library("scatterplot3d")
library("reshape2") #melt
library("ggplot2")
library(hash)
library(plyr)

#---read LAML fusions---#
#KIAA1530 —>UVSSA
#MLL->KMT2A
#ODZ1->TENM1(adjacent gene)
tfusions=read.table("supTable2.tumor_fusion_lists.txt",sep = "\t",header = T,stringsAsFactors = F)
#crch37(hg19):get info from ensembl
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
genetmp<-
  getBM(attributes = c("chromosome_name", "start_position","end_position", "external_gene_name"),
        filters = "external_gene_name",
        values = c(tfusions$GeneA,tfusions$GeneB),
        mart=mart)
genetmp<-genetmp[!grepl('_',genetmp$chromosome_name),]
genetmp$loc<-round((genetmp$start_position+genetmp$end_position)/2,0)

genetmp<-genetmp[,c('external_gene_name','chromosome_name','loc')]
colnames(genetmp)<-c('gene','chr','loc')

tfusions<-merge(merge(tfusions,genetmp,by.x = 'GeneA', by.y = 'gene'),
      genetmp,by.x = 'GeneB',by.y = 'gene')
tfusions<-tfusions[,c("GeneA","GeneB","chrA","loc.x","chrB","loc.y")]
colnames(tfusions)<-c("GeneA","GeneB","leftchr","leftloc","rightchr","rightloc")

#---Location each parters of SVs in dis(dataframe)
 #-left parter
tfusions$Location1<-
  apply(tfusions[,c("leftchr","leftloc")],1,function(x){
    any.chr<-as.numeric(x[1]); any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>dis$start & any.pos<dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )
 #right parter
tfusions$Location2<-
  apply(tfusions[,c("rightchr","rightloc")],1,function(x){
    any.chr<-as.character(x[1]); any.pos<-as.numeric(x[2])
    tmpnum<-which(dis$chr==any.chr & any.pos>dis$start & any.pos<dis$end)[1]
    return(tmpnum);rm(tmpnum) }
  )

# remove chr in 'chr1', change X to 23
#tfusions$leftchr[tfusions$leftchr=='X']<-23  #change X to 23
#tfusions$rightchr[tfusions$rightchr=='X']<-23  #change X to 23
tfusions<-tfusions[tfusions$leftchr!=tfusions$rightchr,]
tfusions.eud<-tfusions

for (i in 1:length(samples)){
  mat.hic=read.table(paste('../20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('../20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(tfusions[,c("Location1","Location2")],1,function(x){
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
  tfusions.eud<-cbind(tfusions.eud,eu.d)
  rm(eu.d)
}
rm(i)

colnames(tfusions.eud)<-c(colnames(tfusions),samples)
#tfusions.eud$min.eud<-as.numeric(tfusions.eud$min.eud)
#tfusions.eud$ave.eud<-as.numeric(tfusions.eud$ave.eud)

#GM12878+PBMC
tfusions.eud$ave.eud<-apply(tfusions.eud[,c(gm12878,pbmc)],1,function(x){mean(x,na.rm = T)})
tfusions.eud$min.eud<-apply(tfusions.eud[,c(gm12878,pbmc)],1,function(x){ min(x,na.rm = T)})  #generate inf value
#tfusions.eud$min.eud[is.infinite(tfusions.eud$min.eud)]<-'NA'  #remove inf value
tfusions.eud$min.eud<-as.numeric(tfusions.eud$min.eud)
tfusions.eud$ave.eud<-as.numeric(tfusions.eud$ave.eud)

#-GM12878
tfusions.eud$ave.gm12878<-apply(tfusions.eud[,gm12878],1,function(x){mean(x,na.rm = T)})
tfusions.eud$min.gm12878<-apply(tfusions.eud[,gm12878],1,function(x){ min(x,na.rm = T)})  #generate inf value
#tfusions.eud$min.gm12878[is.infinite(tfusions.eud$min.gm12878)]<-'NA'  #remove inf value
tfusions.eud$min.gm12878<-as.numeric(tfusions.eud$min.gm12878)
tfusions.eud$ave.gm12878<-as.numeric(tfusions.eud$ave.gm12878)
#PBMC
tfusions.eud$ave.pbmc<-apply(tfusions.eud[,pbmc],1,function(x){mean(x,na.rm = T)})
tfusions.eud$min.pbmc<-apply(tfusions.eud[,pbmc],1,function(x){ min(x,na.rm = T)})  #generate inf value
#tfusions.eud$min.pbmc[is.infinite(tfusions.eud$min.pbmc)]<-'NA'  #remove inf value
tfusions.eud$min.pbmc<-as.numeric(tfusions.eud$min.pbmc)
tfusions.eud$ave.pbmc<-as.numeric(tfusions.eud$ave.pbmc)

###########################################################
#-------#---each fusions combined with three random genes---#
###########################################################
nright<-sample(1:nrow(dis),34*6,replace=F)
tfusctrl<-as.data.frame(cbind(rep(c(tfusions$Location1,tfusions$Location2),3),nright))
colnames(tfusctrl)<-c("Location1","Location2")
which(dis$chr[tfusctrl$Location1]==dis$chr[tfusctrl$Location2]) #5,22,33行相同染色体
#tfusctrl$Location2[c(5,22,33)]<-c(147617,77933,96346) #随机生成3个数字替代
tfusctrl$Location2[c(57,64,74,79,118,129,161,193,194)]<-sample(1:nrow(dis),9,replace=F) #随机生成3个数字替代
#rm(nright)

for (i in 1:length(samples)){
  mat.hic=read.table(paste('../20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('../20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(tfusctrl,1,function(x){
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
  tfusctrl<-cbind(tfusctrl,eu.d)  #ȫ?ֻ???,tfusctrl????һ?μ??㴦????<<-????
}
rm(i,eu.d)

colnames(tfusctrl)<-c(colnames(tfusctrl)[1:2],samples)
tfusctrl$ave.eud<-apply(tfusctrl[,samples],1,function(x){mean(x,na.rm = T)})
tfusctrl$min.eud<-apply(tfusctrl[,samples],1,function(x){ min(x,na.rm = T)})  #generate inf value
tfusctrl$min.eud[is.infinite(tfusctrl$min.eud)]<-'NA'  #remove inf value
tfusctrl$min.eud<-as.numeric(tfusctrl$min.eud)
tfusctrl<-tfusctrl[tfusctrl$ave.eud!='NaN',]

#----tfusctrl.eud:remove rows with much NA----#
na5<-apply(tfusctrl,1,function(x){if(length(which(is.na(x)))<5){return(TRUE)}else{return(FALSE)}}) #romove rows with >=5 NAs
tfusctrl.eud<-tfusctrl[na5,c("Location1","Location2",gm12878,pbmc,'ave.eud','min.eud')]
rm(na5)
#GM12878+PBMC
tfusctrl.eud$ave.eud<-apply(tfusctrl.eud[,c(gm12878,pbmc)],1,function(x){mean(x,na.rm = T)})
tfusctrl.eud$min.eud<-apply(tfusctrl.eud[,c(gm12878,pbmc)],1,function(x){ min(x,na.rm = T)})  #generate inf value
tfusctrl.eud$min.eud[is.infinite(tfusctrl.eud$min.eud)]<-'NA'  #remove inf value
tfusctrl.eud$min.eud<-as.numeric(tfusctrl.eud$min.eud)
#-GM12878
tfusctrl.eud$ave.gm12878<-apply(tfusctrl.eud[,gm12878],1,function(x){mean(x,na.rm = T)})
tfusctrl.eud$min.gm12878<-apply(tfusctrl.eud[,gm12878],1,function(x){ min(x,na.rm = T)})  #generate inf value
tfusctrl.eud$min.gm12878[is.infinite(tfusctrl.eud$min.gm12878)]<-'NA'  #remove inf value
tfusctrl.eud$min.gm12878<-as.numeric(tfusctrl.eud$min.gm12878)
#PBMC
tfusctrl.eud$ave.pbmc<-apply(tfusctrl.eud[,pbmc],1,function(x){mean(x,na.rm = T)})
tfusctrl.eud$min.pbmc<-apply(tfusctrl.eud[,pbmc],1,function(x){ min(x,na.rm = T)})  #generate inf value
tfusctrl.eud$min.pbmc[is.infinite(tfusctrl.eud$min.pbmc)]<-'NA'  #remove inf value
tfusctrl.eud$min.pbmc<-as.numeric(tfusctrl.eud$min.pbmc)
#rm(tfusctrl)
##################################################################################
#---------------make leukemia+ctrl(based on tumorFusion database) for 1B---------#
##################################################################################
#---01---GM12878------#
group<-c(rep('fusion',length(tfusions.eud$ave.gm12878)),rep('ctrl',length(tfusctrl.eud$ave.gm12878)))
dist<-c(tfusions.eud$ave.gm12878,tfusctrl.eud$ave.gm12878)
ave<-data.frame(group,dist,stat=rep('Average',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dist)
####------Using Minimun Distance instead of average between samples-------##
group<-c(rep('fusion',length(tfusions.eud$min.gm12878)),rep('ctrl',length(tfusctrl.eud$min.gm12878)))
dist<-c(tfusions.eud$min.gm12878,tfusctrl.eud$min.gm12878)
mini<-data.frame(group,dist,stat=rep('Minimum',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dist)
tsigfusion.gm12878<-rbind(ave,mini);rm(ave,mini)
colnames(tsigfusion)<-c('Group','Distances','Stat')

#---02---PBMC------#
group<-c(rep('fusion',length(tfusions.eud$ave.pbmc)),rep('ctrl',length(tfusctrl.eud$ave.pbmc)))
dist<-c(tfusions.eud$ave.pbmc,tfusctrl.eud$ave.pbmc)
ave<-data.frame(group,dist,stat=rep('Average',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dist)
####------Using Minimun Distance instead of average between samples-------##
group<-c(rep('fusion',length(tfusions.eud$min.pbmc)),rep('ctrl',length(tfusctrl.eud$min.pbmc)))
dist<-c(tfusions.eud$min.pbmc,tfusctrl.eud$min.pbmc)
mini<-data.frame(group,dist,stat=rep('Minimum',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dist)
tsigfusion.pbmc<-rbind(ave,mini);rm(ave,mini)
colnames(tsigfusion)<-c('Group','Distances','Stat')
tsigfusion<-rbind(tsigfusion.gm12878,tsigfusion.pbmc)
tsigfusion$cell.type<-c(rep('GM12878',nrow(tsigfusion.gm12878)),rep('PBMC',nrow(tsigfusion.pbmc)))
colnames(tsigfusion)<-c('Group','Distances','Stat','Cell.type')
#----Make sigfusion finished----#

#---Plot fig 1B(New)
#pdf("../../../figures/figure.1B.tumorfusion.pdf",height =5,width =7)
ggplot(tsigfusion, aes(x=Group, y=Distances)) + 
  geom_boxplot(aes(fill = Group), alpha = 0.5,show.legend = F) + 
  facet_grid(Cell.type~Stat,scales = "free")
dev.off()
library(lsr)
wilcox.test(tfusctrl.eud$ave.eud,tfusions.eud$ave.eud)$p.value #0.002
wilcox.test(tfusctrl.eud$min.eud,tfusions.eud$min.eud)$p.value #0.0008
wilcox.test(tfusctrl.eud$ave.gm12878,tfusions.eud$ave.gm12878)$p.value #0.0018
wilcox.test(tfusctrl.eud$min.gm12878,tfusions.eud$min.gm12878)$p.value #0.0120
wilcox.test(tfusctrl.eud$ave.pbmc,tfusions.eud$ave.pbmc)$p.value #0.0049
wilcox.test(tfusctrl.eud$min.pbmc,tfusions.eud$min.pbmc)$p.value #0.0048

