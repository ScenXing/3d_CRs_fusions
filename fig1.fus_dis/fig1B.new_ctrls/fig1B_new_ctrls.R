###################################################
#######----Random positions as control---######
###################################################
#nleft<-sample(1:nrow(dis),58,replace=F)
nright<-sample(1:nrow(dis),58*2,replace=F)

ctrlfus<-as.data.frame(cbind(c(dfusions$Location1,dfusions$Location2),nright))
colnames(ctrlfus)<-c("Location1","Location2")
which(dis$chr[ctrlfus$Location1]==dis$chr[ctrlfus$Location2]) #22,79 same
ctrlfus$Location2[c(55,70,95)]<-c(121600,110492,32953) #random
#rm(nright)

for (i in 1:length(samples)){
  mat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.mat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  #colnames(mat.hic) <- c('chr','start','end','chr2','start2','end2','x','y','z')
  pat.hic=read.table(paste('20kb_bin/',samples[i],'.3dg.pat.20k.bed',sep=""),sep="\t",header=F,stringsAsFactors=FALSE)
  
  eu.d<-
    apply(ctrlfus,1,function(x){
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
  ctrlfus<-cbind(ctrlfus,eu.d)
}
rm(i,eu.d)

colnames(ctrlfus)<-c(colnames(ctrlfus)[1:2],samples)
ctrlfus$ave.eud<-apply(ctrlfus[,samples],1,function(x){mean(x,na.rm = T)})
ctrlfus$min.eud<-apply(ctrlfus[,samples],1,function(x){ min(x,na.rm = T)})  #generate inf value
ctrlfus$min.eud[is.infinite(ctrlfus$min.eud)]<-'NA'  #remove inf value
ctrlfus$min.eud<-as.numeric(ctrlfus$min.eud)
ctrlfus<-ctrlfus[ctrlfus$ave.eud!='NaN',]

#----ctrlfus.eud:remove rows with much NA----#
#--ctrlfus.eud:33 coloum, location1, location2, 11 GM12878 and 18 PBMC.
na5<-apply(ctrlfus,1,function(x){if(length(which(is.na(x)))<5){return(TRUE)}else{return(FALSE)}}) #romove rows with >=5 NAs
ctrlfus.eud<-ctrlfus[na5,c("Location1","Location2",gm12878,pbmc,'ave.eud','min.eud')]
rm(na5)
#GM12878+PBMC
ctrlfus.eud$ave.eud<-apply(ctrlfus.eud[,c(gm12878,pbmc)],1,function(x){mean(x,na.rm = T)})
ctrlfus.eud$min.eud<-apply(ctrlfus.eud[,c(gm12878,pbmc)],1,function(x){ min(x,na.rm = T)})  #generate inf value
ctrlfus.eud$min.eud[is.infinite(ctrlfus.eud$min.eud)]<-'NA'  #remove inf value
ctrlfus.eud$min.eud<-as.numeric(ctrlfus.eud$min.eud)
#-GM12878
ctrlfus.eud$ave.gm12878<-apply(ctrlfus.eud[,gm12878],1,function(x){mean(x,na.rm = T)})
ctrlfus.eud$min.gm12878<-apply(ctrlfus.eud[,gm12878],1,function(x){ min(x,na.rm = T)})  #generate inf value
ctrlfus.eud$min.gm12878[is.infinite(ctrlfus.eud$min.gm12878)]<-'NA'  #remove inf value
ctrlfus.eud$min.gm12878<-as.numeric(ctrlfus.eud$min.gm12878)
#PBMC
ctrlfus.eud$ave.pbmc<-apply(ctrlfus.eud[,pbmc],1,function(x){mean(x,na.rm = T)})
ctrlfus.eud$min.pbmc<-apply(ctrlfus.eud[,pbmc],1,function(x){ min(x,na.rm = T)})  #generate inf value
ctrlfus.eud$min.pbmc[is.infinite(ctrlfus.eud$min.pbmc)]<-'NA'  #remove inf value
ctrlfus.eud$min.pbmc<-as.numeric(ctrlfus.eud$min.pbmc)
#rm(ctrlfus)

#---------------new ctrl finished-----------------#
#---------------make leukemia+ctrl(sigfusion) for 1B---------#
#---01---GM12878------#
group<-c(rep('fusion',length(dfusions.eud$ave.gm12878)),rep('ctrl',length(ctrlfus.eud$ave.gm12878)))
dis<-c(dfusions.eud$ave.gm12878,ctrlfus.eud$ave.gm12878)
ave<-data.frame(group,dis,stat=rep('Average',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dis)
####------Using Minimun Distance instead of average between samples-------##
group<-c(rep('fusion',length(dfusions.eud$min.gm12878)),rep('ctrl',length(ctrlfus.eud$min.gm12878)))
dis<-c(dfusions.eud$min.gm12878,ctrlfus.eud$min.gm12878)
mini<-data.frame(group,dis,stat=rep('Minimum',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dis)
sigfusion.gm12878<-rbind(ave,mini);rm(ave,mini)
colnames(sigfusion)<-c('Group','Distances','Stat')

#---02---PBMC------#
group<-c(rep('fusion',length(dfusions.eud$ave.pbmc)),rep('ctrl',length(ctrlfus.eud$ave.pbmc)))
dis<-c(dfusions.eud$ave.pbmc,ctrlfus.eud$ave.pbmc)
ave<-data.frame(group,dis,stat=rep('Average',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dis)
####------Using Minimun Distance instead of average between samples-------##
group<-c(rep('fusion',length(dfusions.eud$min.pbmc)),rep('ctrl',length(ctrlfus.eud$min.pbmc)))
dis<-c(dfusions.eud$min.pbmc,ctrlfus.eud$min.pbmc)
mini<-data.frame(group,dis,stat=rep('Minimum',length(group)))  # data.fram[1-group + 2-distance]
rm(group,dis)
sigfusion.pbmc<-rbind(ave,mini);rm(ave,mini)
colnames(sigfusion)<-c('Group','Distances','Stat')
sigfusion<-rbind(sigfusion.gm12878,sigfusion.pbmc)
sigfusion$cell.type<-c(rep('GM12878',nrow(sigfusion.gm12878)),rep('PBMC',nrow(sigfusion.pbmc)))
colnames(sigfusion)<-c('Group','Distances','Stat','Cell.type')
#----Make sigfuion finished----#
#----Add solid tumor fusions---#
ofusionStat<-rbind(ofusionStat,sigfusion)
ofusionStat$Group<-factor(ofusionStat$Group,levels = c('ctrl','Other fusion','fusion'))

#---Plot fig 1B(New)
#pdf("../../figures/figure.1B_newctrl.pdf",height =5,width =7)
ggplot(ofusionStat, aes(x=Group, y=Distances)) + 
  geom_boxplot(aes(fill = Group), alpha = 0.5,show.legend = F) + 
  facet_grid(Stat ~ Cell.type,scales = "free")
#dev.off()
wilcox.test(ctrlfus.eud$ave.eud,dfusions.eud$ave.eud)$p.value
wilcox.test(ctrlfus.eud$min.eud,dfusions.eud$min.eud)$p.value
wilcox.test(ctrlfus.eud$ave.gm12878,dfusions.eud$ave.gm12878)$p.value
wilcox.test(ctrlfus.eud$min.gm12878,dfusions.eud$min.gm12878)$p.value
wilcox.test(ctrlfus.eud$ave.pbmc,dfusions.eud$ave.pbmc)$p.value
wilcox.test(ctrlfus.eud$min.pbmc,dfusions.eud$min.pbmc)$p.value
