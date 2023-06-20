##################################begin definitions section######################################################
makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}
#my modified heatmap.2 source code, with OLO
source("/media/rs_DATA/s3/20210111_varner_actri_voucher/heatmap.3_original.R")
####################################end definitions section######################################################


features<-c("HLADR","CXCR4","CD28","CXCR3","CXCR5","CD127","CD69","TIM3","CD45RO","PD1","CTLA4","CD25","ICOS","LAG3","CD62L")
tetras<-c("12-20","13-21","CP11","CP18")
#these subjects are ordered by cohort
subjects<-c("NBD-BRI001","NBD-BRI002","NBD-BRI003","NBD-BRI004","NBD-BRI005",
            "T030","T031","T033","T034","T035",
            "RAD007","RAD009","RAD010","RAD011",
            "TNET006","TNET007","TNET008","TNET009","TNET010","TNET011","TNET012","TNET013","TNET014","TNET015")
cohorts<-c("NBD","ES","JD","MON")

longcohorts<-rep(cohorts,c(5,5,4,10))
names(longcohorts)<-subjects


#define "clean" subjects, etc.   
ccohorts<-c("NBD","ES","JD")
csubjects<-c("NBD-BRI001","NBD-BRI002","NBD-BRI003","NBD-BRI004","NBD-BRI005",
            "T030","T031","T033","T034","T035",
            "RAD007","RAD009","RAD010","RAD011")

ncs<-length(csubjects)
clongcohorts<-rep(ccohorts,c(5,5,5))

colors3<-c("grey20","red","pink")
names(colors3)<-ccohorts

colors7<-c("grey70","darkolivegreen","brown","purple","magenta","blue","cyan","grey20","orange","red")
clcolors<-makeTransparent(colors7,alpha=0.08) #this covers A, Aa under one color

colors5<-c("grey70","deepskyblue4","deeppink3","goldenrod1","green")


load("tetras_V.RData") #loads V
names(V)
# [1] "UMAP 1"        "UMAP 2"        "HLADR"         "CXCR4"        
# [5] "CD28"          "CXCR3"         "CXCR5"         "CD127"        
# [9] "CD69"          "TIM3"          "CD45RO"        "PD1"          
#[13] "CTLA4"         "CD25"          "ICOS"          "LAG3"         
#[17] "CD62L"         "tetraPos"      "tetraID"       "subjectID"    
#[21] "cohortID"      "qTetramer"     "qtetraPos"     "UMAP.1"       
#[25] "UMAP.2"        "classes"       "classlabels"   "pointcolors"  
#[29] "quantTetramer" "tetraclasses" 


cl<-unique(V$tetraclasses)
ncl<-length(cl) #4
rge<-range(V[,1:2])
labels<-paste0("T",1:ncl)

pdf("tUMAP.pdf")
   plot(V[,1:2],pch=NA,xlim=rge,ylim=rge,main="") #tUMAP
   points(V[,1:2],pch=19,cex=.2,col=makeTransparent(colors5[V$tetraclasses],alpha=0.8))

   rge<-range(V[,25:26]) #global UMAP
   plot(V[,25:26],pch=NA,xlim=rge,ylim=rge,main="")
   points(V[,25:26],pch=19,cex=.2,col=colors5[V$tetraclasses])
dev.off()

df.mass<-data.frame(table(V$cohortID))
mass<-df.mass$Freq
names(mass)<-df.mass$Var1
mmass<-median(mass)
mscale<-1/sqrt(mass/mmass)
pdf("tUMAPs_by_cohort.pdf",width=12,height=4)
   par(mar=c(4, 1, 2, 1), mfrow=c(1,3), oma = c(4, 4, 2, 2),pty="s")
   rge<-range(V[,1:2])
   for (t in ccohorts) {
      cat(t,"\n")
      plot(V[V$cohortID==t,1:2],pch=19,cex=.2*mscale[t],xlim=rge,ylim=rge,main=t,col=colors5[V$tetraclasses[V$cohortID==t]])
   }
dev.off()


newfeatures<-c("CD25","ICOS","LAG3","PD1","CD69","HLADR","CD127","CD62L","CD45RO","CD28","TIM3","CXCR4","CXCR3","CXCR5","CTLA4")
#ridgeplots please
trcolors<-makeTransparent(colors5,alpha=0.8)
names(trcolors)<-labels
pdf("tetra_cl_ridgeplot.pdf",width=16,height=16)
par(mar=c(4, 1, 2, 1), mfrow=c(4,4), oma = c(4, 4, 2, 2),pty="s")

   #rowFtests
   partition<-labels[V$tetraclasses]
   part.df<-as.data.frame(table(partition))
   nbig<-dim(part.df)[1]
   cls<-part.df$partition
   x<-V[,1:15+2]
   #the following is good for F test (ANOVA), if desired. Here we just want to plot all genes and their distribution in clusters 
#  Ftest<-colFtests(x=data.matrix(x[,features]),fac=factor(x$classes,levels=bigcl))
#  Fstat<-Ftest$statistic
#  o<-order(-Fstat)
#  Ftest[o,]
   #reshape data into a list
   for (f in newfeatures) {
      dat.list<-list()
      for (cl in cls) dat.list[[cl]]<-x[partition==cl,f]
      #get density values
      vals <- Map(function(x, g, i) { with( density(x,n=512,bw=.16), data.frame(x,y=y+(i-1), g) ) }, dat.list, names(dat.list), seq_along(dat.list))
      #now plot ridgeline plot
      xrange <- range(unlist(lapply(vals, function(d) range(d$x))))
      yrange <- range(unlist(lapply(vals, function(d) range(d$y))))
      #plot(0,0, type="n", xlim=c(-.4,4.), ylim=yrange, yaxt="n", ylab="", xlab="standardized intensity",main=f)
      plot(0,0, type="n", xlim=c(-1,14), ylim=c(0,5), yaxt="n", ylab="", xlab="standardized intensity",main=f)
      abline(v=0,col="grey")
      for (cl in rev(cls)) { #this guarantees that the clusters are picked in the right order on the y axis of ridgeplot
         #with(vals[[cl]], polygon(x, y, col=trcolors[transtable.list[[t]][cl]],border=NA))
         with(vals[[cl]], polygon(x, y, col=trcolors[cl],border=NA))
      }
      axis(2, at=seq_along(dat.list)-1, names(dat.list))
   }
   #while(!par('page')) plot.new()
dev.off()


cl_signatures<-list()
for (t in ccohorts) {
   cat(t,"\n")
   comms<-sort(unique(partition))
   Y<-V[,1:15+2]

   #store cluster expression signatures
   cl_means<-matrix(0,nrow=length(newfeatures),ncol=length(comms),dimnames=list(newfeatures,comms))
   for (l in comms) cl_means[features,l]<- colMeans(Y[partition==l & V$cohortID==t,features],na.rm=TRUE)
   cl_means[is.na(cl_means)]<-0 #if no presence, set to 0
   cl_signatures[[t]]<-cl_means
}

save(cl_signatures,file="tetra_CD4_signatures.RData")

nf<-length(features)
cl_max<-rep(0,nf)
for (t in ccohorts) cl_max<-pmax(cl_max,apply(cl_signatures[[t]],1,max))
#cl_max<-pmax(cl_max,2)

#add heatmap of expression signatures
pdf("tetra_CD4_heatmaps.pdf")
   for (t in ccohorts) {
      cl_means<-cl_signatures[[t]]
      cl_means[cl_means<0]<-0 #set negative values to 0 for purpose of visualization
      cl_means<-(cl_means/cl_max)#^0.5 #gamma factor for visualization
      #this heatmap has olo
      ColSideColors<-trcolors[colnames(cl_means)]
      #colnames(cl_means)<-paste0(t,"_",colnames(cl_means))
      heatmap.2.1(cl_means[newfeatures,],key=FALSE,Rowv=FALSE,Colv=FALSE,margins=c(12,8),cexCol=2,cexRow=1.2,lhei=c(0.05,5),dendrogram="none",trace="none",ColSideColors=ColSideColors,col=hcl.colors(50,"Reds",rev=TRUE))
      title(sub=t,line=2,cex.sub=2)
   }
dev.off()


#now cluster fractions by cohort
#now just red cells per cluster
ntc<-length(comms) #4
redCountSummary<-matrix(0,nrow=ntc,ncol=ncs) #ncs is 15
rownames(redCountSummary)<-comms
colnames(redCountSummary)<-csubjects
for (s in csubjects) {
   for (icl in 1:ntc) {
      cl<-comms[icl]
      redCountSummary[cl,s]<-sum(V$subjectID==s & V$tetraclasses==icl,na.rm=TRUE)
   }
}
fractionRedCountSummary<-t(t(redCountSummary)/colSums(redCountSummary))

cohortRedCountSummary<-matrix(0,nrow=ntc,ncol=3)
rownames(cohortRedCountSummary)<-comms
colnames(cohortRedCountSummary)<-ccohorts
for (s in ccohorts) {
   for (icl in 1:ntc) {
      cl<-comms[icl]
      cohortRedCountSummary[cl,s]<-sum(V$cohortID==s & V$tetraclasses==icl,na.rm=TRUE)
   }
}
cohortFractionRedCountSummary<-t(t(cohortRedCountSummary)/colSums(cohortRedCountSummary))

pdf("tetra_fractions_by_cluster.pdf",width=10,height=5)
par(mfrow=c(1,2),pty="s",mar=c(8, 1, 2, 1))
barplot(fractionRedCountSummary[,c(1:5,11:14,6:10)],col=colors5,las=2,border=NA,main=expression(tetramer^"+"~cells),ylab="population fraction")
barplot(cohortFractionRedCountSummary[,c(1,3,2)],col=colors5,las=2,border=NA,main=expression(tetramer^"+"~cells),ylab="population fraction")
dev.off()


pdf("tetra_cluster_fractions.pdf",width=5,height=4)
dx<-0:2
dx<-c(0,2,1) #this order was requested
names(dx)<-ccohorts
barplot(t(cohortFractionRedCountSummary[,c(1,3,2)]),beside=TRUE,border=NA,col=colors3[c(1,3,2)],ylab="population fraction",xlab="cluster",ylim=c(0,.9))
legend("topright",legend=ccohorts[c(1,3,2)],pch=15,col=colors3[c(1,3,2)])
for (coh in ccohorts) {
   fr<-fractionRedCountSummary[,clongcohorts==coh]
   for (i in 1:ntc) {
      x<-(i-1)*4+1+dx[coh]+1/2+rnorm(dim(fr)[2],sd=.1)
      y<-fr[i,]
      points(x,y,cex=0.5,col="grey")
   }
}
box()
dev.off()
save(redCountSummary,fractionRedCountSummary,cohortRedCountSummary,cohortFractionRedCountSummary,file="tetra_Summaries.RData")

library(genefilter)
#get p-values for comparison between cohorts and one cluster at a time, T1-T4
p<-matrix(1,nrow=4,ncol=3) 
colnames(p)<-c("NBD_ES","NBD_JD","ES_JD")
rownames(p)<-rownames(fractionRedCountSummary)
t<-p #for size and names
fct<-factor(rep(0:1,each=5))
#NBD_ES
x<-fractionRedCountSummary[,1:10]
x<-log(x/(1-x)) #compositional variables
rtt<-rowttests(x,fac=fct,tstatOnly=FALSE)
t[,1]<-rtt$statistic
p[,1]<-rtt$p.value
#NBD_JD
x<-fractionRedCountSummary[,c(1:5,11:14)]
x<-log(x/(1-x)) #compositional variables
rtt<-rowttests(x,fac=factor(c(0,0,0,0,0,1,1,1,1)),tstatOnly=FALSE)
t[,2]<-rtt$statistic
p[,2]<-rtt$p.value
#ES_JD
x<-fractionRedCountSummary[,6:14]
x<-log(x/(1-x)) #compositional variables
rtt<-rowttests(x,fac=factor(c(0,0,0,0,0,1,1,1,1)),tstatOnly=FALSE)
t[,3]<-rtt$statistic
p[,3]<-rtt$p.value
write.csv(p,file="p.csv")

#get Bonferroni-adjusted p
p.Bon<-pmin(p*9,1)
write.csv(p.Bon,file="p.Bon.csv")

p
#       NBD_ES      NBD_JD      ES_JD
#T1 0.06902475 0.003853443 0.01508095
#T2 0.01187061 0.015374247 0.29988370
#T3 0.76995124 0.419248295 0.43382791
#T4        NaN         NaN 0.67164249

p.Bon
#      NBD_ES     NBD_JD     ES_JD
#T1 0.6212227 0.03468098 0.1357286
#T2 0.1068355 0.13836822 1.0000000
#T3 1.0000000 1.00000000 1.0000000
#T4       NaN        NaN 1.0000000

#T1 distinguishes between JD and the other two
#T2 distingueshes between NBD and the other two
#this will be used in the building of the classifier
