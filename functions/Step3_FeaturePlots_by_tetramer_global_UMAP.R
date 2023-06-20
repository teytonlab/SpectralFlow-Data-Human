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
#source("/media/rs_DATA/s3/20210111_varner_actri_voucher/heatmap.3_original.R")
source("/Users/smsharma/Documents/OneDrive - Scripps Research/Laboratory Work/Projects T1D/Coding-CCBB-Scripts/heatmap.3_original.R")
####################################end definitions section######################################################


features<-c("HLADR","CXCR4","CD28","CXCR3","CXCR5","CD127","CD69","TIM3","CD45RO","PD1","CTLA4","CD25","ICOS","LAG3","CD62L")
newfeatures<-c("CD62L","CD25","CD45RO","CD28","CD69","LAG3","TIM3","HLADR","CXCR4","CXCR3","CXCR5","CD127","PD1","CTLA4","ICOS")
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
clongcohorts<-rep(ccohorts,c(5,5,4))

colors3<-c("grey20","red","pink")
names(colors3)<-ccohorts

colors5<-c("grey70","deepskyblue4","deeppink3","goldenrod1","green")

colors7<-c("grey70","darkolivegreen","brown","purple","magenta","blue","cyan","grey20","orange","red")
clcolors<-makeTransparent(colors7,alpha=0.08) #this covers A, Aa under one color

names(clcolors)<-c("A","B","C","D","E","F","G","H","I","J")
transtable.list<-list(
"NBD"=c("A","B","C","D","E","F","G","H"),
"ES"=c("A","B","C","D","E","F","G","H"),
"JD"=c("A","B","C","D","E","F","G","H")
)


L<-data.frame()
for (t in ccohorts) {
   load(paste0(t,"_cnormCD4_W_classlabels_pointcolors.RData"))
   L<-rbind(L,W)
}
save(L,file="L.RData")
names(L)
# [1] "HLADR"         "CXCR4"         "CD28"          "CXCR3"        
# [5] "CXCR5"         "CD127"         "CD69"          "TIM3"         
# [9] "CD45RO"        "PD1"           "CTLA4"         "CD25"         
#[13] "ICOS"          "LAG3"          "CD62L"         "tetraPos"     
#[17] "tetraID"       "subjectID"     "cohortID"      "qTetramer"    
#[21] "qtetraPos"     "quantTetramer" "UMAP 1"        "UMAP 2"       
#[25] "classes"       "classlabels"   "pointcolors"  

dim(L)
#[1] 604280     26  



 
trcolors<-makeTransparent(colors7,alpha=0.8)
names(trcolors)<-c("A","B","C","D","E","F","G","H","I","J")

nf<-length(features) #15
nt<-length(ccohorts) #4

#do feature plots
#this takes a while, but it works
mycolrange<-colorRampPalette(c("white","blue4"))
col100<-mycolrange(100) #100 colors from white to blue4
pdf("UMAPs_w_features.pdf",width=4*nf+4,height=4*nt+4) #mice in columns
   par(mar=c(4, 1, 2, 1), mfrow=c(nt+1,nf+1), oma = c(4, 4, 4, 2),pty="s")
   plot.new()
   for (f in features) {
      plot.new()
      text(.5,.5,label=f,cex=4) #make top legend
   }
   for (t in ccohorts) {
      cat(t,"\n")
      load(paste0(t,"_cnormCD4_W_classlabels_pointcolors.RData"))
      rge<-range(W[,23:24])
      plot.new()
      text(.5,.5,label=t,cex=4) #make lhs legend
      for (f in features) {
         zall<-L[,f]
         z99<-quantile(zall,.99)
         z<-pmin(1,W[,f]/z99)
         oz<-order(z)
         zcol<-col100[pmax(1,round(z*100))]
         plot(W[oz,23:24],pch=".",cex=0.1,xlim=rge,ylim=rge,col=zcol[oz],xlab="UMAP 1",ylab="UMAP 2")
      }
   }
dev.off()


#now heatmaps
#for heatmap normalization, find maximum of protein expression in any cluster, any cohort
load("CD4_signatures.RData") #cl_signatures
cl_max<-rep(0,nf)
for (t in ccohorts) cl_max<-pmax(cl_max,apply(cl_signatures[[t]],1,max))
#cl_max<-pmax(cl_max,2)
   
#add heatmap of expression signatures
pdf("CD4_heatmaps.pdf")
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


#Summary object contains mean expression, etc., per marker, per mouse, per cohort, per cell type, also number of cells, as well as number of pos cells, also mgdl per mouse/cohort.
classes<-c("A","B","C","D","E","F","G","H") #main classes
nc<-length(classes)
load("L.RData")



#what is quantile for every subject and tetramer?
L$quantTetramer<-rep(NA,dim(L)[1])
for (i in 1:ncs) {
   for (tetra in tetras) {
      selection<-L$subjectID==csubjects[i] & L$tetraID==tetra
      if (length(selection)>0) L[selection,"quantTetramer"]<-rank(L[selection,"qTetramer"])/sum(selection)
      cat(csubjects[i],tetra,mean(L[L$subjectID==csubjects[i] & L$tetraID==tetra,"qTetramer"]>5),"\n") #one of these is NaN, missing data
   }
}
f<-"ICOS"
plot(ecdf(L[L$cohortID=="JD",f]),pch=NA)
for (q in (995:999)/1000) lines(ecdf(L[L$cohortID=="JD" & L$quantTetramer>q,f]),col="pink",pch=NA)

save(L,file="L.RData")


#now I need a count summary by cluster and cohort
countSummary<-matrix(0,nrow=nc,ncol=ncs) #ncs is 15
rownames(countSummary)<-classes
colnames(countSummary)<-csubjects
cmass<-rep(0,ncs)
names(cmass)<-csubjects
for (s in csubjects) {
   cmass[s]<-sum(L$subjectID==s)
   for (cl in classes) {
      countSummary[cl,s]<-sum(L$subjectID==s & L$classlabels==cl,na.rm=TRUE)
   }
}
fractionSummary<-t(t(countSummary)/colSums(countSummary))

cohortSummary<-matrix(0,nrow=nc,ncol=3)
rownames(cohortSummary)<-classes
colnames(cohortSummary)<-ccohorts 
for (s in ccohorts) {
   for (cl in classes) {
      cohortSummary[cl,s]<-sum(L$cohortID==s & L$classlabels==cl,na.rm=TRUE)
   }
}
cohortFractionSummary<-t(t(cohortSummary)/colSums(cohortSummary))


#TL<-L[L$qtetraPos,]
TL<-L[L$quantTetramer>0.997,]
save(TL,file="TL.RData")
#now just red cells per cluster
redCountSummary<-matrix(0,nrow=nc,ncol=ncs) #ncs is 15
rownames(redCountSummary)<-classes
colnames(redCountSummary)<-csubjects
for (s in csubjects) {
   for (cl in classes) {
      redCountSummary[cl,s]<-sum(TL$subjectID==s & TL$classlabels==cl,na.rm=TRUE)
   }
}
fractionRedCountSummary<-t(t(redCountSummary)/colSums(redCountSummary))

cohortRedCountSummary<-matrix(0,nrow=nc,ncol=3)
rownames(cohortRedCountSummary)<-classes
colnames(cohortRedCountSummary)<-ccohorts
for (s in ccohorts) {
   for (cl in classes) {
      cohortRedCountSummary[cl,s]<-sum(TL$cohortID==s & TL$classlabels==cl,na.rm=TRUE)
   }
}
cohortFractionRedCountSummary<-t(t(cohortRedCountSummary)/colSums(cohortRedCountSummary))





pdf("cluster_fractions.pdf",width=12,height=4)
dx<-0:2
names(dx)<-ccohorts
barplot(t(cohortFractionSummary),beside=TRUE,border=NA,col=colors3,ylab="population fraction",xlab="cluster",ylim=c(0,.7))
legend("topright",legend=ccohorts,pch=15,col=colors3)
for (coh in ccohorts) {
   fr<-fractionSummary[,clongcohorts==coh]
   for (i in 1:nc) {
      x<-(i-1)*4+1+dx[coh]+1/2+rnorm(dim(fr)[2],sd=.1)
      y<-fr[i,]
      points(x,y,cex=0.5,col="grey")
   }
}
box()
dev.off()


pdf("fractions_by_cluster_all_clean_subjects.pdf",width=10,height=5)
par(mfrow=c(1,2),pty="s",mar=c(8, 1, 2, 1))
barplot(fractionSummary,col=colors7,las=2,border=NA,main=expression(all~cells),ylab="population fraction")
barplot(fractionRedCountSummary,col=colors7,las=2,border=NA,main=expression(tetramer^"+"~cells),ylab="population fraction")
dev.off()


pdf("fractions_by_cluster_all_clean_cohorts.pdf",width=10,height=5)
par(mfrow=c(1,2),pty="s",mar=c(8, 1, 2, 1))
barplot(cohortFractionSummary,col=colors7,las=2,border=NA,main=expression(all~cells),ylab="population fraction")
barplot(cohortFractionRedCountSummary,col=colors7,las=2,border=NA,main=expression(tetramer^"+"~cells),ylab="population fraction")
dev.off()

#now just tetras
library(uwot)
library(igraph)
library(leiden)

k<-15
xx<-data.matrix(TL[,1:k])
set.seed(2357123456)
system.time( W<-uwot::umap(X=xx,n_neighbors=k,nn_method="annoy",ret_extra=c("nn","fgraph")) )
save(W,file="tetraW.RData")
U<-W$embedding
colnames(U)<-c("UMAP 1","UMAP 2") # g for global UMAP
rownames(U)<-rownames(TL)


rge<-range(c(U))
plot(U,pch=19,cex=.4,xlim=rge,ylim=rge,col=makeTransparent(TL$pointcolors,alpha=0.9))
plot(U,pch=19,cex=.4,xlim=rge,ylim=rge,col=colors3[TL$cohortID])

gnn<-W$nn[[1]][[1]]

el<-gnn[,c(1,2)] #start edgelist
for (j in 3:k) el<-rbind(el,gnn[,c(1,j)])
swap<-el[,1]>el[,2]
el[swap,]<-el[swap,2:1] #this sorts edge end points so that smaller index goes first
dup<-duplicated(el) #we will remove these duplicated edges when we build the graph
elist<-matrix(as.character(el),nrow=dim(el)[1]) #convert node names to character
g0<-graph_from_edgelist(elist[!dup,], directed = FALSE)
save(g0, file="tetras_g0.RData")
vnames<-as.integer(V(g0)$name)
o<-order(vnames)

#in Constant Potts Model, gamma has a direct interpretation: internal density (number of internal edges/choose(nc,2)) of communities is >= gamma
#see Traag, V. A., Van Dooren, P., & Nesterov, Y. (2011). Narrow scope for resolution-limit-free community detection. Physical Review E, 84(1), 016114.
#this is also resolution-limit-free! i.e., able to detect smaller clusters in big networks
system.time( partition <- leiden(g0,partition_type="CPMVertexPartition",resolution_parameter=0.0044) )
tetraclasses<-partition[o]
print(table(tetraclasses))
#  1   2   3   4 
#901 432 280 140 


plot(U,pch=19,cex=.4,xlim=rge,ylim=rge,col=colors7[tetraclasses])

V<-data.frame(U,TL,tetraclasses,check.names=FALSE)
save(V,file="tetras_V.RData")
